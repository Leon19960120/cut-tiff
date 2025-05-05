# -*- coding: utf-8 -*-
'''
/***************************************************************************
QWaterModel - a QGIS plugin
QWaterModel is a simple tool to calculate Evapotranspiration from thermal images.
                              -------------------
        begin                : 2020-01-11
        git sha              : $Format:%H$
        copyright            : (C) 2020 by Florian Ellsäßer
        email                : el-flori@gmx.de

This QGis plugin is based on:

Ellsäßer et al. (2020), Introducing QWaterModel, a QGIS plugin for predicting
evapotranspiration from land surface temperatures,
Environmental Modelling & Software, https://doi.org/10.1016/j.envsoft.2020.104739.

The DATTUTDUT energy-balance model is based on:
Timmermans, W.J., Kustas, W.P., Andreu, A., 2015. Utility of an automated
thermal-based approach for monitoring evapotranspiration.
Acta Geophys. 63, 1571–1608. https://doi.org/10.1515/acgeo-2015-0016.

More information: https://github.com/FloEll/QWaterModel
 ***************************************************************************/
For the generation of this plugin, I used:
Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
'''


# Import Libraries for the EB-Model
import os
import glob
import csv
from osgeo import gdal
import numpy as np
import math
import datetime
import utm
gdal.UseExceptions()

# Define global variables
sb_const = 5.6704 * 10 ** (-8)  # Stefan Bolzmann constant
sw_exo = 1361.5  # exo-atmospheric short wave radiation
air_temp = 299  # 记得修改这里的空气温度，单位k
time_period = 3600  # 小时蒸腾量的意思
UTC_TIME = "2024-04-10T02:50:00"  # 记得改，UTC时间都改成对应的
atm_trans = 0.7  # 大气透射率，就定死吧
atm_emis = 0.8  # 大气发射率也定死啦
sw_irr = 708.5886  #短波辐射，记得改!!!!!
surf_emis = 1 #表面反射率吧好像是

# 定义输入文件夹路径
input_directory = "D:\\UAV\\0410test-cut\\目视重度感染nodata"  # 将此路径替换为你的输入文件夹路径

def get_raster_centroid(dataset):
    """计算并返回栅格的中心点经纬度。"""
    geotransform = dataset.GetGeoTransform()
    x_center = geotransform[0] + geotransform[1] * (dataset.RasterXSize / 2)
    y_center = geotransform[3] + geotransform[5] * (dataset.RasterYSize / 2)
    return x_center, y_center

#def get_lon_lat(x_center, y_center, zone_number, zone_letter):
#    latitude, longitude = utm.to_latlon(x_center, y_center, zone_number, zone_letter)
#    print(latitude,longitude)
#    return latitude, longitude

def get_raster_min_max(dataset):
    """获取栅格的最小值和最大值，忽略NoData值。"""
    band = dataset.GetRasterBand(1)
    no_data_value = band.GetNoDataValue()
    data = band.ReadAsArray()
    # 将NoData值设置为nan以便使用numpy的nan函数
    data = np.where(data == no_data_value, np.nan, data)
    min_val = np.nanmin(data)
    max_val = np.nanmax(data)
    return min_val, max_val


def get_albedo(lst, min_val, max_val):
    """
    Determines surface albedo based on land surface temperatures (lst),
    minimum (tmin) and maximum (tmax) temperatures.
    """
    # 确保最小值和最大值有效，并且不相等
    if min_val is None or max_val is None or min_val == max_val:
        return None  # 如果条件不满足，返回None或某个默认值

    # 避免分母为零的情况，进行albedo的计算
    diff = max_val - min_val
    if diff == 0:
        return None  # 或者返回一个合理的默认值，如0.0

    albedo = abs(0.05 + ((lst - min_val) / (max_val - min_val)) * 0.2)

    if np.isnan(np.sum(albedo)):
        return None  # 或者某个默认值

    albedo[albedo > 1.0] = 0.25
    albedo[albedo < 0.0] = 0.05

    return albedo


def read_lst_img(in_file, na_val=273.15):
    dataset = gdal.Open(in_file, gdal.GA_ReadOnly)
    band = dataset.GetRasterBand(1)
    lst = band.ReadAsArray().astype(float)
    if np.all(np.isnan(lst)):
        # 如果数据全是NaN，则认为没有有效数据
        return None
    lst = np.where(np.isnan(lst), na_val, lst).copy()
    # 使用np.isnan()检测数组中的nan值
    lst[np.isnan(lst)] = na_val
    return lst


def get_sol_elev_ang(utc_str, y_center, x_center):
    """Calculate the solar elevation angle using provided UTC time and geolocation."""

    # 解析UTC时间，获取日序和时间
    utc = datetime.datetime.strptime(utc_str, '%Y-%m-%dT%H:%M:%S')
    doy = float(utc.strftime('%j'))
    daytime = float(utc.hour + utc.minute / 60)

    # 计算天文参数
    g = (360 / 365.25) * (doy + daytime / 24)
    g_radians = math.radians(g)
    declination = (0.396372 - 22.91327 * math.cos(g_radians) + 4.02543 * math.sin(g_radians)
                   - 0.387205 * math.cos(2 * g_radians) + 0.051967 * math.sin(2 * g_radians)
                   - 0.154527 * math.cos(3 * g_radians) + 0.084798 * math.sin(3 * g_radians))

    time_correction = (0.004297 + 0.107029 * math.cos(g_radians) - 1.837877 * math.sin(g_radians)
                       - 0.837378 * math.cos(2 * g_radians) - 2.340475 * math.sin(2 * g_radians))

    SHA = (daytime - 12) * 15 + x_center + time_correction
    SHA = SHA - 360 if SHA > 180 else SHA + 360 if SHA < -180 else SHA

    # 计算太阳高度角
    lat_radians = math.radians(y_center)
    d_radians = math.radians(declination)
    SHA_radians = math.radians(SHA)

    SZA_radians = math.acos(math.sin(lat_radians) * math.sin(d_radians) +
                            math.cos(lat_radians) * math.cos(d_radians) * math.cos(SHA_radians))

    SZA = math.degrees(SZA_radians)
    SEA = 90 - SZA

    return SEA

def get_evap_frac(lst, min_val, max_val):
    '''蒸发分数计算'''
    ef = (max_val - lst) / (max_val - min_val)
    ef = np.where(np.isnan(ef), 0, ef)  # 替换NaN为0
    ef = np.clip(ef, 0, 1)  # 限制ef的值在0到1之间
    return ef

def get_rn(albedo, sw_irr, surf_emis, atm_emis, air_temp, lst, rn_specified=None):
    '''净辐射计算'''
    if rn_specified is None:
        rn = ((1 - albedo) * sw_irr +
              surf_emis * atm_emis * sb_const * (air_temp ** 4) -
              surf_emis * sb_const * (lst ** 4))
    else:
        rn = rn_specified
    return rn

def get_g(rn, lst, min_val, max_val, g_percentage=None):
    '''地面热通量计算'''
    if g_percentage is None:
        g = rn * (0.05 + ((lst - min_val) / (max_val - min_val)) * 0.4)
    else:
        g = rn * (g_percentage / 100)
    g = np.clip(g, 0.05, 0.45)  # 限制g的值在0.05到0.45之间，这里的范围可以根据实际情况调整
    return g

def get_h_le(rn, g, ef):
    '''感热通量、潜热通量计算'''
    le = (rn - g) * ef
    h = (rn - g) - le
    return h, le

def get_water(le, air_temp, time_period, water_specified=None):
    '''蒸腾计算'''
    if water_specified is None:
        water = (le * time_period / 1000000) / (2.501 - 0.002361 * (air_temp - 273.15))
    else:
        water = water_specified
    return water


def process_tiff(file_path):
    """处理单个TIFF文件，获取中心点经纬度和最大最小值。"""
    dataset = gdal.Open(file_path)
    if dataset is None:
        print(f"无法打开文件：{file_path}")
        return None
    #centroid = get_raster_centroid(dataset)
    lon, lat = get_raster_centroid(dataset)  # 获取中心点坐标
    #lon, lat = get_lon_lat(*centroid, 50, 'N')  # 将中心点坐标转换为经纬度
    lst = read_lst_img(file_path)
    if lst is None:
        print(f"栅格数据全为NoData值：{file_path}。跳过此文件。")
        return None

    min_val, max_val = get_raster_min_max(dataset)
    # 在这里添加检查以确保存在有效数据
    if min_val is None or max_val is None or np.isnan(min_val) or np.isnan(max_val):
        print(f"Raster {file_path} contains no valid data. Skipping...")
        return None  # 无有效数据，返回None

    albedo = get_albedo(lst, min_val, max_val)
    if albedo is None:
        print(f"反照率计算失败：{file_path}。跳过此文件。")
        return None

    ef = get_evap_frac(lst, min_val, max_val)
    rn = get_rn(albedo, sw_irr, surf_emis, atm_emis, air_temp, lst)
    g = get_g(rn, lst, min_val, max_val)
    h,le = get_h_le(rn, g, ef)
    water = get_water(le, air_temp, time_period)
    albedo = get_albedo(lst, min_val, max_val)
    stats = {
        'file_path': file_path,
        'centroid_lon_lat': (lon, lat),
        'min_val': min_val,
        'max_val': max_val,
    }
    for var_name in [ 'albedo','rn', 'ef', 'g', 'le', 'water']:
        var = locals()[var_name]  # 从本地变量中获取变量值
        stats[var_name + '_mean'] = np.mean(var[~np.isnan(var)])
        stats[var_name + '_min'] = np.min(var[~np.isnan(var)])
        stats[var_name + '_max'] = np.max(var[~np.isnan(var)])

    return stats

def process_all_tiffs(directory):
    """处理指定目录下的所有TIFF文件，并为每个文件生成一个CSV文件。"""
    tif_files = glob.glob(os.path.join(directory, '*.tif'))
    for file_path in tif_files:
        raster_info = process_tiff(file_path)

        # 只有在raster_info非None时才继续
        if raster_info is not None:
            # 继续生成CSV文件路径并写入数据到CSV
            base_name = os.path.basename(file_path)
            csv_file_name = f"{os.path.splitext(base_name)[0]}_stats.csv"
            output_csv_path = os.path.join(directory, csv_file_name)
            write_tiff_data_to_csv(output_csv_path, raster_info)
        else:
            # 如果raster_info是None，记录一条消息并跳过这个文件
            print(f"由于数据无效，跳过文件: {file_path} 的CSV写入。")


def write_tiff_data_to_csv(output_csv_path, tiff_data):
    """将处理TIFF文件得到的数据写入到CSV文件中。"""
    with open(output_csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # 写入TIFF文件和地理信息
        writer.writerow(['Parameter', 'Value'])
        writer.writerow(['File Path', tiff_data['file_path']])
        writer.writerow(['Centroid Lon/Lat', tiff_data['centroid_lon_lat']])
        writer.writerow(['Min Val', tiff_data['min_val']])
        writer.writerow(['Max Val', tiff_data['max_val']])

        # 写入每个变量的统计数据
        for key in tiff_data:
            if key.endswith('_mean') or key.endswith('_min') or key.endswith('_max'):
                writer.writerow([key, tiff_data[key]])

if __name__ == "__main__":
    process_all_tiffs(input_directory)
