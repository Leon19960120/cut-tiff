# -*- coding: utf-8 -*-
import arcpy
import os
from arcpy.sa import ExtractByMask
# 允许覆盖现有的输出文件
arcpy.env.overwriteOutput = True
# 检查空间分析扩展是否可用
arcpy.CheckOutExtension("Spatial")
# 设置输入的 Shapefile 和 TIFF 文件
input_shp = "D:\\UAV\\0410-jm-LJL\\轻度感染shp\\Export_Output.shp"
input_tiff = "D:\\UAV\\0410test-cut\\MWS-rhw_index_temperature [°c]273.tif"

# 设置输出目录
output_dir = "D:\\UAV\\0410test-cut\\目视轻度感染"
arcpy.env.workspace = "D:\\UAV\\0410test-cut\\目视轻度感染"
# 获取栅格的范围
raster_extent = arcpy.Describe(input_tiff).extent

# 读取 Shapefile
with arcpy.da.SearchCursor(input_shp, ["FID", "SHAPE@"]) as cursor:
    for row in cursor:
        fid = row[0]
        geometry = row[1]  # 获取当前要素的几何信息

        # 检查几何对象是否与栅格范围相交
        if geometry.extent.overlaps(raster_extent) or raster_extent.contains(geometry.extent):
            # 创建掩膜操作的输出文件名
            output_tiff = os.path.join(output_dir, "output_{}.tif".format(fid))

            # 执行掩膜操作
            outExtractByMask = ExtractByMask(input_tiff, geometry)

            # 保存结果
            outExtractByMask.save(output_tiff)

            print("Processed FID {} and saved to {}".format(fid, output_tiff))
        else:
            print("FID {} is outside the raster extent and will not be processed.".format(fid))

# 释放空间分析扩展
arcpy.CheckInExtension("Spatial")

