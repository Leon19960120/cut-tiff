from osgeo import gdal
import numpy as np
import glob
import os
gdal.UseExceptions()

#在分割完了单个树冠后要进行这一步，记得修改路径，不然最大温度值读取不出来
def process_tiff(input_image_path, output_image_path):
    dataset = gdal.Open(input_image_path)
    band = dataset.GetRasterBand(1)
    nodata_value = band.GetNoDataValue()
    # 读取栅格数据为数组
    data_array = band.ReadAsArray().astype(float)
    # 如果没有找到有效的NoData值，则设置一个
    if nodata_value is None:
        nodata_value = -32767  # 你可能需要手动设置这个正确的NoData值。

    # 创建一个掩码，将NoData值标记为True
    nodata_mask = data_array == nodata_value
    # 然后将这些NoData值的位置设置为np.nan
    data_array[nodata_mask] = np.nan

    driver = gdal.GetDriverByName('GTiff')
    out_dataset = driver.Create(output_image_path, dataset.RasterXSize, dataset.RasterYSize, 1, band.DataType)
    out_dataset.SetGeoTransform(dataset.GetGeoTransform())
    out_dataset.SetProjection(dataset.GetProjection())

    out_band = out_dataset.GetRasterBand(1)
    out_band.WriteArray(data_array)
    out_band.SetNoDataValue(np.nan)

    out_band.FlushCache()
    out_dataset = None
    dataset = None

# 定义输入和输出目录
input_directory = "D:\\UAV\\0410test-cut\\目视重度感染"
output_directory = "D:\\UAV\\0410test-cut\\目视重度感染nodata"

# 创建输出目录如果它不存在
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# 获取所有TIF文件的路径
tiff_files = glob.glob(os.path.join(input_directory, '*.tif'))

# 批量处理所有TIF文件
for input_image_path in tiff_files:
    # 创建输出文件的路径
    filename = os.path.basename(input_image_path)
    output_image_path = os.path.join(output_directory, filename)

    # 处理每个TIF文件
    process_tiff(input_image_path, output_image_path)

print("Processing complete.")
