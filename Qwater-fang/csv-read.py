import pandas as pd
import os
import glob

# 定义目录和输出文件路径
directory = 'D:\\UAV\\0410test-cut\\目视重度感染nodata'
output_csv_path = 'D:\\UAV\\0410test-cut\\目视重度感染nodata\\mintemp.csv'  # 确保这是一个CSV文件路径,记得改路径，命名就自己起一个好了

# 获取目录下所有CSV文件的路径
csv_files = glob.glob(os.path.join(directory, '*.csv'))

# 准备一个列表来收集数据
data_frames = []

for file_path in csv_files:
    # 读取CSV文件
    df = pd.read_csv(file_path, encoding='ISO-8859-1')
    # 检查DataFrame是否至少有三行
    if df.shape[0] >= 20:
        # 提取倒数第三行的数据
        # 提取特定列的数据，例如提取名为'Parameter'和'Value'的列
        specific_data = df[['Parameter', 'Value']].iloc[2]  # 如果你需要特定列
        specific_data = specific_data.to_frame().T  # 转换为DataFrame

        # 添加文件名标识符（可选）
        file_name = os.path.basename(file_path)
        specific_data['file_name'] = file_name
        # 将DataFrame添加到列表中
        data_frames.append(specific_data)
    else:
        print(f"Warning: {file_path} has less than 3 rows and will be skipped.")

# 使用concat合并所有DataFrame
combined_data = pd.concat(data_frames, ignore_index=True)

# 保存汇总后的数据到新的CSV文件
combined_data.to_csv(output_csv_path, index=False)

