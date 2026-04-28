"""
create_groundtruth_csv.py - 创建csv文件,从所有参考文件中提取位置和时间戳等信息
"""

import os
from os.path import join, isfile, dirname, abspath

DATA_SET = 'dataset'  # 选择数据文件夹
datapath = '../data/' # 相对python脚本的路径

# Get absolute path based on script location
SCRIPT_DIR = dirname(abspath(__file__))
datapath = abspath(join(SCRIPT_DIR, datapath, DATA_SET))

# 选择要处理的手机型号，留空则自动识别数据目录下所有机型文件夹名
PHONES = ['mi8','pixel7pro','sm-g988b','sm-s908b']  # 为空时自动识别(所有机型)

# open output file
os.chdir(datapath)
output_dir = abspath(join(SCRIPT_DIR, '..', 'data'))
fout = open(join(output_dir, 'ground_truths_'+ DATA_SET +'.csv'),'w')
fout.write('tripId,UnixTimeMillis,LatitudeDegrees,LongitudeDegrees, Height\n')

# get list of data sets in data path
datasets = sorted(os.listdir(datapath))

# loop through data set folders
for dataset in datasets:
    if isfile(dataset):
        continue
    dataset_path = join(datapath, dataset)
    if len(PHONES) == 0:
        phones = [d for d in os.listdir(dataset_path) if os.path.isdir(join(dataset_path, d))]
    else:
        phones = sorted(PHONES)
        
    for phone in phones:
        folder = join(datapath, dataset, phone)
        if isfile(folder):
            continue
        csv_file = join(folder, 'ground_truth.csv')
        if not isfile(csv_file):
            continue
        # parse ground truth file
        with open(csv_file) as f:
            lines = f.readlines()[1:]
        flag = 0
        for line in lines:
            if len(line) <= 1:
                continue
            d = line.split(',')
            t = float(d[8]) # get time stamp
            if flag == 0:
                print('%20s,%16s' % (dataset, phone))
                flag = 1
            # write results to combined file
            fout.write('%s/%s,%.0f,%s,%s,%s\n' % ((dataset, phone, t, d[2],
                                                      d[3], d[4][:8])))
fout.close()
