"""
run_multi.py - 将原始Android GNSS日志文件转换为RINEX格式(convert_rnx),
并运行修改后的C版RTKLIB二进制文件进行定位解算(run_rtklib)
"""

import numpy as np
import os
from os.path import join, isdir, isfile, abspath, dirname
from glob import glob
from multiprocessing import Pool
import subprocess
import gnsslogger_to_rnx as rnx
from time import time
import sys
import io

# 设置输入数据集的位置，以及解算结果文件
DATA_SET = 'test'       # 选择数据文件夹
SOL_TAG = 'spp-test'          # 定位解决方案
datapath = '../data/'    # 相对python脚本的路径

# 设置二进制文件和配置文件, 相对python脚本的路径
# binpath_rtklib  = "../rtklib/rnx2rtkp-df.exe"
binpath_rtklib = "D:/Desktop/RTKLIB-VS/Graduation/RTKLIB-PRO/app/consapp/rnx2rtkp/msc/Debug/rnx2rtkp.exe"   
cfgfile_rtklib = "../config/Hybrid-frequency-spp.conf"

# 设置解算选项
OVERWRITE_RINEX = False   # 是否覆盖已存在的rinex文件
NEED_BASE_FILE = False   # 是否需要基站文件 (SPP设为False, PPK设为True)
ENABLE_RTKLIB = True     # 是否使用RTKLIB生成解算结果
OVERWRITE_SOL = True     # 是否覆盖已存在的解算结果文件

# 选择要处理的手机型号，留空则自动识别数据目录下所有机型文件夹名
PHONES = []  # 为空时自动识别(所有机型)

# 设置观测和导航文件的匹配规则
basefiles = '*0.2*o'                      # 观测文件, 支持多种扩展名
navfiles = ['BRDM*MN.rnx', '*0.2*n']      # 导航文件, 支持多种扩展名

# 将相对路径改为绝对路径
SCRIPT_DIR = dirname(abspath(__file__))     # 获取当前脚本所在目录
datadir = abspath(join(SCRIPT_DIR, datapath, DATA_SET))
# binpath_rtklib = abspath(join(SCRIPT_DIR, binpath_rtklib))
cfgfile_rtklib = abspath(join(SCRIPT_DIR, cfgfile_rtklib))

# 解决Windows控制台中文编码问题
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', line_buffering=True)
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', line_buffering=True)

# rinex转换的输入参数结构
class args:
    def __init__(self):
        # rinex转换相关输入参数
        self.slip_mask = 0 # 后续会被覆盖
        self.fix_bias = True
        self.timeadj = 1e-7
        self.pseudorange_bias = 0
        self.filter_mode = 'sync'
        # rinex文件可选头信息
        self.marker_name = ''
        self.observer = ''
        self.agency = ''
        self.receiver_number = ''
        self.receiver_type = ''
        self.receiver_version = ''
        self.antenna_number = ''
        self.antenna_type = ''

# 单个（单线程）rinex文件转换函数
def convert_rnx(args_tuple):
    folder, rawFile, obsFile, slipMask = args_tuple
    argsIn = args()
    argsIn.input_log = rawFile
    argsIn.output = obsFile
    argsIn.slip_mask = slipMask
    rnx.convert2rnx(argsIn)
    return True

# 单个（单线程）RTKLIB解算函数
def run_rtklib(args_tuple):
    binpath_rtklib, cfgfile_rtklib, folder, obsfile, basefile, navfile, solfile = args_tuple
    # 构建命令，只传存在的文件，空文件用 None 跳过
    rtkcmd = [binpath_rtklib, '-k', cfgfile_rtklib, '-o', solfile, obsfile]
    if basefile:
        rtkcmd.append(basefile)
    rtkcmd.append(navfile)
    subprocess.run(rtkcmd, cwd=folder, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return True

####### 主程序入口 ##########################

def main():
    datasets = np.sort(os.listdir(datadir))

    print(f'找到 {len(datasets)} 个数据集: {list(datasets)}')

    rinexIn = []     # 存储需要进行rinex转换的传参列表
    rtklibIn = []    # 存储需要进行RTKLIB解算的传参列表
    for dataset in datasets:

        # 自动识别需要处理的机型
        if len(PHONES) == 0:    # 处理全部机型，根据数据集目录下的文件夹名自动识别
            dataset_path = join(datadir, dataset)
            phones = [d for d in os.listdir(dataset_path) if isdir(join(dataset_path, d))]
        else:                   
            phones = PHONES     # 处理指定机型，根据PHONES列表中的机型名自动识别

        print(f'扫描 {dataset}: {len(phones)} 个手机')

        # 遍历每个手机文件夹
        for phone in phones:
            folder = join(datadir, dataset, phone)      # dataset的根目录下的手机文件夹
            if not isdir(folder):                       # 如果没有此手机的文件夹则跳过
                continue
            rawFile = join(folder, 'supplemental', 'gnss_log.txt')
            obsFile = join(folder, 'supplemental', 'gnss_log.obs')

            rinex = False
            # 检查是否需要进行rinex转换（如果需要覆盖已存在的rinex文件或rinex文件不存在，则需要）
            if OVERWRITE_RINEX or not isfile(obsFile):
                slipMask = 0
                rinexIn.append((folder, rawFile, obsFile, slipMask))
                print('convert_rnx: ', rawFile, '->', obsFile)
                rinex = True

            #--------------------下面进行RTKLIB解算的参数构建--------------------------#

            dataset_root = join(datadir, dataset)                   # dataset的根目录
            base_candidates = glob(join(dataset_root, basefiles))   # 查找基站文件
            nav_candidates = []                                     # 查找导航文件
            for p in navfiles:
                nav_candidates += glob(join(dataset_root, p))

            # 检查文件是否存在（分为是否需要基站文件）
            if NEED_BASE_FILE:
                # PPK模式：必须同时有基站文件和导航文件
                if not base_candidates or not nav_candidates:
                    print(folder, '  错误: 缺少base文件或nav文件')
                    print('   查找目录:', dataset_root)
                    continue
                baseFile = base_candidates[0]
                navFile  = nav_candidates[0]
            else:
                # SPP模式：只需要导航文件，基站文件可为空
                baseFile = ''
                if not nav_candidates:
                    print(folder, '  错误: 缺少nav文件')
                    print('   查找目录:', dataset_root)
                    continue
                navFile = nav_candidates[0]

            solFile = obsFile[:-4] + '_' + SOL_TAG + '.pos'  # 解算结果文件标签

            # 检查是否需要进行RTKLIB解算（用户要求、解算结果文件不存在、需要覆盖已存在的解算结果文件，则需要）
            if ENABLE_RTKLIB and (OVERWRITE_SOL == True or
                    len(glob(solFile)) == 0 or rinex == True):
                print('Run_rtklib: ', join(dataset, phone))
                rtklibIn.append((binpath_rtklib, cfgfile_rtklib,
                                 folder, obsFile, baseFile, navFile, solFile))

    if len(rinexIn) > 0:
        print(f'\n开始并行转换rinex文件 (共{len(rinexIn)}个)...')
        sys.stdout.flush()
        total = len(rinexIn)
        completed = 0
        with Pool() as pool:
            for _ in pool.imap_unordered(convert_rnx, rinexIn, chunksize=1):
                completed += 1
                print(f'\r转换进度: {completed}/{total}', end='', flush=True)
        print('\nrinex转换全部完成!')
        sys.stdout.flush()

    if len(rtklibIn) > 0:
        print(f'\n开始并行计算RTKLIB解算结果 (共{len(rtklibIn)}个)...')
        sys.stdout.flush()
        total = len(rtklibIn)
        completed = 0
        # 单线程模式
        for item in rtklibIn:
            run_rtklib(item)
            completed += 1
            print(f'\r完成: {completed}/{total}', end='', flush=True)
        print('\nRTKLIB解算全部完成!')
        sys.stdout.flush()

if __name__ == '__main__':
    t0 = time()
    main()
    print('Runtime=%.1f' % (time() - t0))