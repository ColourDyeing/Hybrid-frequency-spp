""" 
get_base_data.py - retrieve base observation and navigation data for the
    2023 GSDC competition 
Modified to use IGS/WHU data centers instead of CDDIS
"""

import os
from os.path import join
from datetime import datetime
import numpy as np
import requests
import gzip
from glob import glob
import subprocess
import multiprocessing
from multiprocessing import Pool, cpu_count

# Input parameters
datadir =  'D:\\Desktop\\RTKLIB-VS\\Reference\\GSDC_2023\\data\\train' # relative to python script
# List of CORS stations to use
stas = ['slac', 'vdcy', 'p222']  # Bay Area, LA, backup for Bay Area

# site to retrieve base observation data
obs_url_base = 'https://geodesy.noaa.gov/corsdata/rinex'  

# Make sure you have downloaded this executable before running this code
crx2rnx_bin = "D:\\Desktop\\RTKLIB-VS\\Reference\\GSDC_2023\\rtklib\\crx2rnx.exe"

def download_nav_data_igs_bkg(year, doy, dataset_dir):
    """从IGS BKG数据中心下载导航数据"""
    try:
        fname = f'BRDM00DLR_S_{year}{doy}0000_01D_MN.rnx.gz'
        url = f'https://igs.bkg.bund.de/root_ftp/IGS/BRDC/{year}/{doy}/{fname}'
        
        print(f'尝试从IGS BKG下载: {url}')
        response = requests.get(url, timeout=30)
        
        if response.status_code == 200:
            nav_data = gzip.decompress(response.content)
            output_file = join(dataset_dir, fname[:-3])  # 去掉.gz后缀
            with open(output_file, "wb") as f:
                f.write(nav_data)
            print(f'IGS BKG下载成功: {fname}')
            return True
        else:
            print(f'IGS BKG HTTP错误: {response.status_code}')
            return False
            
    except Exception as e:
        print(f'IGS BKG下载失败: {e}')
        return False

def download_nav_data_ign(year, doy, dataset_dir):
    """从IGN数据中心下载导航数据 (使用requests处理FTP)"""
    try:
        fname = f'BRDM00DLR_S_{year}{doy}0000_01D_MN.rnx.gz'
        # 转换为HTTP访问
        url = f'http://igs.ign.fr/pub/igs/data/{year}/{doy}/{fname}'
        
        print(f'尝试从IGN下载: {url}')
        response = requests.get(url, timeout=30)
        
        if response.status_code == 200:
            nav_data = gzip.decompress(response.content)
            output_file = join(dataset_dir, fname[:-3])  # 去掉.gz后缀
            with open(output_file, "wb") as f:
                f.write(nav_data)
            print(f'IGN下载成功: {fname}')
            return True
        else:
            print(f'IGN HTTP错误: {response.status_code}')
            return False
            
    except Exception as e:
        print(f'IGN下载失败: {e}')
        return False

def download_nav_data_fallback(year, doy, year2, dataset_dir):
    """备用下载方案 - 尝试旧格式的BRDC文件"""
    try:
        # 尝试传统的brdc格式
        fname = f'brdc{doy}0.{year2}n.gz'
        
        # 多个备用源
        backup_urls = [
            f'https://igs.bkg.bund.de/root_ftp/IGS/BRDC/{year}/{doy}/{fname}',
            f'http://igs.ign.fr/pub/igs/data/{year}/{doy}/{fname}',
        ]
        
        for url in backup_urls:
            try:
                print(f'尝试备用源: {url}')
                response = requests.get(url, timeout=30)
                
                if response.status_code == 200:
                    nav_data = gzip.decompress(response.content)
                    output_file = join(dataset_dir, fname[:-3])  # 去掉.gz后缀
                    with open(output_file, "wb") as f:
                        f.write(nav_data)
                    print(f'备用源下载成功: {fname}')
                    return True
                    
            except Exception as e:
                print(f'备用源 {url} 失败: {e}')
                continue
                
        return False
        
    except Exception as e:
        print(f'备用下载方案失败: {e}')
        return False

def process_navigation_data(args):
    """处理单个数据集的导航数据下载"""
    dataset, year, doy, year2 = args
    dataset_path = join(datadir, dataset)
    
    # 检查导航文件是否已存在
    nav_files = glob(join(dataset_path,'*.rnx')) + glob(join(dataset_path,'*.*n'))
    if len(nav_files) > 0:
        print(f'导航文件已存在: {nav_files[0]}')
        return True
    
    print(f'开始下载导航数据: {dataset} ({year}/{doy})')
    
    # 按优先级尝试不同的数据源
    success = False
    
    # 1. 尝试IGS BKG数据中心
    if not success:
        success = download_nav_data_igs_bkg(year, doy, dataset_path)
    
    # 2. 尝试IGN数据中心
    if not success:
        success = download_nav_data_ign(year, doy, dataset_path)
    
    # 3. 最后尝试备用方案
    if not success:
        success = download_nav_data_fallback(year, doy, year2, dataset_path)
    
    if not success:
        print(f'所有数据源都失败，导航数据下载失败: {dataset}')
        return False
    
    return True

def process_observation_data(dataset, year, doy, year2):
    """处理观测数据下载和转换"""
    dataset_path = join(datadir, dataset)
    
    # 处理观测数据部分 (保持原有逻辑)
    if len(glob(join(dataset_path,'*.*o'))) == 0:
        # get obs data
        i = 1 if '-LAX-' in dataset.upper() else 0  # use different base for LA
        fname = stas[i] + doy + '0.' + year2 + 'd.gz'
        url = '/'.join([obs_url_base, year, doy, stas[i], fname])
        try:
            obs = gzip.decompress(requests.get(url).content) # get obs and decompress
            # write obs data
            open(join(dataset_path, fname[:-3]), "wb").write(obs)
        except:
            # try backup CORS station
            print('Try backup CORS:', dataset)
            i += 2
            fname = stas[i] + doy + '0.' + year2 + 'd.gz'
            url = '/'.join([obs_url_base, year, doy, stas[i], fname])
            try:
                obs = gzip.decompress(requests.get(url).content) # get obs and decompress
                # write obs data
                open(join(dataset_path, fname[:-3]), "wb").write(obs)
            except:
                print('Fail obs: %s' % dataset)
            
        # convert compact rinex to rinex
        crx_files = glob(join(dataset_path,'*.*d'))
        if len(crx_files) > 0:
            crx_file = crx_files[0]
            # 预期的输出文件名（.d改为.o）
            expected_output = crx_file[:-1] + 'o'
            
            ret_code = subprocess.call([crx2rnx_bin, '-f', crx_file])
            
            # 检查转换是否成功且输出文件存在
            if ret_code == 0 and os.path.exists(expected_output):
                os.remove(crx_file)
                print(f'转换成功，已删除: {crx_file}')
            else:
                print(f'转换失败或输出文件不存在，保留: {crx_file}')

# 主处理逻辑
if __name__ == '__main__':
    # Loop through data sets in the data directory
    os.chdir(datadir)
    datasets = []
    nav_tasks = []
    
    for dataset in np.sort(os.listdir()):
        if not os.path.isdir(join(dataset)):
            continue
        print(f"\n处理数据集: {dataset}")
        ymd = dataset.split('-')
        doy = datetime(int(ymd[0]), int(ymd[1]), int(ymd[2])).timetuple().tm_yday # get day of year
        doy = str(doy).zfill(3)
        year = ymd[0]
        year2 = ymd[0][2:4]  # 两位年份
        
        datasets.append((dataset, year, doy, year2))
        nav_tasks.append((dataset, year, doy, year2))
        
        # 处理观测数据（保持串行处理）
        process_observation_data(dataset, year, doy, year2)
    
    # 使用多进程处理导航数据下载
    if nav_tasks:
        print(f"\n开始多进程下载导航数据，使用 {cpu_count()} 个CPU核心")
        with Pool(processes=cpu_count()) as pool:
            results = pool.map(process_navigation_data, nav_tasks)
        
        # 统计结果
        success_count = sum(results)
        total_count = len(results)
        print(f"\n导航数据下载完成: {success_count}/{total_count} 成功")
    else:
        print("没有找到需要处理的数据集")