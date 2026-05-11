""" 
create_baseline_csv_from_pos.py - 创建csv文件,使用参考文件中的时间戳,从解算结果文件中提取位置和质量等信息
"""

import os
from os.path import join, isfile, dirname, abspath
import numpy as np
from datetime import datetime

########### Input parameters ###############################

DATA_SET = 'train' # 选择数据文件夹
SOL_TAG = 'SPP-HFUC-RM-IGG3' # 解算文件标签
datapath = '../data/' # 相对python脚本的路径
rovfile = 'gnss_log'
hdrlen = 25    # 25 for RTKLIB, 1 for RTKLIB-py，表示跳过前多少行

outThresh = 100   # max horizontal accuracy estimate（本次不使用这个参数，后续可以根据需要添加）
# Select all phones to process
# PHONES = []
PHONES = ['mi8','pixel6pro','pixel7pro','s20ultra','s21ultra','s22ultra']  # 为空时自动识别(所有机型)

############################################################

# Get absolute path based on script location
SCRIPT_DIR = dirname(abspath(__file__)) # 获取脚本所在目录
datapath = abspath(join(SCRIPT_DIR, datapath)) # 获取数据所在目录

GPS_TO_UTC = 315964782  # second

def create_csv(datapath, DATA_SET, SOL_TAG):
    # get timestamps from existing baseline file
    os.chdir(datapath)
    baseline_file = 'ground_truths_' + DATA_SET + '.csv'

    # read data from baseline file
    base_txt = np.genfromtxt(baseline_file, delimiter=',',invalid_raise=False, 
                             skip_header=1, dtype=str)
    msecs_base = base_txt[:,1].astype(np.int64)
    phones_base = base_txt[:,0]
    pos_base = base_txt[:,2:4].astype(float) # baseline positions
    
    # open output file
    output_filename = 'locations_' + SOL_TAG + '_' + DATA_SET + '_' + datetime.now().strftime("%m_%d") + '.csv'
    fout =open(output_filename,'w') # 保存的文件名字
    fout.write('tripId,UnixTimeMillis,LatitudeDegrees,LongitudeDegrees,Height,Quality,NumSatellites,Sde,Sdu,Sdn,Sdne,Sdeu,Sdun\n')
    
    # get list of data sets in data path
    os.chdir(join(datapath, DATA_SET))
    trips = np.sort(os.listdir())
    
    # loop through data set folders
    ix_b, npts = [], 0
    for trip in trips:
        if isfile(trip):
            continue
        phones = os.listdir(trip)
        # loop through phone folders
        for phone in phones:
            if isinstance(phone, bytearray):
                phone = phone.decode('utf-8')
            # check for valid folder and file
            folder = join(trip, phone)
            if isfile(folder):
                continue
            if PHONES != [] and phone not in PHONES:
                continue
            trip_phone = trip + '/' + phone
            #print(trip_phone)
    
            ix_b = np.where(phones_base == trip_phone)[0]
            sol_path = join(folder, 'supplemental', rovfile + '_' + SOL_TAG + '.pos')
            fields = []
            if isfile(sol_path):
                # parse solution file
                fields = np.genfromtxt(sol_path, invalid_raise=False, skip_header=hdrlen)
            if len(fields) > 1:
                if int(fields[0,1]) > int(fields[-1,1]): # invert if backwards solution
                    fields = fields[::-1]
                pos = fields[:,2:5]
                qs = fields[:,5].astype(int)
                nss = fields[:,6].astype(int)
                acc = fields[:,7:10]
                # Extract covariance: sdne (col 10), sdeu (col 11), sdun (col 12)
                cov = fields[:,10:13] if fields.shape[1] >= 13 else np.zeros((len(fields), 3))
                msecs = (1000 * (fields[:,0] * 7 * 24 * 3600 + fields[:,1])).astype(np.int64)
                msecs += GPS_TO_UTC * 1000
            # if no data, use baseline data
            if not isfile(sol_path) or len(fields) == 0:
                print('Warning: data substitution: ', sol_path)
                msecs = msecs_base[ix_b].copy()
                pos = acc = cov = np.zeros((len(msecs), 3))
                pos[:,:2] = pos_base[ix_b].copy()
                qs = nss = np.zeros(len(msecs))
           
            # Find matching ground truth for each solution timestamp (no interpolation)
            # Instead of interpolating solutions to all ground truth times,
            # we match each solution time to the nearest ground truth time
            llhs = [[] for _ in range(3)]
            stds = [[] for _ in range(3)]
            # Add covariance lists: sdne, sdeu, sdun
            sdne_list, sdeu_list, sdun_list = [], [], []
            qsi_list = []
            nssi_list = []
            output_msecs = []
            
            # For each solution timestamp, find the nearest ground truth timestamp
            for j in range(len(msecs)):
                sol_time = msecs[j]
                # Find the nearest ground truth time for this trip
                if len(ix_b) > 0:
                    # Get ground truth times for this trip
                    trip_gt_times = msecs_base[ix_b]
                    # Find nearest index
                    nearest_idx = np.argmin(np.abs(trip_gt_times - sol_time))
                    nearest_time = trip_gt_times[nearest_idx]
                    # Only use if within a reasonable threshold (e.g., 1 second = 1000ms)
                    if abs(nearest_time - sol_time) < 1000:
                        # Use the ground truth timestamp for output (so it matches during merge)
                        output_msecs.append(nearest_time)
                        for k in range(3):
                            llhs[k].append(pos[j, k])
                            stds[k].append(acc[j, k])
                        # Add covariance: sdne, sdeu, sdun
                        sdne_list.append(cov[j, 0])
                        sdeu_list.append(cov[j, 1])
                        sdun_list.append(cov[j, 2])
                        qsi_list.append(qs[j])
                        nssi_list.append(nss[j])
            
            # Convert to arrays
            if len(output_msecs) > 0:
                output_msecs = np.array(output_msecs)
                llhs = np.array(llhs).T  # shape: (n, 3)
                stds = np.array(stds).T  # shape: (n, 3)
                sdne_arr = np.array(sdne_list)
                sdeu_arr = np.array(sdeu_list)
                sdun_arr = np.array(sdun_list)
                qsi_arr = np.array(qsi_list)
                nssi_arr = np.array(nssi_list)
    
            # write results to combined file
            if len(output_msecs) > 0:
                for i in range(len(output_msecs)):
                    fout.write('%s,%d,%.9f,%.9f,%.4f,%.0f,%.0f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n' % 
                            (trip_phone, output_msecs[i], llhs[i,0], llhs[i,1],
                             llhs[i,2], qsi_arr[i], nssi_arr[i], stds[i,0], stds[i,1], 
                             stds[i,2], sdne_arr[i], sdeu_arr[i], sdun_arr[i]))
                try:
                    npts += len(output_msecs)
                except:
                    pass
    
    fout.close()
    return npts

if __name__ == '__main__':
    create_csv(datapath, DATA_SET, SOL_TAG)
