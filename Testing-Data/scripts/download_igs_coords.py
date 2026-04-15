#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从IGS官网自动下载测站真值坐标
使用IGS每周SINEX文件获取ITRF2020框架下的精确坐标
"""

import os
import sys
import re
import urllib.request
import urllib.error
from datetime import datetime, timedelta
from pathlib import Path

# 香港IGS测站列表
HK_STATIONS = [
    'hkcl', 'hkks', 'hkkt', 'hklm', 'hklt', 'hkmw', 'hknp',
    'hkoh', 'hkpc', 'hkqt', 'hksc', 'hksl', 'hkss', 'hkst',
    'hktk', 'hkws', 'wuhn'
]

# IGS SINEX文件URL模板
# IGS提供多种SINEX格式：
# 1. igs.snx - 最新综合SINEX
# 2. igsWWW0.snx - 按周命名
IGS_SINEX_URLS = [
    'https://cddis.nasa.gov/archive/gnss/snx/igs.snx',
    'https://cddis.nasa.gov/archive/gnss/snx/igs21P03160.snx',  # 2021年第316周
    'https://cddis.nasa.gov/archive/gnss/snx/igs22P03510.snx',  # 2022年
    'https://cddis.nasa.gov/archive/gnss/snx/igs23P03960.snx',  # 2023年
    'https://cddis.nasa.gov/archive/gnss/snx/igs24P04080.snx',  # 2024年
    'https://cddis.nasa.gov/archive/gnss/snx/igs25P04200.snx',  # 2025年
]


def get_latest_igs_snx_url():
    """获取最新的IGS SINEX文件URL"""
    current_year = datetime.now().year
    current_week = datetime.now().isocalendar()[1]

    # 尝试获取最新的URL
    for year in range(current_year, 2019, -1):
        for week in range(52, 0, -1):
            week_str = f"P{year % 100:02d}{week:02d}0"
            url = f'https://cddis.nasa.gov/archive/gnss/snx/igs{week_str}.snx'
            yield url
        week = 52  # 继续当年的其他周


def parse_snx_file(snx_file):
    """解析SINEX文件，提取测站坐标

    SINEX格式说明：
    %COMMENT
    ...
    %MARKER_NAME
    ...
    %SOLUTION/ESTIMATE
    参数说明：
    1: 参数类型 (STA, POS, VEL, etc.)
    2: 站点代码
    3: 描述
    4: 约束标志
    5: 约束值
    6: 约束sigma
    7: 单位

    站点坐标以XYZ格式存储：
    TYPE  CODE  PT  description  con  X(m)  y-sig  unit
    TYPE  CODE  PT  description  con  Y(m)  y-sig  unit
    TYPE  CODE  PT  description  con  Z(m)  z-sig  unit
    """
    stations_coords = {}

    try:
        with open(snx_file, 'r', errors='ignore') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"[错误] 文件不存在: {snx_file}")
        return {}

    in_solution = False
    current_station = None
    current_coords = {}

    for line in lines:
        line = line.strip()

        # 找到解决方案开始标记
        if line.startswith('%SOLUTION/ESTIMATE'):
            in_solution = True
            continue

        # 解决方案结束
        if line.startswith('%ENDSOLUTION'):
            in_solution = False
            continue

        if not in_solution:
            continue

        # 解析数据行
        if line.startswith('TYPE'):
            parts = line.split()
            if len(parts) >= 7:
                param_type = parts[0]
                code = parts[1]
                pt = parts[2]
                desc = parts[3]

                if param_type in ['TYPE', 'STA']:  # STA用于站点名
                    code_upper = code.upper()
                    if code_upper not in stations_coords:
                        stations_coords[code_upper] = {'desc': desc}
                    current_station = code_upper
                    current_coords = stations_coords[code_upper]

                elif param_type == 'POS':
                    code_upper = code.upper()
                    if code_upper not in stations_coords:
                        stations_coords[code_upper] = {'desc': desc}

                    # 检查是否是坐��参数
                    if pt == '   X':
                        try:
                            stations_coords[code_upper]['X'] = float(parts[5])
                        except (ValueError, IndexError):
                            pass
                    elif pt == '   Y':
                        try:
                            stations_coords[code_upper]['Y'] = float(parts[5])
                        except (ValueError, IndexError):
                            pass
                    elif pt == '   Z':
                        try:
                            stations_coords[code_upper]['Z'] = float(parts[5])
                        except (ValueError, IndexError):
                            pass

    return stations_coords


def llh_to_xyz(lat, lon, hgt, a=6378137.0, f=1/298.257223563):
    """大地坐标转地心地固坐标"""
    import math
    e2 = 2 * f - f * f
    sin_lat = math.sin(math.radians(lat))
    cos_lat = math.cos(math.radians(lat))
    N = a / math.sqrt(1 - e2 * sin_lat ** 2)

    x = (N + hgt) * cos_lat * math.cos(math.radians(lon))
    y = (N + hgt) * cos_lat * math.sin(math.radians(lon))
    z = (N * (1 - e2) + hgt) * sin_lat
    return x, y, z


def xyz_to_llh(x, y, z, a=6378137.0, f=1/298.257223563):
    """地心地固坐标转大地坐标"""
    import math
    e2 = 2 * f - f * f

    lon = math.atan2(y, x)
    p = math.sqrt(x**2 + y**2)

    # 迭代计算纬度
    lat = math.atan2(z, p * (1 - e2))  # 初值
    for _ in range(10):
        sin_lat = math.sin(lat)
        N = a / math.sqrt(1 - e2 * sin_lat**2)
        lat = math.atan2(z + e2 * N * sin_lat, p)

    sin_lat = math.sin(lat)
    N = a / math.sqrt(1 - e2 * sin_lat**2)
    hgt = p / math.cos(lat) - N

    return math.degrees(lat), math.degrees(lon), hgt


def download_snx(url, output_file):
    """下载SINEX文件"""
    print(f"正在下载: {url}")
    try:
        # 设置请求头，模拟浏览器访问
        request = urllib.request.Request(
            url,
            headers={
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
            }
        )

        with urllib.request.urlopen(request, timeout=60) as response:
            content = response.read()
            with open(output_file, 'wb') as f:
                f.write(content)
            print(f"  下载完成: {output_file} ({len(content)} bytes)")
            return True

    except urllib.error.HTTPError as e:
        print(f"  [HTTP错误] {e.code}: {e.reason}")
        return False
    except urllib.error.URLError as e:
        print(f"  [网络错误] {e.reason}")
        return False
    except Exception as e:
        print(f"  [错误] {e}")
        return False


def main():
    import argparse

    parser = argparse.ArgumentParser(description='从IGS下载测站真值坐标')
    parser.add_argument('--output', '-o', default=r'D:\Desktop\RTKLIB-VS\Graduation\Testing-Data\data\igs\true_crd.true_crd',
                        help='输出真值文件路径')
    parser.add_argument('--snx-dir', default=r'D:\Desktop\RTKLIB-VS\Graduation\Testing-Data\data\igs',
                        help='SINEX文件保存目录')
    parser.add_argument('--stations', '-s', nargs='+', default=HK_STATIONS,
                        help='测站列表')
    args = parser.parse_args()

    print("="*60)
    print(" IGS测站真值坐标自动下载程序")
    print("="*60)

    # 确保目录存在
    os.makedirs(args.snx_dir, exist_ok=True)

    # 尝试下载SINEX文件
    snx_file = os.path.join(args.snx_dir, 'igs.snx')
    success = False

    # 方法1：尝试下载igs.snx
    print("\n[方法1] 尝试下载最新IGS综合SINEX文件...")
    if download_snx(IGS_SINEX_URLS[0], snx_file):
        success = True
    else:
        # 方法2：尝试按周下载
        print("\n[方法2] 尝试按周下载SINEX文件...")
        for url in IGS_SINEX_URLS[1:]:
            if download_snx(url, snx_file):
                success = True
                break

        # 方法3：生成周URL并尝试
        if not success:
            print("\n[方法3] 生成周URL...")
            for url in get_latest_igs_snx_url():
                if download_snx(url, snx_file):
                    success = True
                    break

    # 解析SINEX文件
    if success and os.path.exists(snx_file):
        print("\n正在解析SINEX文件...")
        stations = parse_snx_file(snx_file)
        print(f"共解析到 {len(stations)} 个测站")

        # 提取目标测站坐标
        output_lines = []
        output_lines.append("% IGS测站真值坐标 (ITRF2020)")
        output_lines.append("% 格式: 测站名  X(m)  Y(m)  Z(m)  纬度(deg)  经度(deg)  高程(m)")
        output_lines.append("")

        found_stations = []
        for station in args.stations:
            station_upper = station.upper()
            if station_upper in stations:
                coord = stations[station_upper]
                if 'X' in coord and 'Y' in coord and 'Z' in coord:
                    x, y, z = coord['X'], coord['Y'], coord['Z']
                    lat, lon, hgt = xyz_to_llh(x, y, z)
                    desc = coord.get('desc', '')

                    output_lines.append(
                        f"{station_upper:6} {x:>12.3f}  {y:>12.3f}  {z:>12.3f}  "
                        f"{lat:>12.6f}  {lon:>12.6f}  {hgt:>8.3f}"
                    )
                    found_stations.append(station_upper)
                    print(f"  {station_upper}: lat={lat:.6f} lon={lon:.6f} hgt={hgt:.3f}m")
                else:
                    print(f"  [缺失坐标] {station_upper}")
            else:
                print(f"  [未找到] {station_upper}")

        # 保存结果
        with open(args.output, 'w', encoding='utf-8') as f:
            f.write('\n'.join(output_lines))

        print(f"\n真值文件已保存: {args.output}")
        print(f"成功获取 {len(found_stations)}/{len(args.stations)} 个测站")

        if found_stations:
            print("\n已获取坐标的测站:")
            print(', '.join(found_stations))

    else:
        print("\n[警告] 无法下载SINEX文件")
        print("请手动从以下地址下载:")
        print("  https://cddis.nasa.gov/archive/gnss/snx/igs.snx")
        print("  或访问 https://cddis.nasa.gov/archive/gnss/snx/")

        # 生成手动下载命令
        print("\n手动下载命令示例:")
        print(f"curl -o {snx_file} https://cddis.nasa.gov/archive/gnss/snx/igs.snx")

        return 1

    return 0


if __name__ == '__main__':
    sys.exit(main())
