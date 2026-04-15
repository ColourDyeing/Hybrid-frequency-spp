#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IGS测站精度详细分析脚本
分析误差分布、定位质量指标
"""

import numpy as np
from pathlib import Path


# 测站真值坐标 (ITRF2020)
TRUE_COORDS = {
    'HKSL': {'X': -2393383.179, 'Y': 5393860.917, 'Z': 2412592.153,
             'lat': 22.372002407, 'lon': 113.927989182, 'hgt': 95.268},
    'HKWS': {'X': -2430579.769, 'Y': 5374285.420, 'Z': 2418956.073,
             'lat': 22.434279919, 'lon': 114.335385105, 'hgt': 63.767},
    'CUSV': {'X': -1132915.066, 'Y': 6092528.520, 'Z': 1504633.141,
             'lat': 13.735913485, 'lon': 100.533924257, 'hgt': 74.247},
    'GUAM': {'X': -5071312.678, 'Y': 3568363.664, 'Z': 1488904.422,
             'lat': 13.589330048, 'lon': 144.868359886, 'hgt': 201.951},
    'YARR': {'X': -2389025.730, 'Y': 5043315.573, 'Z': -3078532.817,
             'lat': -29.046580107, 'lon': 115.346981097, 'hgt': 241.355},
}


def llh_to_xyz(lat, lon, hgt, a=6378137.0, f=1/298.257223563):
    """大地坐标转地心地固坐标"""
    e2 = 2 * f - f * f
    sin_lat = np.sin(np.radians(lat))
    cos_lat = np.cos(np.radians(lat))
    N = a / np.sqrt(1 - e2 * sin_lat ** 2)
    x = (N + hgt) * cos_lat * np.cos(np.radians(lon))
    y = (N + hgt) * cos_lat * np.sin(np.radians(lon))
    z = (N * (1 - e2) + hgt) * sin_lat
    return x, y, z


def parse_pos(pos_file):
    """解析pos文件"""
    solutions = []
    try:
        with open(pos_file, 'r', errors='ignore') as f:
            for line in f:
                if line.startswith('%') or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 8:
                    try:
                        solutions.append({
                            'week': int(parts[0]),
                            'sow': float(parts[1]),
                            'lat': float(parts[2]),
                            'lon': float(parts[3]),
                            'hgt': float(parts[4]),
                            'Q': int(parts[5]),
                            'ns': int(parts[6])
                        })
                    except ValueError:
                        continue
    except:
        return []
    return solutions


def analyze_station(station, data_dir, date_str='070', year='23'):
    """详细分析单个测站"""
    station_upper = station.upper()
    if station_upper not in TRUE_COORDS:
        return None

    true = TRUE_COORDS[station_upper]
    true_x, true_y, true_z = llh_to_xyz(true['lat'], true['lon'], true['hgt'])

    iflc_file = data_dir / f"{station}{date_str}0_iflc.pos"
    klb_file = data_dir / f"{station}{date_str}0_klobuchar.pos"

    iflc_sols = parse_pos(iflc_file)
    klb_sols = parse_pos(klb_file)

    print(f"\n{'='*70}")
    print(f"测站: {station_upper} (纬度: {true['lat']:.2f}°)")
    print(f"{'='*70}")

    for mode, sols in [('IFLC', iflc_sols), ('Klobuchar', klb_sols)]:
        if not sols:
            continue

        # 计算误差
        h_errors = []
        v_errors = []
        dlat_errors = []
        dlon_errors = []
        dhgt_errors = []
        ns_list = []

        for sol in sols:
            if sol['Q'] != 5:  # 只看单点定位结果
                continue

            x, y, z = llh_to_xyz(sol['lat'], sol['lon'], sol['hgt'])
            h_err = np.sqrt((x - true_x)**2 + (y - true_y)**2)
            v_err = abs(z - true_z)
            dlat = (sol['lat'] - true['lat']) * 111000  # 转换为米
            dlon = (sol['lon'] - true['lon']) * 111000 * np.cos(np.radians(true['lat']))

            h_errors.append(h_err)
            v_errors.append(v_err)
            dlat_errors.append(dlat)
            dlon_errors.append(dlon)
            dhgt_errors.append(sol['hgt'] - true['hgt'])
            ns_list.append(sol['ns'])

        h_errors = np.array(h_errors)
        v_errors = np.array(v_errors)
        dlat_errors = np.array(dlat_errors)
        dlon_errors = np.array(dlon_errors)
        dhgt_errors = np.array(dhgt_errors)

        print(f"\n{mode} 模式 (有效历元: {len(h_errors)})")
        print("-"*50)
        print(f"{'指标':<15} {'均值':<12} {'中位数':<12} {'STD':<12} {'RMS':<12}")
        print("-"*50)

        def stats(arr):
            return f"{np.mean(arr):.4f}m", f"{np.median(arr):.4f}m", f"{np.std(arr):.4f}m", f"{np.sqrt(np.mean(arr**2)):.4f}m"

        print(f"{'水平误差(2D)':<15} {' '.join(stats(h_errors))}")
        print(f"{'北向误差':<15} {' '.join(stats(dlat_errors))}")
        print(f"{'东向误差':<15} {' '.join(stats(dlon_errors))}")
        print(f"{'天向误差':<15} {' '.join(stats(dhgt_errors))}")
        print(f"{'卫星数':<15} 均值: {np.mean(ns_list):.1f}, 范围: {min(ns_list)}-{max(ns_list)}")

        # 计算百分比
        h_mean_iflc = np.mean(h_errors) if mode == 'IFLC' else None
        v_mean_iflc = np.mean(v_errors) if mode == 'IFLC' else None

    return True


def main():
    # 获取脚本所在目录
    script_dir = Path(__file__).parent.resolve()
    base_dir = script_dir.parent
    data_dir = base_dir / 'data' / 'igs'

    print("="*70)
    print(" IGS测站 IFLC vs Klobuchar 详细精度分析")
    print("="*70)

    stations = ['hksl', 'hkws', 'cusv', 'guam', 'yarr']

    for station in stations:
        analyze_station(station, data_dir)

    print("\n" + "="*70)
    print(" 分析结论")
    print("="*70)
    print("""
1. IFLC vs Klobuchar 适用场景:
   - 低纬度/高纬度（电离层活跃区）: IFLC明显优于Klobuchar
   - 中纬度（电离层平静区）: 两者差异不大，Klobuchar可能更好

2. YARR站点Klobuchar更好的原因:
   - 位于南纬29°，属于中纬度
   - 电离层延迟本身较小
   - Klobuchar模型在中纬度拟合效果好
   - IFLC组合在电离层延迟小时反而放大噪声

3. 建议:
   - 对于低纬度区域，优先使用IFLC
   - 对于中纬度区域，两者均可
   - 对于高精度需求，使用精密星历代替广播星历
""")


if __name__ == '__main__':
    main()
