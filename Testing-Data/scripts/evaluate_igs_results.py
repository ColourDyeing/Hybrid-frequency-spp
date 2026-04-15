#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IGS测站 IFLC与Klobuchar精度评定脚本
对批量解算的结果进行精度对比分析
"""

import os
import sys
import glob
import argparse
import numpy as np
from datetime import datetime
from pathlib import Path

# 测站真值坐标 (ITRF2020)
# 格式: 测站名 -> [X(m), Y(m), Z(m), lat(deg), lon(deg), hgt(m)]
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

# 所有待处理的测站列表
ALL_STATIONS = ['hksl', 'hkws', 'cusv', 'guam', 'yarr']


def parse_rtklib_pos(pos_file):
    """解析RTKLIB输出的pos文件，返回定位结果列表"""
    solutions = []
    try:
        with open(pos_file, 'r', errors='ignore') as f:
            for line in f:
                line = line.strip()
                # 跳过注释行和表头
                if line.startswith('%') or not line:
                    continue
                parts = line.split()
                # 检查是否是有效数据行
                # 格式: GPST_week sow lat lon height Q ns sdn sde sdu ...
                if len(parts) >= 8:
                    try:
                        week = int(parts[0])
                        sow = float(parts[1])
                        lat = float(parts[2])
                        lon = float(parts[3])
                        hgt = float(parts[4])
                        q = int(parts[5])
                        ns = int(parts[6])
                        solutions.append({
                            'week': week,
                            'sow': sow,
                            'lat': lat,
                            'lon': lon,
                            'hgt': hgt,
                            'Q': q,
                            'ns': ns
                        })
                    except ValueError:
                        continue
    except FileNotFoundError:
        print(f"  [错误] 文件不存在: {pos_file}")
        return []
    except Exception as e:
        print(f"  [错误] 解析失败: {e}")
        return []
    return solutions


def llh_to_xyz(lat, lon, hgt, a=6378137.0, f=1/298.257223563):
    """将大地坐标转换为地心地固坐标"""
    e2 = 2 * f - f * f
    sin_lat = np.sin(np.radians(lat))
    cos_lat = np.cos(np.radians(lat))
    N = a / np.sqrt(1 - e2 * sin_lat ** 2)

    x = (N + hgt) * cos_lat * np.cos(np.radians(lon))
    y = (N + hgt) * cos_lat * np.sin(np.radians(lon))
    z = (N * (1 - e2) + hgt) * sin_lat
    return x, y, z


def calculate_errors(solutions, true_coords, q_filter=None):
    """计算定位误差

    Args:
        solutions: 定位结果列表
        true_coords: 真值坐标字典
        q_filter: Q值过滤器，如[5]表示只统计单点定位结果

    Returns:
        errors: 误差数组 [horizontal_2d, vertical_up, total_3d]
    """
    if not solutions:
        return np.array([]), np.array([]), np.array([])

    # 提取真值坐标
    true_x = true_coords['X']
    true_y = true_coords['Y']
    true_z = true_coords['Z']

    horizontal_errors = []
    vertical_errors = []
    total_errors = []

    for sol in solutions:
        # 过滤Q值
        if q_filter and sol['Q'] not in q_filter:
            continue

        # 转换解算结果到XYZ
        x, y, z = llh_to_xyz(sol['lat'], sol['lon'], sol['hgt'])

        # 计算误差
        dx = x - true_x
        dy = y - true_y
        dz = z - true_z

        # 水平误差（二维）
        h_error = np.sqrt(dx**2 + dy**2)
        # 垂直误差
        v_error = abs(dz)
        # 总误差（三维）
        t_error = np.sqrt(dx**2 + dy**2 + dz**2)

        horizontal_errors.append(h_error)
        vertical_errors.append(v_error)
        total_errors.append(t_error)

    return (np.array(horizontal_errors),
            np.array(vertical_errors),
            np.array(total_errors))


def calculate_statistics(errors):
    """计算误差统计指标"""
    if len(errors) == 0:
        return None

    return {
        'count': len(errors),
        'mean': np.mean(errors),
        'median': np.median(errors),
        'std': np.std(errors),
        'min': np.min(errors),
        'max': np.max(errors),
        'p50': np.percentile(errors, 50),
        'p95': np.percentile(errors, 95),
        'p99': np.percentile(errors, 99),
        'rmse': np.sqrt(np.mean(errors**2))
    }


def format_stats(stats, name="误差"):
    """格式化统计结果"""
    if stats is None:
        return f"  {name}: 无数据\n"
    return (
        f"  {name}:\n"
        f"    样本数: {stats['count']}\n"
        f"    均值:   {stats['mean']:.4f} m\n"
        f"    中位数: {stats['median']:.4f} m\n"
        f"    标准差: {stats['std']:.4f} m\n"
        f"    50分位: {stats['p50']:.4f} m\n"
        f"    95分位: {stats['p95']:.4f} m\n"
        f"    99分位: {stats['p99']:.4f} m\n"
        f"    最大值: {stats['max']:.4f} m\n"
        f"    最小值: {stats['min']:.4f} m\n"
        f"    RMSE:   {stats['rmse']:.4f} m\n"
    )


def evaluate_station(station, data_dir, date_str='070', year='23', output_file=None):
    """评估单个测站的定位精度"""
    obs_file = f"{data_dir}\\{station}{date_str}0.{year}o"
    if not os.path.exists(obs_file):
        print(f"[跳过] {station.upper()}: 观测文件不存在")
        return None

    # 检查真值是否存在
    station_upper = station.upper()
    if station_upper not in TRUE_COORDS:
        print(f"[跳过] {station.upper()}: 真值未定义")
        return None

    true_coords = TRUE_COORDS[station_upper]
    print(f"\n{'='*60}")
    print(f"测站: {station.upper()}")
    print(f"真值: lat={true_coords['lat']:.6f} lon={true_coords['lon']:.6f} hgt={true_coords['hgt']:.3f}")
    print(f"{'='*60}")

    # 解析IFLC结果
    iflc_file = f"{data_dir}\\{station}{date_str}0_iflc.pos"
    klb_file = f"{data_dir}\\{station}{date_str}0_klobuchar.pos"

    print("\n解析IFLC结果...")
    iflc_sols = parse_rtklib_pos(iflc_file)
    print(f"  IFLC有效历元: {len(iflc_sols)}")

    print("解析Klobuchar结果...")
    klb_sols = parse_rtklib_pos(klb_file)
    print(f"  Klobuchar有效历元: {len(klb_sols)}")

    # 计算误差
    print("计算误差统计...")
    iflc_h, iflc_v, iflc_t = calculate_errors(iflc_sols, true_coords, q_filter=[5])
    klb_h, klb_v, klb_t = calculate_errors(klb_sols, true_coords, q_filter=[5])

    # 计算统计指标
    iflc_h_stats = calculate_statistics(iflc_h)
    iflc_v_stats = calculate_statistics(iflc_v)
    iflc_t_stats = calculate_statistics(iflc_t)

    klb_h_stats = calculate_statistics(klb_h)
    klb_v_stats = calculate_statistics(klb_v)
    klb_t_stats = calculate_statistics(klb_t)

    # 打印结果
    print("\n" + "-"*60)
    print("IFLC 无电离层组合")
    print("-"*60)
    print(format_stats(iflc_h_stats, "水平误差 (2D)"))
    print(format_stats(iflc_v_stats, "垂直误差 (Up)"))
    print(format_stats(iflc_t_stats, "三维误差 (3D)"))

    print("-"*60)
    print("Klobuchar 广播模型")
    print("-"*60)
    print(format_stats(klb_h_stats, "水平误差 (2D)"))
    print(format_stats(klb_v_stats, "垂直误差 (Up)"))
    print(format_stats(klb_t_stats, "三维误差 (3D)"))

    # 对比分析
    print("-"*60)
    print("精度对比 (IFLC vs Klobuchar)")
    print("-"*60)
    print(f"{'指标':<15} {'IFLC(m)':<12} {'Klobuchar(m)':<14} {'改善(%)':<10}")
    print("-"*60)

    improvements = {}
    if iflc_h_stats and klb_h_stats:
        imp = (klb_h_stats['mean'] - iflc_h_stats['mean']) / klb_h_stats['mean'] * 100
        improvements['水平均值'] = imp
        print(f"{'水平均值':<15} {iflc_h_stats['mean']:<12.4f} {klb_h_stats['mean']:<14.4f} {imp:>+8.2f}%")

        imp = (klb_h_stats['std'] - iflc_h_stats['std']) / klb_h_stats['std'] * 100
        improvements['水平STD'] = imp
        print(f"{'水平STD':<15} {iflc_h_stats['std']:<12.4f} {klb_h_stats['std']:<14.4f} {imp:>+8.2f}%")

        imp = (klb_h_stats['p95'] - iflc_h_stats['p95']) / klb_h_stats['p95'] * 100
        improvements['水平P95'] = imp
        print(f"{'水平P95':<15} {iflc_h_stats['p95']:<12.4f} {klb_h_stats['p95']:<14.4f} {imp:>+8.2f}%")

    if iflc_v_stats and klb_v_stats:
        imp = (klb_v_stats['mean'] - iflc_v_stats['mean']) / klb_v_stats['mean'] * 100
        improvements['垂直均值'] = imp
        print(f"{'垂直均值':<15} {iflc_v_stats['mean']:<12.4f} {klb_v_stats['mean']:<14.4f} {imp:>+8.2f}%")

        imp = (klb_v_stats['std'] - iflc_v_stats['std']) / iflc_v_stats['std'] * 100
        improvements['垂直STD'] = imp
        print(f"{'垂直STD':<15} {iflc_v_stats['std']:<12.4f} {klb_v_stats['std']:<14.4f} {imp:>+8.2f}%")

        imp = (klb_v_stats['p95'] - iflc_v_stats['p95']) / klb_v_stats['p95'] * 100
        improvements['垂直P95'] = imp
        print(f"{'垂直P95':<15} {iflc_v_stats['p95']:<12.4f} {klb_v_stats['p95']:<14.4f} {imp:>+8.2f}%")

    if iflc_t_stats and klb_t_stats:
        imp = (klb_t_stats['mean'] - iflc_t_stats['mean']) / klb_t_stats['mean'] * 100
        improvements['3D均值'] = imp
        print(f"{'3D均值':<15} {iflc_t_stats['mean']:<12.4f} {klb_t_stats['mean']:<14.4f} {imp:>+8.2f}%")

    return {
        'station': station.upper(),
        'iflc': {'horizontal': iflc_h_stats, 'vertical': iflc_v_stats, 'total': iflc_t_stats},
        'klb': {'horizontal': klb_h_stats, 'vertical': klb_v_stats, 'total': klb_t_stats},
        'improvements': improvements
    }


def generate_summary_report(results):
    """生成汇总报告"""
    print("\n" + "="*100)
    print(" " * 30 + "汇总报告 - 所有测站精度对比")
    print("="*100)

    # 统计各测站的改善情况
    valid_results = [r for r in results if r is not None]

    if not valid_results:
        print("\n无有效结果！")
        return

    print(f"\n处理测站数: {len(valid_results)}/{len(ALL_STATIONS)}\n")

    # 表格表头
    print(f"{'测站':<8} {'IFLC水平(m)':<14} {'Klb水平(m)':<14} {'水平改善(%)':<14} {'IFLC垂直(m)':<14} {'Klb垂直(m)':<14} {'垂直改善(%)':<14}")
    print("-"*106)

    total_h_iflc, total_h_klb = [], []
    total_v_iflc, total_v_klb = [], []

    for r in valid_results:
        station = r['station']
        iflc_h = r['iflc']['horizontal']
        klb_h = r['klb']['horizontal']
        iflc_v = r['iflc']['vertical']
        klb_v = r['klb']['vertical']

        if iflc_h and klb_h and iflc_h['mean'] > 0:
            h_imp = (klb_h['mean'] - iflc_h['mean']) / klb_h['mean'] * 100
            total_h_iflc.append(iflc_h['mean'])
            total_h_klb.append(klb_h['mean'])
        else:
            h_imp = float('nan')

        if iflc_v and klb_v and iflc_v['mean'] > 0:
            v_imp = (klb_v['mean'] - iflc_v['mean']) / klb_v['mean'] * 100
            total_v_iflc.append(iflc_v['mean'])
            total_v_klb.append(klb_v['mean'])
        else:
            v_imp = float('nan')

        iflc_h_str = f"{iflc_h['mean']:.4f}" if iflc_h else "N/A"
        klb_h_str = f"{klb_h['mean']:.4f}" if klb_h else "N/A"
        iflc_v_str = f"{iflc_v['mean']:.4f}" if iflc_v else "N/A"
        klb_v_str = f"{klb_v['mean']:.4f}" if klb_v else "N/A"

        print(f"{station:<8} {iflc_h_str:<14} {klb_h_str:<14} {h_imp:>+12.2f}% {iflc_v_str:<14} {klb_v_str:<14} {v_imp:>+12.2f}%")

    # 计算平均改善
    print("-"*106)
    if total_h_iflc:
        avg_h_iflc = np.mean(total_h_iflc)
        avg_h_klb = np.mean(total_h_klb)
        avg_h_imp = (avg_h_klb - avg_h_iflc) / avg_h_klb * 100

        avg_v_iflc = np.mean(total_v_iflc)
        avg_v_klb = np.mean(total_v_klb)
        avg_v_imp = (avg_v_klb - avg_v_iflc) / avg_v_klb * 100

        print(f"{'平均':<8} {avg_h_iflc:<14.4f} {avg_h_klb:<14.4f} {avg_h_imp:>+12.2f}% {avg_v_iflc:<14.4f} {avg_v_klb:<14.4f} {avg_v_imp:>+12.2f}%")


def main():
    # 获取脚本所在目录的父目录
    script_dir = Path(__file__).parent.resolve()
    base_dir = script_dir.parent  # Testing-Data目录

    parser = argparse.ArgumentParser(description='IGS测站IFLC与Klobuchar精度评定')
    parser.add_argument('--data-dir', '-d', default=str(base_dir / 'data' / 'igs'),
                        help='数据目录路径')
    parser.add_argument('--date', default='070', help='年日期(DOY),如070表示第70天')
    parser.add_argument('--year', '-y', default='23', help='年份后两位,如23表示2023年')
    parser.add_argument('--stations', '-s', nargs='+', default=None,
                        help='指定测站列表,如 hksl hkws')
    parser.add_argument('--all', '-a', action='store_true',
                        help='处理所有已定义的测站')
    args = parser.parse_args()

    print("="*100)
    print(" IGS测站 IFLC vs Klobuchar 精度评定程序")
    print("="*100)
    print(f"数据目录: {args.data_dir}")
    print(f"日期: 20{args.year}年第{args.date}天")

    # 确定要处理的测站列表
    if args.all:
        stations = ALL_STATIONS
    elif args.stations:
        stations = args.stations
    else:
        # 只处理有真值的测站
        stations = ALL_STATIONS

    print(f"待处理测站: {', '.join([s.upper() for s in stations])}")

    # 批量评估
    results = []
    for station in stations:
        result = evaluate_station(station, args.data_dir, args.date, args.year)
        results.append(result)

    # 生成汇总报告
    generate_summary_report(results)

    print("\n" + "="*100)
    print(" 评定完成！")
    print("="*100)


if __name__ == '__main__':
    main()
