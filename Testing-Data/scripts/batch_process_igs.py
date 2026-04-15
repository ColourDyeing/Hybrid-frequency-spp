#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IGS测站批量解算脚本
对多个测站同时运行IFLC和Klobuchar模式解算
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path


def main():
    # 获取脚本所在目录的父目录
    script_dir = Path(__file__).parent.resolve()
    base_dir = script_dir.parent  # Testing-Data目录

    parser = argparse.ArgumentParser(description='IGS测站批量解算')
    parser.add_argument('--data-dir', '-d',
                        default=str(base_dir / 'data' / 'igs'),
                        help='数据目录路径')
    parser.add_argument('--rtklib-dir', '-r',
                        default=str(base_dir / 'rtklib'),
                        help='RTKLIB目录')
    parser.add_argument('--config-dir', '-c',
                        default=str(base_dir / 'config'),
                        help='配置文件目录')
    parser.add_argument('--stations', '-s', nargs='+',
                        default=['hksl', 'hkws', 'cusv', 'guam', 'yarr'],
                        help='测站列表')
    parser.add_argument('--year', '-y', default='23', help='年份后两位')
    parser.add_argument('--doy', default='070', help='年日期(DOY)')
    args = parser.parse_args()

    # 路径设置
    data_dir = Path(args.data_dir)
    rtklib_dir = Path(args.rtklib_dir)
    config_dir = Path(args.config_dir)

    # 可执行文件路径
    exe = rtklib_dir / 'rnx2rtkp-test.exe'
    if not exe.exists():
        print(f"[错误] 找不到RTKLIB程序: {exe}")
        return 1

    # 配置文件路径
    iflc_config = config_dir / 'Hybrid-frequency-spp-iono-free.conf'
    klb_config = config_dir / 'Hybrid-frequency-spp-klobuchar.conf'

    if not iflc_config.exists():
        print(f"[错误] 找不到IFLC配置文件: {iflc_config}")
        return 1
    if not klb_config.exists():
        print(f"[错误] 找不到Klobuchar配置文件: {klb_config}")
        return 1

    print("=" * 60)
    print(" IGS测站 IFLC vs Klobuchar 批量解算")
    print("=" * 60)
    print(f"RTKLIB路径: {exe}")
    print(f"数据目录: {data_dir}")
    print(f"IFLC配置: {iflc_config.name}")
    print(f"Klobuchar配置: {klb_config.name}")
    print(f"测站列表: {', '.join([s.upper() for s in args.stations])}")
    print(f"日期: 20{args.year}年第{args.doy}天")
    print()

    # 统计
    total = len(args.stations)
    success = 0
    failed = 0

    for i, station in enumerate(args.stations, 1):
        obs_file = data_dir / f"{station}{args.doy}0.{args.year}o"
        nav_file = data_dir / f"brdc{args.doy}0.{args.year}p"

        print(f"[{i}/{total}] 处理测站: {station.upper()}")
        print(f"  观测文件: {obs_file.name}")

        # 检查观测文件
        if not obs_file.exists():
            print(f"  [跳过] 观测文件不存在")
            failed += 1
            print()
            continue

        # 检查导航文件
        if not nav_file.exists():
            print(f"  [跳过] 导航文件不存在")
            failed += 1
            print()
            continue

        # 解算IFLC模式
        iflc_out = data_dir / f"{station}{args.doy}0_iflc.pos"
        print(f"  解算IFLC模式...")

        cmd_iflc = [
            str(exe),
            '-k', str(iflc_config),
            '-o', str(iflc_out),
            str(obs_file),
            str(nav_file)
        ]

        try:
            result = subprocess.run(cmd_iflc, capture_output=True, text=True, timeout=300)
            if result.returncode == 0 and iflc_out.exists():
                print(f"    [完成] IFLC解算")
            else:
                print(f"    [失败] IFLC解算失败")
                print(f"    {result.stderr[:200] if result.stderr else ''}")
                failed += 1
                print()
                continue
        except subprocess.TimeoutExpired:
            print(f"    [失败] IFLC解算超时")
            failed += 1
            print()
            continue
        except Exception as e:
            print(f"    [失败] {e}")
            failed += 1
            print()
            continue

        # 解算Klobuchar模式
        klb_out = data_dir / f"{station}{args.doy}0_klobuchar.pos"
        print(f"  解算Klobuchar模式...")

        cmd_klb = [
            str(exe),
            '-k', str(klb_config),
            '-o', str(klb_out),
            str(obs_file),
            str(nav_file)
        ]

        try:
            result = subprocess.run(cmd_klb, capture_output=True, text=True, timeout=300)
            if result.returncode == 0 and klb_out.exists():
                print(f"    [完成] Klobuchar解算")
                success += 1
            else:
                print(f"    [失败] Klobuchar解算失败")
                print(f"    {result.stderr[:200] if result.stderr else ''}")
                failed += 1
        except subprocess.TimeoutExpired:
            print(f"    [失败] Klobuchar解算超时")
            failed += 1
        except Exception as e:
            print(f"    [失败] {e}")
            failed += 1

        print(f"  [成功] {station.upper()} 解算完成")
        print()

    # 汇总
    print("=" * 60)
    print(" 批量解算完成！")
    print("=" * 60)
    print(f"总测站数: {total}")
    print(f"成功: {success}")
    print(f"失败/跳过: {failed}")
    print()
    print(f"解算结果保存在: {data_dir}")
    print(f"  * *_iflc.pos     - IFLC无电离层组合结果")
    print(f"  * *_klobuchar.pos - Klobuchar模型结果")

    return 0


if __name__ == '__main__':
    sys.exit(main())
