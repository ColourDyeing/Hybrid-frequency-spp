"""
create_overview_kml.py - 生成展示6组机型共52组数据集运动轨迹的 Google Earth KML 文件。

每个数据集单独一个 Folder，不同手机型号用不同颜色区分。
"""

import os
import sys
import io
from os.path import join, isfile, dirname, abspath
from glob import glob

# 配置参数
DATA_SET = 'train'
datapath = '../data/'
OUTPUT_FILE = '../data/all_trajectories_overview.kml'

PHONES = ['mi8', 'pixel6pro', 'pixel7pro', 's20ultra', 's21ultra', 's22ultra']

# 手机型号 -> KML 颜色 (aabbggrr格式) + 横向偏移(度)
# 颜色两两一组: 红色 ff, 绿色 00, 蓝色 00 -> ff0000ff
# 透明度 aa: ff=不透明
# 横向偏移：对轨迹整体平移一个小角度量，让同路线的多条轨迹在2D平铺显示
# 1度纬度 ≈ 111km，0.0001度 ≈ 11m，向北偏移
PHONE_COLORS = {
    'mi8':        {'color': 'ff1966ff', 'offset': 0.00000},    # 橙红色
    'pixel6pro':  {'color': 'ff00cc47', 'offset': 0.00015},    # 绿色，向北 +17m
    'pixel7pro':  {'color': 'ff0080ff', 'offset': -0.00015},   # 蓝色，向南 -17m
    's20ultra':  {'color': 'ff00e5cc', 'offset': 0.00030},    # 青色，向北 +33m
    's21ultra':  {'color': 'ffe06619', 'offset': -0.00030},    # 橙黄色，向南 -33m
    's22ultra':  {'color': 'ffcc00ff', 'offset': 0.00045},    # 紫色，向北 +50m
}
PHONE_LINE_WIDTH = 3.5

# 抽稀：每隔 N 个点取 1 个（1 = 不抽稀）
DECIMATE = 1

# 最小距离阈值（单位：度），去除过近的冗余点
# 0.00005 度 ≈ 5.5m
MIN_DIST_DEG = 0.00005

# 三次样条插值：每对相邻点之间插入 N 个中间点（0 = 不插值）
SPLINE_INTERP = 3

# 海拔模式：回到贴地模式，只做横向偏移区分颜色
ALTITUDE_MODE = 'clampToGround'


SCRIPT_DIR = dirname(abspath(__file__))
datapath = abspath(join(SCRIPT_DIR, datapath))

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', line_buffering=True)
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', line_buffering=True)


def read_ground_truth(gt_csv_path):
    """读取单个 ground_truth.csv，返回 [(lat, lon, alt, ts), ...]"""
    points = []
    with open(gt_csv_path, encoding='utf-8') as f:
        lines = f.readlines()[1:]
    for line in lines:
        if not line.strip():
            continue
        d = line.strip().split(',')
        if len(d) < 9:
            continue
        try:
            lat = float(d[2])
            lon = float(d[3])
            alt = float(d[4])
            ts = int(float(d[8]))
            points.append((lat, lon, alt, ts))
        except (ValueError, IndexError):
            continue
    return points


def dist_sq(p1, p2):
    return (p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2


def thin_points(points, min_dist_sq):
    if not points:
        return []
    thinned = [points[0]]
    for pt in points[1:]:
        if dist_sq(pt, thinned[-1]) >= min_dist_sq:
            thinned.append(pt)
    return thinned


def spline_interpolate(points, n_interp):
    if len(points) < 3 or n_interp <= 0:
        return points
    lats = [p[0] for p in points]
    lons = [p[1] for p in points]
    alts = [p[2] for p in points]
    ts_vals = [p[3] for p in points]

    try:
        from scipy.interpolate import interp1d
        lat_fn = interp1d(list(range(len(points))), lats, kind='cubic', fill_value='extrapolate')
        lon_fn = interp1d(list(range(len(points))), lons, kind='cubic', fill_value='extrapolate')
        alt_fn = interp1d(list(range(len(points))), alts, kind='cubic', fill_value='extrapolate')
    except Exception:
        return points

    result = []
    for i in range(len(points) - 1):
        for j in range(n_interp + 1):
            frac = j / (n_interp + 1)
            lat = lats[i] + frac * (lats[i + 1] - lats[i])
            lon = lons[i] + frac * (lons[i + 1] - lons[i])
            alt = alts[i] + frac * (alts[i + 1] - alts[i])
            ts = int(ts_vals[i] + frac * (ts_vals[i + 1] - ts_vals[i]))
            result.append((lat, lon, alt, ts))
    result.append(points[-1])
    return result


def process_points(points, decimate=1, min_dist_deg=0.0, spline_n=0):
    if not points:
        return []
    if decimate > 1:
        points = points[::decimate]
    if min_dist_deg > 0:
        points = thin_points(points, min_dist_deg ** 2)
    if spline_n > 0:
        points = spline_interpolate(points, spline_n)
    return points


def coords_string(points, lat_offset=0):
    lines = []
    for lat, lon, alt, _ in points:
        lines.append(f'{lon:.9f},{lat + lat_offset:.9f},{alt:.4f}')
    return '\n'.join(lines)


def build_single_kml(dataset_name, phone, points, style_config):
    """为一个 (dataset, phone) 组合生成一个 Placemark 的 KML 字符串。"""
    color = style_config['color']
    offset = style_config['offset']
    coords = coords_string(points, lat_offset=offset)
    return f'''    <Folder>
      <name>{phone} - {dataset_name}</name>
      <Placemark>
        <styleUrl>#{phone.replace(' ', '_')}_style</styleUrl>
        <LineString>
          <tessellate>1</tessellate>
          <altitudeMode>{ALTITUDE_MODE}</altitudeMode>
          <coordinates>
{coords}
          </coordinates>
        </LineString>
      </Placemark>
    </Folder>'''


def build_style_block(phone, style_config):
    """生成一个手机型号的 Style 定义。"""
    safe_id = phone.replace(' ', '_')
    color = style_config['color']
    width = PHONE_LINE_WIDTH
    return f'''    <Style id="{safe_id}_style">
      <LineStyle>
        <color>{color}</color>
        <width>{width}</width>
      </LineStyle>
    </Style>'''


def main():
    datasets_dir = join(datapath, DATA_SET)
    if not os.path.isdir(datasets_dir):
        print(f'错误：数据集目录不存在 {datasets_dir}')
        return

    datasets = sorted(os.listdir(datasets_dir))
    print(f'扫描数据集目录: {datasets_dir}')
    print(f'目标机型: {PHONES}\n')

    style_blocks = []
    folder_blocks = []
    stats = {p: 0 for p in PHONES}
    skipped = []

    for dataset in datasets:
        dataset_path = join(datasets_dir, dataset)
        if not os.path.isdir(dataset_path):
            continue

        for phone in PHONES:
            gt_file = join(dataset_path, phone, 'ground_truth.csv')
            if not isfile(gt_file):
                skipped.append((dataset, phone))
                continue

            raw_points = read_ground_truth(gt_file)
            if len(raw_points) < 2:
                skipped.append((dataset, phone))
                continue

            points = process_points(raw_points, DECIMATE, MIN_DIST_DEG, SPLINE_INTERP)

            style_config = PHONE_COLORS.get(phone, {'color': 'ffff0000', 'offset': 0})
            folder_blocks.append(build_single_kml(dataset, phone, points, style_config))
            stats[phone] += 1

    for phone in PHONES:
        style_config = PHONE_COLORS.get(phone, {'color': 'ffff0000', 'offset': 0})
        style_blocks.append(build_style_block(phone, style_config))

    # 统计摘要
    total_datasets = sum(stats.values())
    print('=' * 60)
    print('生成摘要:')
    for phone in PHONES:
        if stats[phone] > 0:
            print(f'  {phone:15s}: {stats[phone]:2d} 条轨迹')
    print(f'  {"总计":15s}: {total_datasets:2d} 条轨迹')
    if skipped:
        print(f'\n跳过 ({len(skipped)} 条):')
        for ds, ph in skipped[:10]:
            print(f'  {ds} / {ph}')
        if len(skipped) > 10:
            print(f'  ... 还有 {len(skipped) - 10} 条')

    # 颜色图例说明
    legend_items = '\n'.join(
        f'        <tr><td style="background-color:#{PHONE_COLORS[p]["color"][6:8]}{PHONE_COLORS[p]["color"][4:6]}{PHONE_COLORS[p]["color"][2:4]}{PHONE_COLORS[p]["color"][0:2]};width:20px;"></td><td>{p} ({stats[p]}条)</td></tr>'
        for p in PHONES if stats[p] > 0
    )
    description = f'''<![CDATA[
<h3>研究数据集轨迹概览</h3>
<p><b>6组机型 · {total_datasets}条轨迹 · 2D横向分层显示</b></p>
<table border="0" cellpadding="3" cellspacing="0">
{legend_items}
</table>
<hr size="1"/>
<p style="font-size:small;color:gray;">
每条轨迹代表一个(dataset, phone)组合的运动路径。<br/>
同路线不同机型轨迹横向平铺显示，颜色+偏移双重区分。<br/>
AltitudeMode: {ALTITUDE_MODE}<br/>
抽稀: {DECIMATE}x | 最小间隔: {MIN_DIST_DEG}度 | 插值: {SPLINE_INTERP}段/点
</p>
]]>'''

    kml = f'''<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
  <Document>
    <name>研究轨迹概览 - 6组机型52组数据集</name>
    <description>{description}</description>

{''.join(style_blocks)}

{''.join(folder_blocks)}

  </Document>
</kml>
'''

    output_path = abspath(join(SCRIPT_DIR, OUTPUT_FILE))
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(kml)

    print(f'\n输出文件: {output_path}')
    print(f'文件大小: {os.path.getsize(output_path) / 1024:.1f} KB')


if __name__ == '__main__':
    main()
