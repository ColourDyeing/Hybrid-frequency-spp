"""
create_kml.py - 读取真值文件和算法解算轨迹，生成 Google Earth KML 文件进行可视化。
"""

import os
import csv
from datetime import datetime, timezone
import sys
import io
from os.path import join, isfile, dirname, abspath, isabs
from glob import glob

try:
    from scipy.interpolate import interp1d
    from scipy.signal import savgol_filter
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


########### 配置参数 ##########################################

# 数据集名称，用于拼接文件名
DATA_SET = 'train'
datapath = '../data/'
OUTPUT_SUBDIR = 'kml'

# 真值文件路径，None 则自动查找 ground_truths_{DATA_SET}.csv
GT_FILE = None

# 海拔模式: 'clampToGround' | 'absolute' | 'relativeToGround'
# 推荐 clampToGround，与卫星影像贴合最稳定
ALTITUDE_MODE = 'clampToGround'

# ---- 待比较的算法列表 ----
# 每项格式: (TAG, 颜色, 线宽)
# KML 颜色格式为 aabbggrr (透明度, 蓝色, 绿色, 红色)，例如:
#   红色   ff0000ff   蓝色   ffff0000
#   绿色   ff00ff00   黄色   ff00ffff
#   紫色   ff7f00ff   橙色   ff00a0ff
SOL_ALGOS = [
    ('SPP-BRDC', 'ffff0000', 5),   # 蓝色 - spp-brdc(baseline)
    ('SPP-HFUC', 'ff00ff00', 4),   # 绿色 - spp-uduc
    ('SPP-HFUC-RM',      'ff00a0ff', 6),  # 橙色, 粗, 不透明(底层)
    # ('SPP-HFUC-RM-IGG3', 'ffffff00', 2),  # 青色, 细, 半透明(顶层)
]

# 真值轨迹的独立颜色和线宽（不受 SOL_ALGOS 影响）
GT_COLOR = 'ff0000ff'   # 红色
GT_WIDTH = 5

# 抽稀：每隔 N 个点取 1 个（1 = 不抽稀，保留所有点）
DECIMATE = 1

# 最小距离阈值（单位：度），去除过近的冗余点
# 1 度纬度 ≈ 111km，1 度经度 ≈ 85km（中纬度）
# 0.0001 度 ≈ 11m， 0.0005 度 ≈ 55m， 0.001 度 ≈ 111m
MIN_DIST_DEG = 0.0001

# 每隔 N 个点放置一个标记点（0 = 禁用）
POINT_INTERVAL = 0

# 三次样条插值：每对相邻点之间插入 N 个中间点
# 值越大曲线越密、越平滑，但 KML 文件越大
# 0 = 不插值
SPLINE_INTERP = 5

# Savitzky-Golay 平滑滤波器（0 = 禁用）
# window_length: 窗口大小，必须为奇数，越大越平滑但轨迹越"滞后"
# polyorder: 多项式阶数，必须 < window_length，建议 2-3
SAVGOL_WINDOW = 5
SAVGOL_ORDER = 3

# 使用 gx:Track 代替 LineString（Google Earth 卫星视图兼容性更好）
USE_GX_TRACK = False

#############################################################

SCRIPT_DIR = dirname(abspath(__file__))
datapath = abspath(join(SCRIPT_DIR, datapath))
output_dir = abspath(join(SCRIPT_DIR, '..', 'data', OUTPUT_SUBDIR))

# 解决Windows控制台中文编码问题
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', line_buffering=True)
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', line_buffering=True)


def read_csv(filepath):
    """读取 CSV 文件，返回 {tripId: [(lat, lon, height, timestamp), ...]} 字典。"""
    by_trip = {}
    with open(filepath, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            trip_id = row.get('tripId') or row.get('tripId'.strip())
            if not trip_id:
                continue
            lat_key = next((k for k in row if k.strip() == 'LatitudeDegrees'), None)
            lon_key = next((k for k in row if k.strip() == 'LongitudeDegrees'), None)
            h_key = next((k for k in row if k.strip() == 'Height'), None)
            ts_key = next((k for k in row if k.strip() == 'UnixTimeMillis'), None)
            if not (lat_key and lon_key and h_key):
                continue
            try:
                lat = float(row[lat_key])
                lon = float(row[lon_key])
                height = float(row[h_key])
                timestamp = int(row[ts_key]) if ts_key else None
            except ValueError:
                continue
            if trip_id not in by_trip:
                by_trip[trip_id] = []
            by_trip[trip_id].append((lat, lon, height, timestamp))
    return by_trip


def dist_sq(p1, p2):
    """计算两点在 (lat, lon) 空间的欧氏距离平方（对小范围足够精确）。"""
    return (p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2


def thin_points(points, min_dist_sq):
    """根据最小距离阈值去除相邻过近的点。"""
    if not points:
        return []
    thinned = [points[0]]
    for pt in points[1:]:
        if dist_sq(pt, thinned[-1]) >= min_dist_sq:
            thinned.append(pt)
    return thinned


def spline_interpolate(points, n_interp):
    """在每对相邻点之间均匀插入 n_interp 个中间点。

    有 scipy 时使用三次样条插值，否则回退为线性插值。
    """
    if len(points) < 3 or n_interp <= 0:
        return points

    lats = [p[0] for p in points]
    lons = [p[1] for p in points]
    alts = [p[2] for p in points]
    ts_vals = [p[3] for p in points]
    t_in = list(range(len(points)))

    if HAS_SCIPY:
        lat_fn = interp1d(t_in, lats, kind='cubic', fill_value='extrapolate')
        lon_fn = interp1d(t_in, lons, kind='cubic', fill_value='extrapolate')
        alt_fn = interp1d(t_in, alts, kind='cubic', fill_value='extrapolate')
    else:
        lat_fn = interp1d(t_in, lats, kind='linear', fill_value='extrapolate')
        lon_fn = interp1d(t_in, lons, kind='linear', fill_value='extrapolate')
        alt_fn = interp1d(t_in, alts, kind='linear', fill_value='extrapolate')

    t_out = []
    for i in range(len(points) - 1):
        for j in range(n_interp + 1):
            t_out.append(t_in[i] + j * (1.0 / (n_interp + 1)))
    t_out.append(float(t_in[-1]))

    result = []
    for t_val in t_out:
        idx = max(0, min(int(t_val), len(points) - 2))
        frac = t_val - idx
        lat = lats[idx] + frac * (lats[idx + 1] - lats[idx])
        lon = lons[idx] + frac * (lons[idx + 1] - lons[idx])
        alt = alts[idx] + frac * (alts[idx + 1] - alts[idx])
        result.append((lat, lon, alt, ts_vals[idx]))
    return result


def savgol_smooth(points, window, polyorder):
    """对经纬高序列应用 Savitzky-Golay 平滑滤波器，去除高频抖动。"""
    if len(points) < window:
        return points
    lats = [p[0] for p in points]
    lons = [p[1] for p in points]
    alts = [p[2] for p in points]
    ts_vals = [p[3] for p in points]

    lat_smooth = savgol_filter(lats, window, polyorder, mode='interp')
    lon_smooth = savgol_filter(lons, window, polyorder, mode='interp')
    alt_smooth = savgol_filter(alts, window, polyorder, mode='interp')

    return list(zip(lat_smooth, lon_smooth, alt_smooth, ts_vals))


def process_points(points, decimate=1, min_dist_deg=0.0, spline_n=0,
                   savgol_window=0, savgol_order=2):
    """依次执行：抽稀 -> 距离过滤 -> 平滑 -> 插值。"""
    if not points:
        return []
    if decimate > 1:
        points = points[::decimate]
    if min_dist_deg > 0:
        points = thin_points(points, min_dist_deg ** 2)
    if savgol_window > 3 and len(points) >= savgol_window:
        points = savgol_smooth(points, savgol_window, savgol_order)
    if spline_n > 0:
        points = spline_interpolate(points, spline_n)
    return points


def coords_string_linestring(points):
    """将点列表转换为 KML LineString coordinates 格式：lon,lat[,alt] 每行一个。"""
    lines = []
    for lat, lon, alt, _ in points:
        lines.append(f'{lon:.9f},{lat:.9f},{alt:.4f}')
    return '\n'.join(lines)


def gx_track_block(points, style_url, folder_name):
    """生成包含 gx:Track 的 Folder（Google Earth 卫星视图兼容性更好）。"""
    when_lines = []
    coord_lines = []
    for lat, lon, alt, ts in points:
        coord_lines.append(f'{lon:.9f},{lat:.9f},{alt:.4f}')
        if ts is not None:
            s = ts / 1000.0
            dt = datetime.fromtimestamp(s, tz=timezone.utc)
            when_lines.append(dt.strftime('%Y-%m-%dT%H:%M:%SZ'))
        else:
            when_lines.append('')

    track_elems = []
    for when, coord in zip(when_lines, coord_lines):
        if when:
            track_elems.append(f'          <when>{when}</when>')
        track_elems.append(f'          <gx:coord>{coord}</gx:coord>')

    track_section = '\n'.join(track_elems)

    return f'''    <Folder>
      <name>{folder_name}</name>
      <Placemark>
        <styleUrl>{style_url}</styleUrl>
        <gx:Track xmlns:gx="http://www.google.com/kml/ext/2.2">
          <extrude>1</extrude>
          <tessellate>1</tessellate>
          <altitudeMode>{ALTITUDE_MODE}</altitudeMode>
{track_section}
        </gx:Track>
      </Placemark>
    </Folder>'''


def build_kml(trip_id, gt_points, algo_pts_dict):
    """Build a complete KML string for one trip.

    algo_pts_dict: {tag: [(lat, lon, alt, ts), ...], ...}
    """
    # 处理各算法的轨迹点
    processed = {}
    for tag, pts in algo_pts_dict.items():
        processed[tag] = process_points(pts, DECIMATE, MIN_DIST_DEG,
                                       SPLINE_INTERP, SAVGOL_WINDOW, SAVGOL_ORDER)

    # 真值使用独立配置的颜色
    gt_color = GT_COLOR
    gt_width = GT_WIDTH
    gt_pts = process_points(gt_points, DECIMATE, MIN_DIST_DEG,
                           SPLINE_INTERP, SAVGOL_WINDOW, SAVGOL_ORDER)

    # 标记点（可选）
    gt_points_section = ''
    if POINT_INTERVAL > 0:
        for i, (lat, lon, alt, _) in enumerate(gt_pts):
            if i % POINT_INTERVAL == 0:
                gt_points_section += f'''      <Placemark>
        <name>GT {i}</name>
        <Point><coordinates>{lon:.9f},{lat:.9f},{alt:.4f}</coordinates></Point>
      </Placemark>
'''
        gt_points_section = f'''    <Folder>
      <name>Ground Truth Points</name>
{gt_points_section}    </Folder>
'''

    # 构建样式定义
    style_blocks = f'''    <Style id="gt_line">
      <LineStyle>
        <color>{gt_color}</color>
        <width>{gt_width}</width>
      </LineStyle>
    </Style>
'''
    for tag, color, width in SOL_ALGOS:
        style_blocks += f'''    <Style id="sol_{tag}_line">
      <LineStyle>
        <color>{color}</color>
        <width>{width}</width>
      </LineStyle>
    </Style>
'''

    # 构建轨迹段
    def make_track_section(tag, pts, is_gt=False):
        style = '#gt_line' if is_gt else f'#sol_{tag}_line'
        label = 'Ground Truth' if is_gt else tag
        coords = coords_string_linestring(pts)
        if USE_GX_TRACK:
            block = gx_track_block(pts, style, label)
        else:
            block = f'''    <Folder>
      <name>{label}</name>
      <Placemark>
        <styleUrl>{style}</styleUrl>
        <LineString>
          <tessellate>1</tessellate>
          <altitudeMode>{ALTITUDE_MODE}</altitudeMode>
          <coordinates>
{coords}
          </coordinates>
        </LineString>
      </Placemark>
    </Folder>'''
        return block

    sections = [make_track_section('', gt_pts, is_gt=True), gt_points_section]
    for tag in [t for t, _, _ in SOL_ALGOS]:
        sections.append(make_track_section(tag, processed.get(tag, [])))

    kml = f'''<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
  <Document>
    <name>{trip_id}</name>
    <description>Ground Truth (ref) vs {' / '.join([t for t,_,_ in SOL_ALGOS])}</description>

{style_blocks}

{''.join(sections)}

  </Document>
</kml>
'''
    return kml


def main():
    # 真值文件路径
    if GT_FILE:
        gt_file = abspath(join(datapath, GT_FILE)) if not isabs(GT_FILE) else GT_FILE
    else:
        gt_file = join(datapath, f'ground_truths_{DATA_SET}.csv')

    if not isfile(gt_file):
        print(f'错误：未找到真值文件 {gt_file}')
        return

    # 加载所有算法的 CSV 数据
    algo_data = {}
    for tag, _, _ in SOL_ALGOS:
        pattern = join(datapath, f'locations_{tag}_{DATA_SET}_*.csv')
        files = sorted(glob(pattern), reverse=True)
        if files:
            algo_data[tag] = read_csv(files[0])
            print(f'  {tag}: {files[0]}')
        else:
            print(f'警告：未找到 tag "{tag}" 对应的 CSV 文件')
            algo_data[tag] = {}

    gt_by_trip = read_csv(gt_file)
    print(f'真值文件: {gt_file}')

    # 收集所有行程 ID
    all_trips = set(gt_by_trip.keys())
    for data in algo_data.values():
        all_trips |= set(data.keys())
    all_trips = sorted(all_trips)

    os.makedirs(output_dir, exist_ok=True)

    for trip_id in all_trips:
        algo_pts = {tag: data.get(trip_id, []) for tag, data in algo_data.items()}
        kml_content = build_kml(trip_id, gt_by_trip.get(trip_id, []), algo_pts)
        safe_name = trip_id.replace('/', '_').replace('\\', '_')
        out_path = join(output_dir, f'{safe_name}.kml')
        with open(out_path, 'w', encoding='utf-8') as f:
            f.write(kml_content)

    print(f'\n输出目录: {output_dir}')
    print(f'算法列表: {[t for t,_,_ in SOL_ALGOS]}')
    print(f'生成 KML 文件: {len(all_trips)}')


if __name__ == '__main__':
    main()
