"""
06_download_dcb.py - 下载AIUB CODE_MGEX DCB/BIA参数文件

通过FTP连接 ftp.aiub.unibe.ch 下载卫星偏差参数文件，
自动解压保存供RTKLIB使用。

数据来源与格式（AIUB CODE_MGEX）:
  - 2022年及之前: /CODE_MGEX/CODE/YYYY/COMwwwwD.DCB.Z
  - 2023年及之后: /CODE_MGEX/CODE/YYYY/COD0MGXFIN_YYYYDOY*_01D_01D_OSB.BIA.gz

RTKLIB的readdcb()支持: .DCB, .BIA, .BSX格式文件

使用方法:
    python 06_download_dcb.py                          # 下载所需文件
    python 06_download_dcb.py --force                 # 强制重新下载
    python 06_download_dcb.py --data_dir ../data       # 指定数据目录
    python 06_download_dcb.py --check                 # 仅检查缺失文件
"""

import os
import sys
import re
import ftplib
import io
import gzip
from os.path import join, isdir, isfile, abspath, dirname
from glob import glob
from datetime import datetime, timedelta
import argparse
import io as sysio

# 解决Windows中文编码
if sys.platform == 'win32':
    sys.stdout = sysio.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', line_buffering=True)
    sys.stderr = sysio.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', line_buffering=True)

DEFAULT_DCB_DIR = '../data/dcb'
FTP_HOST = 'ftp.aiub.unibe.ch'
FTP_BASE = '/CODE_MGEX/CODE'

# 用于解压 .DCB.Z (Unix compress格式)
try:
    from unlzw import unlzw as _unlzw
except ImportError:
    _unlzw = None

# 2023年之前的COM DCB格式截止年份
COM_DCB_CUTOFF_YEAR = 2022


def parse_nav_filename(nav_file):
    """从导航文件名解析GPS周和周内日"""
    filename = os.path.basename(nav_file)
    match = re.search(r'BRDM00DLR_S_(\d{4})(\d{3})', filename)
    if match:
        year = int(match.group(1))
        doy = int(match.group(2))
        date = datetime(year, 1, 1) + timedelta(days=doy - 1)
        gps_start = datetime(1980, 1, 6)
        delta = (date - gps_start).days
        gps_week = delta // 7
        gps_day = delta % 7
        return (year, date.month, gps_week, gps_day, date)
    return None


def get_required_gps_weeks(data_dir):
    """扫描数据目录，返回GPS周列表"""
    weeks = set()
    dataset_dir = join(data_dir, 'dataset')
    if not isdir(dataset_dir):
        return list(weeks)
    for pattern in ['BRDM*.rnx', 'BRDM*.RNX', '*.rnx', '*.RNX']:
        for nav_file in glob(join(dataset_dir, '**', pattern), recursive=True):
            result = parse_nav_filename(nav_file)
            if result:
                weeks.add((result[2], result[3]))
    return sorted(list(weeks))


def gps_to_date(gps_week, gps_day):
    """GPS周+周内日 -> datetime"""
    gps_start = datetime(1980, 1, 6)
    return gps_start + timedelta(weeks=gps_week, days=gps_day)


def is_modern_format(gps_week):
    """判断是否使用新的BIA格式（2023年之后）"""
    date = gps_to_date(gps_week, 0)
    return date.year >= COM_DCB_CUTOFF_YEAR + 1


def download_com_dcb(ftp, dcb_dir, gps_week, gps_day, force=False):
    """下载旧的COM DCB格式（2022年及之前）"""
    year = gps_to_date(gps_week, 0).year
    dcb_name = f'COM{gps_week}{gps_day}.DCB'
    dest_path = join(dcb_dir, dcb_name)

    if isfile(dest_path) and not force:
        try:
            with open(dest_path, 'r', errors='ignore') as f:
                if 'BIAS' in f.read(200).upper():
                    print(f'  已存在: {dcb_name}')
                    return True
        except:
            pass

    # 优先日文件
    for filename, ftp_subpath in [
        (f'COM{gps_week}{gps_day}.DCB.Z', f'{FTP_BASE}/{year}/COM{gps_week}{gps_day}.DCB.Z'),
    ]:
        try:
            buf = io.BytesIO()
            ftp.retrbinary(f'RETR {ftp_subpath}', buf.write)
            buf.seek(0)
            raw = buf.read()
            print(f'  下载: {filename} ({len(raw)} bytes)')
        except ftplib.error_perm:
            break  # 日文件不存在
        except Exception as e:
            print(f'  下载失败: {e}')
            return False
        else:
            # 解压并保存
            if _unlzw is None:
                print(f'  unlzw未安装，无法解压')
                return False
            result = _unlzw(raw)
            if result and b'BIAS' in result.upper():
                with open(dest_path, 'wb') as f:
                    f.write(result)
                print(f'  保存: {dest_path} ({len(result)} bytes)')
                return True
            print(f'  解压或验证失败')
            return False

    print(f'  日文件不存在，改用周文件...')

    # 回退到周文件
    year_start = year
    for y in [year, year - 1]:
        week_filename = f'COM{gps_week}.DCB.Z'
        week_ftp = f'{FTP_BASE}/{y}/{week_filename}'
        try:
            buf = io.BytesIO()
            ftp.retrbinary(f'RETR {week_ftp}', buf.write)
            buf.seek(0)
            raw = buf.read()
            print(f'  下载周文件: {week_filename} ({len(raw)} bytes)')
        except ftplib.error_perm:
            continue
        except Exception as e:
            print(f'  周文件下载失败: {e}')
            return False

        if _unlzw is None:
            print(f'  unlzw未安装')
            return False
        result = _unlzw(raw)
        if not result:
            print(f'  周文件解压失败')
            return False

        # 从周文件中提取指定日的数据
        day_content = extract_day_from_weekly(result, gps_day)
        if day_content:
            with open(dest_path, 'wb') as f:
                f.write(day_content)
            print(f'  提取并保存: {dest_path} ({len(day_content)} bytes)')
            return True
        else:
            print(f'  提取失败')
            return False

    print(f'  周文件也不存在')
    return False


def extract_day_from_weekly(data, gps_day):
    """从周DCB文件中提取指定日的数据"""
    text = data.decode('ascii', errors='replace')
    lines = text.split('\n')

    # 找到数据开始位置（跳过标题）
    data_start = 0
    for i, line in enumerate(lines):
        if 'BIAS' in line.upper() or '----' in line[:10]:
            data_start = i + 1
            break

    if data_start == 0:
        return None

    # 收集有效数据行
    data_lines = []
    for line in lines[data_start:]:
        s = line.strip()
        if not s or s.startswith('$$'):
            continue
        if 'FOR DAY' in s or 'CODE BIASES' in s.upper():
            continue
        # 数据行格式: PRN  VALUE  RMS ...
        parts = s.split()
        if len(parts) >= 2:
            try:
                float(parts[1])
                data_lines.append(line)
            except ValueError:
                pass

    if not data_lines:
        return None

    # 构建输出（保留标题，替换日期）
    header_end = data_start - 1
    header = '\n'.join(lines[:header_end + 1])
    header_lines = header.split('\n')
    new_header = []
    for hline in header_lines:
        if 'FOR DAY' in hline:
            hline = re.sub(r'FOR DAY \d+,', f'FOR DAY {gps_day},', hline)
        new_header.append(hline)

    output = '\n'.join(new_header) + '\n' + '\n'.join(data_lines) + '\n'
    return output.encode('ascii', errors='replace')


def download_modern_bia(ftp, dcb_dir, gps_week, gps_day, force=False):
    """下载新的BIA格式（2023年及之后）"""
    date = gps_to_date(gps_week, gps_day)
    year = date.year
    doy = date.timetuple().tm_yday
    bia_name = f'COD0MGXFIN_{year}{doy:03d}0000_01D_01D_OSB.BIA'
    dest_path = join(dcb_dir, bia_name)

    if isfile(dest_path) and not force:
        try:
            with open(dest_path, 'r', errors='ignore') as f:
                if 'BIAS' in f.read(200).upper():
                    print(f'  已存在: {bia_name}')
                    return True
        except:
            pass

    # 构建BIA文件名（通配符版本，用于通配符匹配）
    bia_pattern = f'COD0MGXFIN_{year}{doy:03d}*_01D_01D_OSB.BIA.gz'
    # 实际需要精确的日期范围文件名
    # 从前一天到后一天尝试（覆盖一天的数据）
    for day_offset in range(-1, 2):
        test_date = date + timedelta(days=day_offset)
        test_year = test_date.year
        test_doy = test_date.timetuple().tm_yday
        filename = f'COD0MGXFIN_{test_year}{test_doy:03d}0000_01D_01D_OSB.BIA.gz'
        ftp_path = f'{FTP_BASE}/{test_year}/{filename}'

        print(f'  尝试: {filename}')
        try:
            buf = io.BytesIO()
            ftp.retrbinary(f'RETR {ftp_path}', buf.write)
            buf.seek(0)
            raw = buf.read()
            print(f'  下载成功: {len(raw)} bytes')
        except ftplib.error_perm:
            continue
        except Exception as e:
            print(f'  下载失败: {e}')
            return False

        # 解压gzip
        try:
            decompressed = gzip.decompress(raw)
            print(f'  解压成功: {len(decompressed)} bytes')
        except Exception as e:
            print(f'  解压失败: {e}')
            return False

        # 验证BIA格式
        header = decompressed[:200].decode('ascii', errors='replace')
        if 'BIAS-SINEX' not in header and 'BIAS' not in header.upper():
            print(f'  文件格式无效')
            return False

        # 保存为.BIA文件（RTKLIB的readdcb()支持）
        with open(dest_path, 'wb') as f:
            f.write(decompressed)
        print(f'  保存: {dest_path} ({len(decompressed)} bytes)')
        return True

    print(f'  未找到对应日期的BIA文件')
    return False


def check_and_download(ftp, dcb_dir, gps_week, gps_day, force=False):
    """根据GPS周选择合适的数据源

    注意: RTKLIB配置文件同时配置了COM*.DCB和COD0MGXFIN*.BIA两个路径，
    因此需要同时下载两种格式。2023年之后的BIA文件包含所有卫星系统的DCB，
    而2022年及之前的COM.DCB包含GPS/GLO的DCB。
    """
    success = False
    # 始终尝试下载COM.DCB（适用于所有年份，GPS/GLO DCB）
    if download_com_dcb(ftp, dcb_dir, gps_week, gps_day, force):
        success = True
    # 始终尝试下载BIA（2023年后包含所有卫星系统DCB）
    if download_modern_bia(ftp, dcb_dir, gps_week, gps_day, force):
        success = True
    return success


def download_dcb_files(data_dir, force=False):
    """主函数：下载所有需要的DCB/BIA文件"""
    script_dir = dirname(abspath(__file__))
    dcb_dir = abspath(join(script_dir, DEFAULT_DCB_DIR))

    if not isdir(dcb_dir):
        os.makedirs(dcb_dir)

    print('=' * 60)
    print('DCB/BIA文件下载工具 (FTP直连 AIUB)')
    print('=' * 60)

    print('\n扫描数据集...')
    weeks = get_required_gps_weeks(data_dir)
    if not weeks:
        print('未找到导航文件')
        return False

    # 统计新旧格式分布
    modern = [w for w in weeks if is_modern_format(w[0])]
    legacy = [w for w in weeks if not is_modern_format(w[0])]
    print(f'找到 {len(weeks)} 个GPS周:')
    print(f'  COM DCB格式（≤2022）: {len(legacy)} 个')
    print(f'  BIA格式（≥2023）: {len(modern)} 个')

    print('\n连接FTP服务器...')
    try:
        ftp = ftplib.FTP(FTP_HOST, timeout=60)
        ftp.login()
        print(f'连接成功! 当前目录: {ftp.pwd()}')
    except Exception as e:
        print(f'FTP连接失败: {e}')
        return False

    success_count = 0
    fail_count = 0

    for gps_week, gps_day in weeks:
        date = gps_to_date(gps_week, gps_day)
        print(f'\n[GPS周 {gps_week} 日{gps_day}] = {date.date()}')
        if check_and_download(ftp, dcb_dir, gps_week, gps_day, force):
            success_count += 1
        else:
            fail_count += 1

    ftp.quit()

    print('\n' + '=' * 60)
    print(f'下载完成: 成功 {success_count}, 失败 {fail_count}')
    print(f'文件目录: {dcb_dir}')

    existing = [f for f in os.listdir(dcb_dir)
                if f.endswith(('.DCB', '.BIA'))]
    if existing:
        print(f'\n已下载的文件 ({len(existing)} 个):')
        for f in sorted(existing)[:20]:
            size = os.path.getsize(join(dcb_dir, f))
            print(f'  {f} ({size} bytes)')

    return fail_count == 0


def check_dcb_files(data_dir):
    """检查缺失的文件

    注意: RTKLIB配置文件同时配置了COM*.DCB和COD0MGXFIN*.BIA两个路径，
    因此需要同时检查两种格式是否完整。
    """
    script_dir = dirname(abspath(__file__))
    dcb_dir = abspath(join(script_dir, DEFAULT_DCB_DIR))
    weeks = get_required_gps_weeks(data_dir)
    missing = []

    for gps_week, gps_day in weeks:
        date = gps_to_date(gps_week, gps_day)
        year = date.year
        doy = date.timetuple().tm_yday
        # 检查COM.DCB
        com_filename = f'COM{gps_week}{gps_day}.DCB'
        if not isfile(join(dcb_dir, com_filename)):
            missing.append((gps_week, gps_day, com_filename))
        # 检查BIA
        bia_filename = f'COD0MGXFIN_{year}{doy:03d}0000_01D_01D_OSB.BIA'
        if not isfile(join(dcb_dir, bia_filename)):
            missing.append((gps_week, gps_day, bia_filename))

    return missing


def main():
    parser = argparse.ArgumentParser(
        description='下载AIUB CODE_MGEX DCB/BIA参数文件',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--data_dir', '-d', default='../data',
                        help='数据集根目录 (默认: ../data)')
    parser.add_argument('--force', '-f', action='store_true',
                        help='强制重新下载')
    parser.add_argument('--check', '-c', action='store_true',
                        help='仅检查缺失文件')

    args = parser.parse_args()
    script_dir = dirname(abspath(__file__))
    data_dir = abspath(join(script_dir, args.data_dir))

    if args.check:
        missing = check_dcb_files(data_dir)
        if not missing:
            print('所有文件已完整')
        else:
            print(f'缺少 {len(missing)} 个文件:')
            for w, d, fn in missing:
                print(f'  GPS周{w}日{d}: {fn}')
    else:
        download_dcb_files(data_dir, args.force)


if __name__ == '__main__':
    main()
