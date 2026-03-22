# Hybrid-frequency-SPP

基于 RTKLIB-PRO 的多频点单点定位（SPP）研究项目，包含改进的 RTKLIB 源码、数据处理脚本及实验配置。

## 项目结构

```
Hybrid-frequency-SPP/
├── RTKLIB-PRO/           # 改进版 RTKLIB 源码
│   ├── app/              # 应用程序（consapp/qtapp/winapp）
│   ├── src/              # 核心定位算法源码
│   ├── lib/              # 第三方依赖库（OpenBLAS, SOFA, IERS）
│   ├── data/             # 卫星星历、精密星历、TLE 等数据
│   ├── doc/              # 文档
│   ├── test/             # 测试数据
│   ├── util/             # 工具程序（rinex2rtcm, simobs 等）
│   └── CMakeLists.txt    # CMake 构建配置
│
└── Testing-Data/        # 测试数据与脚本
    ├── config/           # RTKLIB 定位配置文件
    ├── python/            # Python 数据处理脚本
    │   ├── 01_get_base_data.py              # 从 IGS/WHU 下载基站观测数据和导航电文
    │   ├── 02_run_multi.py                  # GNSS 日志转 RINEX 并批量运行 RTKLIB 定位
    │   ├── 03_create_groundtruth_csv.py     # 从参考轨迹文件生成真值 CSV
    │   ├── 04_create_baseline_csv_from_pos.py # 从 RTKLIB 结果提取定位 CSV
    │   ├── 05_Scientific_Drawing.ipynb      # 科学绘图（Jupyter Notebook）
    │   ├── gnsslogger.py                    # GNSS Logger 文件解析
    │   ├── gnsslogger_to_rnx.py             # Android GNSS 日志转 RINEX 格式
    │   └── rinex3.py                        # RINEX 格式工具
    └── rtklib/            # 预编译 RTKLIB 可执行文件（crx2rnx, rnx2rtkp）
```

## RTKLIB-PRO 简介

RTKLIB-PRO 是针对低成本 GNSS 接收机（单频、双频或三频，尤其 u-blox 接收机和 Android 手机）优化的 RTKLIB 版本，基于 RTKLIB 2.4.3 开发。

**主要特性：**

- 支持 GPS、GLONASS、Galileo、BeiDou、QZSS、SBAS、NavIC 多系统
- 支持单频（L1）、双频（L1+L2 / L1+L5）、三频（L1+L2+L5）定位模式
- 支持 Single、Kinematic、Static、PPK、PPP 等多种定位模式
- 提供 GUI（Qt）和命令行（ConsApp）两种界面
- CMake 构建系统，支持 Windows / Linux

### RTKLIB-PRO 与官方 RTKLIB 的主要区别

- 针对 Android GNSS 原始观测值的特殊处理
- 优化了低质量接收机的定位算法
- 支持多频点组合定位（Hybrid-frequency）
- 增强了周跳检测与修复能力

## 快速开始

### 依赖环境

- **C/C++ 编译器**：MSVC (Visual Studio)、GCC、Clang
- **CMake**：3.10+
- **Qt5**（如需 GUI）
- **Python**：3.8+（运行数据处理脚本）
- **RTKLIB 可选依赖**：OpenBLAS（高性能线性代数库）

### 构建 RTKLIB-PRO（CMake）

```bash
cd RTKLIB-PRO
mkdir build && cd build
cmake ..
make
```

生成的可执行文件位于 `RTKLIB-PRO/bin/` 目录。

### 使用 Python 脚本进行批量定位处理

**步骤 1：下载基站数据**

修改 `01_get_base_data.py` 中的数据路径参数，然后运行：

```bash
cd Testing-Data/python
python 01_get_base_data.py
```

**步骤 2：转换 GNSS 日志并运行定位**

修改 `02_run_multi.py` 中的数据集路径和参数，然后运行：

```bash
python 02_run_multi.py
```

**步骤 3：提取定位结果**

```bash
# 生成真值 CSV
python 03_create_groundtruth_csv.py

# 从定位结果提取 CSV
python 04_create_baseline_csv_from_pos.py
```

**步骤 4：绘图分析**

打开 Jupyter Notebook 运行 `05_Scientific_Drawing.ipynb`。

### 使用预编译二进制文件

直接运行 `Testing-Data/rtklib/` 目录下的可执行文件：

```bash
# 格式转换
./crx2rnx.exe input.crx output.rnx

# 定位解算
./rnx2rtkp.exe -k ../config/Hybrid-frequency-spp.conf rov.obs base.obs
```

## 定位配置说明

项目默认使用 **GPS L1 单频单点定位（Single Point Positioning）** 配置，主要参数如下：

| 参数 | 值 | 说明 |
|------|-----|------|
| posmode | single | 单点定位模式 |
| frequency | l1 | L1 单频 |
| elmask | 15° | 截止高度角 |
| navsys | GPS only | 仅使用 GPS 系统 |
| ionoopt | brdc | 使用广播电离层模型 |
| tropopt | saas | 使用 Saastamoinen 对流层模型 |

如需使用多频或 PPK 模式，请参考 `RTKLIB-PRO/data/cmd/` 中的示例命令文件进行配置。

## 数据说明

> **注意**：原始 GNSS 观测数据（train/test 等目录）因其体积较大未纳入仓库管理。请根据需要自行下载或从 [Google GNSS Dataset](https://www.googleapis.com/storage/packages/gnss-raw-dataset)、[GSDC 竞赛数据集](https://github.com/koesie10/GSDC_2023) 等公开数据源获取。

## 致谢

- RTKLIB 由 [Tomoji Takasu](https://www.rtklib.com/) 开发
- RTKLIB-PRO 由 [rtklibexplorer](https://github.com/rtklibexplorer/RTKLIB) 维护
- Android GNSS 数据格式参考 [Google GNSS Logger](https://developer.android.com/guide/topics/sensors/gnss)

## 许可

本项目代码遵循与 RTKLIB 相同的许可协议，具体请参见各源码文件头部注释。
