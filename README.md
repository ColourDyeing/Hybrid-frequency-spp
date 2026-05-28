# Hybrid-frequency-SPP

基于 RTKLIB-PRO 的多频 GNSS 单点定位（Single Point Positioning, SPP）研究项目，实现了多种单点定位算法并进行对比分析，包括：单频 SPP、双频/三频消电离层组合（IFLC）以及多频未组合+抗差估计（HFUC+RM）等方案。

## 项目成果

项目在低成本多频 GNSS 接收机（Android 智能手机）上实现了多频 SPP 算法，主要贡献包括：

- **多频未组合观测模型**：支持 L1/L2/L5 三频未组合 SPP，直接利用原始伪距和载波相位观测值进行解算，避免频率组合导致的信息损失
- **IGG-III 抗差估计**：引入中国科学院测量与地球物理研究所提出的 IGG-III 权函数方案，有效抑制粗差和周跳的影响
- **多种算法对比**：提供单频（SF）、未组合（HFUC）、抗差（HFUC-RM）、消电离层（IFLC）等多种定位算法的完整对比框架

## 项目结构

```
Hybrid-frequency-SPP/
|
├── RTKLIB-PRO/                      # 改进版 RTKLIB 源码
│   ├── app/                         # 应用程序
│   │   ├── consapp/                # 命令行应用程序
│   │   │   ├── rnx2rtkp/           # 定位解算主程序
│   │   │   ├── convbin/            # RINEX 格式转换
│   │   │   ├── pos2kml/            # 定位结果转 KML
│   │   │   ├── str2str/            # 数据流转发
│   │   │   ├── rtkrcv/             # RTK 服务器
│   │   │   └── CMakeLists.txt
│   │   ├── winapp/                 # Windows GUI 应用程序
│   │   └── qtapp/                  # Qt GUI 应用程序
│   ├── src/                        # 核心定位算法源码
│   │   ├── rtkcmn.c/h              # 公共函数库
│   │   ├── rtklib.h                # 主头文件
│   │   ├── pntpos.c/h              # 单点定位算法
│   │   ├── rtkpos.c/h              # RTK 定位算法
│   │   ├── ppp.c/h                 # 精密单点定位
│   │   ├── rinex.c/h               # RINEX 文件读写
│   │   ├── rtcm3.c/h               # RTCM3 协议
│   │   ├── ephemeris.c/h           # 星历处理
│   │   ├── ionex.c/h               # IONEX 电离层模型
│   │   ├── preceph.c/h             # 精密星历
│   │   ├── solution.c/h            # 定位结果输出
│   │   ├── convrnx.c/h             # RINEX 格式转换
│   │   ├── convkml.c/h             # KML 格式转换
│   │   ├── lambda.c/h              # LAMBDA 周跳修复
│   │   ├── rcv/                    # 接收机驱动
│   │   └── CMakeLists.txt
│   ├── lib/                        # 第三方依赖库
│   │   ├── openblas/               # OpenBLAS 高性能线性代数
│   │   ├── sofa/                   # SOFA 天文库
│   │   └── iers/                   # IERS 地球定向库
│   ├── data/                       # 配置与参考数据
│   │   ├── cmd/                    # 接收机命令文件
│   │   ├── ant/                    # 天线相位中心校正
│   │   ├── config/                 # 定位配置模板
│   │   └── pos/                    # 精密星历、IGS 站点列表
│   ├── test/                       # 单元测试与测试数据
│   ├── util/                       # 工具程序（rinex2rtcm, simobs 等）
│   └── CMakeLists.txt
│
└── Testing-Data/                   # 测试数据与处理脚本
    ├── config/                      # RTKLIB 定位配置文件
    │   └── Hybrid-frequency-spp.conf
    ├── python/                      # Python 数据处理脚本
    │   ├── 01_get_base_data.py     # 从 IGS/WHU 下载基站观测数据和导航电文
    │   ├── 02_run_multi.py         # GNSS 日志转 RINEX 并批量运行 RTKLIB 定位
    │   ├── 03_create_groundtruth_csv.py    # 从参考轨迹生成真值 CSV
    │   ├── 04_create_baseline_csv_from_pos.py  # 从定位结果提取 CSV
    │   ├── 05_Scientific_Drawing.ipynb     # 科学绘图分析（Jupyter Notebook）
    │   ├── 06_create_kml.py         # 生成轨迹 KML 文件
    │   ├── 07_Observation_Quality_Analysis.ipynb  # 观测数据质量分析
    │   ├── 08_create_overview_kml.py   # 生成轨迹总览 KML
    │   ├── gnsslogger.py            # GNSS Logger 文件解析器
    │   ├── gnsslogger_to_rnx.py     # Android GNSS 日志转 RINEX 格式
    │   └── rinex3.py                # RINEX 格式工具函数
    ├── rtklib/                      # 预编译 RTKLIB 可执行文件
    │   ├── rnx2rtkp-sf.exe          # 单频 SPP（基线对比）
    │   ├── rnx2rtkp-df.exe          # 双频 SPP
    │   ├── rnx2rtkp-hfuc.exe        # 多频未组合 SPP（HFUC）
    │   ├── rnx2rtkp-hfuc-rm.exe     # HFUC + IGG-III 抗差估计
    │   ├── rnx2rtkp-hfuc-rm-igg3.exe  # HFUC + IGG-III（别名）
    │   └── rnx2rtkp-iflc.exe        # 消电离层线性组合 SPP（IFLC）
    └── data/                        # 测试数据集
        └── train/                   # 训练数据集
            └── <date>-<location>/    # 按日期和地点组织
                └── <phone_model>/    # 按手机型号组织（mi8, pixel6pro 等）
                    ├── supplemental/   # 定位解算结果
                    └── gnss_log.txt    # 原始 GNSS 日志
```

## RTKLIB-PRO 简介

RTKLIB-PRO 是针对低成本 GNSS 接收机优化的 RTKLIB 版本，基于 RTKLIB 2.4.3 开发。

**主要特性：**

- 支持 GPS、GLONASS、Galileo、BeiDou、QZSS、SBAS、NavIC 多系统
- 支持单频（L1）、双频（L1+L2 / L1+L5）、三频（L1+L2+L5）定位模式
- 支持 Single、Kinematic、Static、PPK、PPP 等多种定位模式
- 提供 Qt GUI、Console 和 Windows GUI 三种界面
- CMake 构建系统，支持 Windows / Linux / macOS

### 本项目对 RTKLIB-PRO 的主要改进

- 新增多频未组合 SPP 定位模式（`posopt7=on`）
- 集成 IGG-III 抗差估计方案（`posopt8=on`）
- 支持智能手机多频 GNSS 原始观测值处理
- 优化 L1/L2/L5 信噪比阈值和误差模型参数

## 快速开始

### 依赖环境

| 依赖 | 版本要求 | 说明 |
|------|----------|------|
| C/C++ 编译器 | MSVC / GCC / Clang | Visual Studio 2019+ 或 GCC 7+ |
| CMake | 3.10+ | 构建工具 |
| Python | 3.8+ | 运行数据处理脚本 |
| Jupyter Notebook | - | 运行科学绘图（可选） |
| Qt5 | 5.15+ | 构建 GUI 应用程序（可选） |

### 构建 RTKLIB-PRO

```bash
cd RTKLIB-PRO
mkdir build && cd build
cmake ..
make
```

### 使用 Python 脚本进行批量定位处理

**步骤 1：下载基站数据**

编辑 `01_get_base_data.py` 中的数据路径参数后运行：

```bash
cd Testing-Data/python
python 01_get_base_data.py
```

**步骤 2：转换 GNSS 日志并批量运行定位**

编辑 `02_run_multi.py` 中的数据集路径和参数后运行：

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

### 使用预编译可执行文件

直接运行 `Testing-Data/rtklib/` 目录下的二进制文件：

```bash
# 单频 SPP（基线）
./rnx2rtkp-sf.exe -k ../config/Hybrid-frequency-spp.conf rov.obs nav.obs

# 多频未组合 SPP
./rnx2rtkp-hfuc.exe -k ../config/Hybrid-frequency-spp.conf rov.obs nav.obs

# 多频未组合 + IGG-III 抗差估计
./rnx2rtkp-hfuc-rm.exe -k ../config/Hybrid-frequency-spp.conf rov.obs nav.obs

# 消电离层组合 SPP
./rnx2rtkp-iflc.exe -k ../config/Hybrid-frequency-spp.conf rov.obs nav.obs
```

## 定位配置说明

默认配置文件 `Hybrid-frequency-spp.conf` 主要参数如下：

| 参数 | 值 | 说明 |
|------|-----|------|
| `pos1-posmode` | `single` | 单点定位模式 |
| `pos1-frequency` | `3` | 三频（L1+L2+L5） |
| `pos1-elmask` | `15°` | 截止高度角 |
| `pos1-posopt7` | `on` | 启用多频未组合 SPP 模型 |
| `pos1-posopt8` | `on` | 启用 IGG-III 抗差估计 |
| `pos1-ionoopt` | `brdc` | 广播电离层模型 |
| `pos1-tropopt` | `saas` | Saastamoinen 对流层模型 |
| `pos1-navsys` | `45` | GPS+GLONASS+Galileo+BeiDou |
| `pos1-snrmask_r` | `on` | 启用移动站信噪比阈值 |
| `stats-eratio1` | `15` | L1 伪距/载波噪声比（智能手机） |
| `stats-eratio5` | `10` | L1-L5 伪距/载波噪声比 |

> 注：配置文件中 `pos1-posopt7`（多频未组合 SPP）和 `pos1-posopt8`（IGG-III 抗差估计）为自定义扩展参数，仅在使用本项目提供的 `rnx2rtkp-hfuc*.exe` 系列可执行文件时生效。

## 算法说明

本项目实现了以下 5 种 SPP 算法：

| 算法 | 可执行文件 | 观测模型 | 抗差方案 | 说明 |
|------|-----------|----------|----------|------|
| SPP-SF | `rnx2rtkp-sf.exe` | 单频 L1 伪距 | 无 | 基线对照组 |
| SPP-DF | `rnx2rtkp-df.exe` | 双频 L1+L2 | 无 | 双频消电离层组合 |
| SPP-HFUC | `rnx2rtkp-hfuc.exe` | 多频未组合 | 无 | 保留原始观测信息 |
| SPP-HFUC-RM | `rnx2rtkp-hfuc-rm.exe` | 多频未组合 | IGG-III | 抑制粗差影响 |
| SPP-IFLC | `rnx2rtkp-iflc.exe` | 消电离层组合 | 无 | 经典双频/三频组合 |

### IGG-III 抗差估计

IGG-III 方案是中国科学院测量与地球物理研究所提出的渐变截断型权函数，根据验后残差大小自适应调整观测值权重：

$$
w_i = \begin{cases} 1 & |v_i| \leq k_0 \\ \dfrac{k_0}{|v_i|} \cdot \left( \dfrac{k_1 - |v_i|}{k_1 - k_0} \right)^2 & k_0 < |v_i| \leq k_1 \\ 0 & |v_i| > k_1 \end{cases}
$$

其中 $k_0$ 和 $k_1$ 为经验阈值（通常取 $k_0 = 2.5\sigma$，$k_1 = 3.5\sigma$）。

## 数据说明

> 原始 GNSS 观测数据因其体积较大未纳入仓库管理。请从以下公开数据源获取：
>
> - [Google GNSS Dataset](https://www.googleapis.com/storage/packages/gnss-raw-dataset)
> - [GSDC 竞赛数据集](https://github.com/koesie10/GSDC_2023)
> - [Google Developer GNSS Logger](https://developer.android.com/guide/topics/sensors/gnss)

## 致谢

- [RTKLIB](https://www.rtklib.com/) 由 Tomoji Takasu 开发
- [RTKLIB-PRO](https://github.com/rtklibexplorer/RTKLIB) 由 rtklibexplorer 维护
- Android GNSS 数据格式参考 [Google GNSS Logger](https://developer.android.com/guide/topics/sensors/gnss)

## 许可

本项目代码遵循与 RTKLIB 相同的 BSD 2-Clause License，具体请参见源码文件头部注释。
