# 混合频率智能手机单点定位算法设计

本仓库是本科毕业设计《混合频率智能手机单点定位算法设计》的实现代码与实验处理工具，GitHub 仓库名为 **Hybrid-frequency-spp**。

研究面向智能手机 GNSS 观测中单频与双频数据并存、L5 观测不完整以及观测噪声较大的问题，基于 RTKLIB 二次开发，实现 **L1/L5 混合频率非组合伪距单点定位（SPP）**，并进一步引入频率—高度角—信噪比联合定权和自适应 IGG-III 抗差估计。实验使用 GSDC 2023 公开数据集。

这里的“混合频率”指不同卫星的单频、双频可用伪距可以在同一历元共同参与解算。定位使用广播星历和广播电离层改正，无需基站观测；真值轨迹用于事后评估。仓库保留了 RTKLIB 的其他功能，但论文研究对象是智能手机伪距 SPP，并非载波相位 RTK/PPP 或 Android 实时定位应用。

## 论文方法与实验方案

| 论文方案 | 含义 | 仓库内对应实验程序 |
| --- | --- | --- |
| SPP-SF | 传统 L1 单频 SPP，作为基线 | `rnx2rtkp-sf.exe` |
| SPP-HFUC | 混合频率非组合伪距模型 | `rnx2rtkp-hfuc.exe` |
| SPP-HFUC-RM | HFUC 加频率、高度角、信噪比联合定权的随机模型 | `rnx2rtkp-hfuc-rm.exe` |
| SPP-HFUC-RM-IGG3 | HFUC-RM 加自适应 IGG-III 抗差估计 | `rnx2rtkp-hfuc-rm-igg3.exe` |

**RM 对应随机模型改进，不等同于 IGG-III；`hfuc-rm` 和 `hfuc-rm-igg3` 是不同实验方案。** `Testing-Data/rtklib/` 还保留了 `rnx2rtkp-df.exe`、`rnx2rtkp-iflc.exe` 和未带方案后缀的 `rnx2rtkp.exe`，它们不作为上述四组论文主实验的替代入口。

论文方法分为三个部分：

1. **混合频率非组合模型**：各频率伪距分别建立观测方程，加入 L5 频间偏差及相对 GPS 的系统间偏差；按实际可用观测联合估计接收机位置与钟差。
2. **联合定权随机模型**：同时考虑频率差异、卫星高度角和信噪比，调整不同观测的方差。
3. **自适应 IGG-III**：依据标准化残差及四分位统计确定阈值，结合分段权函数、阈值保护和最小权值控制处理异常观测。

论文第 5.2 节比较前三种方案，第 5.3 节另行评估加入 IGG-III 后的效果。不要将不同章节、不同场景或不同基准的提升百分比合并为一个通用精度指标。完整复现时应同时记录精度与有效历元数量。

## 仓库内容

```text
Hybrid-frequency-spp/
├── RTKLIB-PRO/                         # RTKLIB 二次开发源码
│   ├── src/
│   │   ├── pntpos.c                    # 伪距模型、随机模型、IGG-III、位置估计
│   │   ├── rtklib.h                    # 状态/选项定义及公共接口
│   │   └── options.c                   # 配置项注册与读取
│   ├── app/consapp/rnx2rtkp/            # 命令行解算程序
│   │   └── msc/msc.sln                 # Visual Studio 工程
│   ├── CMakeLists.txt                  # 上游保留的 CMake 构建入口
│   ├── license.txt                     # RTKLIB 许可
│   └── readme.txt                      # 上游说明
└── Testing-Data/
    ├── config/Hybrid-frequency-spp.conf # 实验配置模板
    ├── rtklib/                         # Windows 实验程序和 crx2rnx.exe
    ├── python/                         # 日志转换、批处理、统计和绘图
    └── data/                           # 本地数据与历史处理结果，Git 忽略
```

仓库已跟踪 Python 脚本、Notebook、配置和上述 Windows 可执行文件。**原始实验数据及 `Testing-Data/data/` 下的本地结果不会随 `git clone` 下载**，需自行准备。论文原件也不随本 README 自动发布。

## 快速开始

### 1. 获取代码与准备环境

以下命令使用 Windows PowerShell，从仓库根目录开始：

```powershell
git clone https://github.com/ColourDyeing/Hybrid-frequency-spp.git
cd Hybrid-frequency-spp
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install numpy requests pandas scipy matplotlib SciencePlots jupyter
```

使用 Python 3 环境。`numpy` 用于批处理及 CSV 提取，`requests` 用于下载辅助数据；其余依赖用于 Notebook 统计绘图。当前仓库没有锁定依赖版本，正式实验应保存实际环境的版本清单。只运行已有 RINEX 文件的命令行定位时，不需要 Python 或 Qt。

### 2. 准备观测、广播星历与真值

从 GSDC 2023 数据发布渠道取得训练数据，保留原始日志与官方真值。按批处理脚本要求组织为：

```text
Testing-Data/data/train/
└── <trip>/                                  # 如 2022-02-24-18-29-us-ca-lax-o
    ├── BRDM00DLR_S_20220550000_01D_MN.rnx    # 对应观测日期的多系统广播星历
    └── <phone>/                             # 如 pixel6pro
        ├── ground_truth.csv                 # 官方真值，仅评估时需要
        └── supplemental/
            ├── gnss_log.txt                 # 原始 Android GNSS Logger 日志
            └── gnss_log.obs                 # 转换生成的 RINEX 观测文件
```

广播星历放在 **轨迹目录**，日志放在 **手机目录的 `supplemental/`**。`02_run_multi.py` 默认查找 `BRDM*MN.rnx`、`*0.2*n`、`*0.2*p`；其他命名方式需调整 `navfiles`。它只取首个匹配文件，跨日观测需自行准备覆盖全部观测时段的星历并确认选择正确。

论文选取的机型代号为 `mi8`、`pixel6pro`、`pixel7pro`、`s20ultra`、`s21ultra`、`s22ultra`。公开数据中的原始目录名可能不同，例如 `sm-g988b`、`samsungs21ultra`、`samsungs22ultra` 或 `sm-s908b`。应修改脚本的 `PHONES` 以匹配实际目录，或在工作副本中统一命名；保留原始数据。`PHONES = []` 表示自动识别全部手机目录。

广播星历也可通过 [01_get_base_data.py](Testing-Data/python/01_get_base_data.py) 下载，但运行前必须修改 `datadir` 和 `crx2rnx_bin`，它们目前是作者电脑的绝对路径。该脚本还会下载 CORS 观测并解压转换；**SPP 只需要广播星历，不必执行基站观测下载**。`crx2rnx.exe` 用于解压紧凑 RINEX，不负责将 Android 日志转为 RINEX；日志转换由 `gnsslogger_to_rnx.py` 及其配套模块完成。

### 3. 先运行一条轨迹

已有 `gnss_log.obs` 时，可直接调用实验程序。以下示例中的轨迹与手机路径必须换成自己已准备的数据：

```powershell
$trip = 'Testing-Data/data/train/2022-02-24-18-29-us-ca-lax-o'
$obs = "$trip/pixel6pro/supplemental/gnss_log.obs"
$nav = "$trip/BRDM00DLR_S_20220550000_01D_MN.rnx"
$run = 'outputs/first-run-igg3'
New-Item -ItemType Directory -Force $run | Out-Null
Copy-Item 'Testing-Data/config/Hybrid-frequency-spp.conf' "$run/solution.conf"

& './Testing-Data/rtklib/rnx2rtkp-hfuc-rm-igg3.exe' `
  -k "$run/solution.conf" `
  -o "$run/solution.pos" `
  $obs $nav
```

此示例使用带 IGG-III 的程序和模板中已开启的非组合、抗差选项。成功后检查 `solution.pos`：除 `%` 开头的头部外，应有数值解算记录。默认输出为 GPS 周、周内秒、纬度、经度、椭球高、解状态 `Q`、卫星数及精度信息。必要时加 `-x 3` 生成 trace，诊断无有效解的原因。

每次实验使用新的 `$run` 目录，保存配置副本、程序版本或哈希及输入清单。不要覆盖原始日志；旧结果可归入 `archive/`，可复现结果集中保存至 `outputs/`。

### 4. 从原始日志批量解算

编辑 [02_run_multi.py](Testing-Data/python/02_run_multi.py) 顶部参数，再运行。脚本参数在文件中设置，尚未提供对应的命令行参数接口。

| 参数 | 使用方式 |
| --- | --- |
| `DATA_SET`、`datapath` | 默认读取 `Testing-Data/data/train/`；路径相对脚本文件解析 |
| `PHONES` | 与实际手机目录名一致，或设为 `[]` |
| `binpath_rtklib` | 选择方案表中的一个可执行文件 |
| `SOL_TAG` | 与方案对应，例如 `SPP-HFUC-RM-IGG3` |
| `cfgfile_rtklib` | 指向为本次实验保存的配置副本 |
| `NEED_BASE_FILE` | SPP 保持 `False` |
| `USE_DCB_BIA` | 默认 `False`；开启时需相应 `.BIA` 文件，改变此项应单独记录实验条件 |
| `ENABLE_RTKLIB` | `True` 为转换并解算，`False` 仅转换 |
| `OVERWRITE_RINEX`、`OVERWRITE_SOL` | 仓库当前均为 `True`；复用现有文件前改为 `False` |

例如，要运行抗差方案，将 `SOL_TAG` 设为 `SPP-HFUC-RM-IGG3`，`binpath_rtklib` 设为 `../rtklib/rnx2rtkp-hfuc-rm-igg3.exe`。切换方案时同步修改标签、程序和配置，避免把结果错误归类。当前默认标签和程序为 **SPP-SF**，直接运行不会自动完成四种方案对比。

```powershell
.\.venv\Scripts\python.exe Testing-Data/python/02_run_multi.py
```

脚本将结果写为 `<trip>/<phone>/supplemental/gnss_log_<SOL_TAG>.pos`。建议在独立数据工作副本中批处理，完成后将 `.pos`、配置和统计结果归档至本次 `outputs/` 目录。即使 `OVERWRITE_SOL=False`，发生 RINEX 转换时脚本仍会重新解算；需要保留旧结果时应使用新工作副本或新的 `SOL_TAG`。

批处理当前屏蔽了解算进程的标准输出和错误信息，且未检查返回码，因此“解算全部完成”不代表每个文件都成功。遇到空文件或缺少结果时，用上一节的单文件命令查看报错。

## 配置与源码的关系

配置模板为 [Hybrid-frequency-spp.conf](Testing-Data/config/Hybrid-frequency-spp.conf)。下表描述**当前源码和模板**，不把历史二进制文件名当作编译版本证明。

| 配置项 | 模板值 | 含义 |
| --- | --- | --- |
| `pos1-posmode` | `single` | 单点定位 |
| `pos1-frequency` | `3` | 允许读取 L1/L2/L5 三个频率槽，包含 L5；不表示每颗卫星都有三频观测 |
| `pos1-posopt7` | `on` | 开启混合频率非组合处理，源码还要求 `nf >= 2` |
| `pos1-posopt8` | `on` | 开启 IGG-III 抗差处理 |
| `pos1-elmask` | `15` | 截止高度角，单位度 |
| `pos1-snrmask_L1/L2/L5` | `28 / 34 / 20` | 各高度角分段使用相同的频段 C/N0 门限，单位 dB·Hz |
| `pos1-navsys` | `45` | GPS + GLONASS + Galileo + BDS |
| `pos1-sateph`、`pos1-ionoopt` | `brdc`、`1` | 广播星历、广播电离层模型 |
| `pos1-tropopt` | `saas` | Saastamoinen 对流层模型 |
| `stats-eratio1/2/5` | `15 / 20 / 10` | 各频率槽的误差比例参数，参与观测方差计算 |
| `out-solformat`、`out-height` | `llh`、`ellipsoidal` | 大地坐标、椭球高 |
| `out-timesys`、`out-timeform` | `gpst`、`tow` | GPS 周与周内秒，匹配 CSV 提取脚本的读取方式 |

在当前源码中，`posopt7=on`、`posopt8=off` 可关闭抗差并保留非组合模型及现有随机模型；`posopt7=off` 退回首频率处理。**这些开关不能单独重建论文所有历史消融方案**：联合定权已写入 `varerr()`，未提供独立的随机模型切换开关。仓库未提供每个历史 EXE 与源码提交的一一对应构建清单；重新编译当前源码生成的是当前实现，不会自动生成四个历史算法版本。

核心代码入口：

- [pntpos.c](RTKLIB-PRO/src/pntpos.c)：`rescode()` 建立观测方程，`varerr()` 计算观测方差，`estpos()` 执行位置估计，`calc_adaptive_igg_params()` 和 `calc_classic_igg_weight()` 实现抗差阈值及权函数。
- [rtklib.h](RTKLIB-PRO/src/rtklib.h)：`POSOPT_UNCOMB`、`POSOPT_IGG3` 等选项定义。
- [options.c](RTKLIB-PRO/src/options.c)：`pos1-posopt7`、`pos1-posopt8` 等配置项映射。
- [rnx2rtkp.c](RTKLIB-PRO/app/consapp/rnx2rtkp/rnx2rtkp.c)：命令行输入与后处理入口。

## 结果提取、统计与绘图

1. 编辑并运行 `03_create_groundtruth_csv.py`，汇总各手机的 `ground_truth.csv`，生成 `Testing-Data/data/ground_truths_<DATA_SET>.csv`。
2. 在 `04_create_baseline_csv_from_pos.py` 中设置相同的 `DATA_SET`、`PHONES` 和对应 `SOL_TAG`，提取该方案的 `.pos`，生成 `locations_<SOL_TAG>_<DATA_SET>_<月_日>.csv`。每个方案分别执行一次。
3. 从 `Testing-Data/python/` 启动 Jupyter，打开 `05_Scientific_Drawing.ipynb`，检查数据目录、方案列表及筛选范围后运行。

```powershell
.\.venv\Scripts\python.exe Testing-Data/python/03_create_groundtruth_csv.py
.\.venv\Scripts\python.exe Testing-Data/python/04_create_baseline_csv_from_pos.py
cd Testing-Data/python
..\..\.venv\Scripts\python.exe -m jupyter notebook
```

Notebook 当前默认选择 `SPP-SF`、`SPP-HFUC`、`SPP-HFUC-RM`，抗差方案一行被注释；分析抗差效果时需启用 `SPP-HFUC-RM-IGG3`。`SELECTED_DATASETS=[]` 表示分析所有匹配数据，也可填入指定轨迹名。图表和报告的目录由 `DATA_DIR`、`DRAW_DIR`、`REPO_DIR` 等变量控制；归档时保存对应 CSV、筛选条件及 Notebook。

其他工具：

| 文件 | 用途 |
| --- | --- |
| [07_Observation_Quality_Analysis.ipynb](Testing-Data/python/07_Observation_Quality_Analysis.ipynb) | 观测完整性、卫星跟踪与信噪比分析，对应论文第 2 章 |
| [06_create_kml.py](Testing-Data/python/06_create_kml.py) | 将定位结果导出为 KML 轨迹 |
| [08_create_overview_kml.py](Testing-Data/python/08_create_overview_kml.py) | 生成数据集轨迹总览 |

### 统计前必须核对

- **排除填充记录**：`04_create_baseline_csv_from_pos.py` 在 `.pos` 缺失或为空时，用真值经纬度填充，并将 `Quality`、`NumSatellites` 设为 0。这些记录不是定位结果，必须排除 `Quality=0` 的行后再统计，否则可能产生虚假的低水平误差。
- **确认时间与列格式**：CSV 提取使用 GPS 周/周内秒，并包含适用于实验年代的 GPS–UTC 18 秒差；按轨迹匹配最近真值，时间差要求小于 1000 ms。改用 UTC、年月日时间格式或其他年代数据时，需要同步修改解析逻辑。脚本还固定使用 `hdrlen=25`，应与实际 `.pos` 头部核对。
- **避免混入旧实验**：真值汇总文件重跑会覆盖；同一标签的结果 CSV 在同一天重跑也会覆盖。Notebook 按文件名模式寻找 CSV，应检查它实际选中的文件，并将旧实验移入 `archive/`。
- **同时报告覆盖率**：抗差前后有效历元数可能不同，精度比较应明确使用共同历元还是各方案自身有效历元，并报告失解/丢弃情况。论文选取的子集不等于下载到的全部训练数据，需保留轨迹、手机和剔除清单，才能核对论文表格。

## 从源码构建

论文使用 Windows、Visual Studio 和 MSVC，生成 Win32 控制台程序。当前 [Visual Studio 工程](RTKLIB-PRO/app/consapp/rnx2rtkp/msc/msc.sln) 使用 **v145** 工具集和 **Windows SDK 10.0**。

1. 安装包含 C++ 桌面开发组件及对应工具集的 Visual Studio。
2. 打开 `RTKLIB-PRO/app/consapp/rnx2rtkp/msc/msc.sln`，选择 `Release | Win32` 后生成。使用其他工具集时需先重定向项目并验证编译结果。
3. 输出为 `RTKLIB-PRO/app/consapp/rnx2rtkp/msc/Release/rnx2rtkp.exe`。使用独立名称保存新程序，并记录源码提交与配置，避免覆盖历史实验 EXE。

工程启用了 `NFREQ=3`、`NEXOBS=3` 及多个星座宏。仓库同时保留了要求 CMake 3.16 及以上的上游构建入口，可用于适配其他构建环境：

```powershell
cmake -S RTKLIB-PRO -B outputs/build
cmake --build outputs/build --config Release --target rnx2rtkp
```

CMake 路线需按目标平台配置编译器及依赖；本 README 更新未验证跨平台构建。Linux/macOS 不能直接运行仓库提供的 Windows EXE，重新编译后还需修改批处理程序路径。

## 论文与复现记录

本科毕业设计：《混合频率智能手机单点定位算法设计》。

使用本项目开展实验时，建议在每个 `outputs/<实验名>/` 中保留：输入文件清单、原始数据来源、Git 提交号、可执行文件 SHA-256、配置、处理脚本参数、有效历元筛选规则及最终统计表。论文结论应能够追溯到对应数据和代码，单凭 EXE 文件名或汇总百分比无法确认复现一致性。

## 致谢与许可

本项目基于 RTKLIB 及 rtklibexplorer 分支进行二次开发；Android 日志转换使用仓库内保留的 GNSS Logger 转换模块，实验数据来源为 GSDC 2023。上游项目说明见 [RTKLIB-PRO/readme.txt](RTKLIB-PRO/readme.txt)。

RTKLIB 部分的 BSD 2-Clause 条款见 [RTKLIB-PRO/license.txt](RTKLIB-PRO/license.txt)。再分发时应保留版权与许可声明；附带工具和第三方组件遵循各自许可，GSDC 数据遵循其发布条款。
