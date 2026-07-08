# TDE Calculator

**阈值位移能 + DPA 计算工具包 — 面向单质金属、合金材料**

集成 **LAMMPS 分子动力学**（Byggmästar knock-on atom 方法）进行阈值位移能（Threshold displacement energy，TDE）计算、**收敛性可视化**、以及 **DPA 对应具体级联次数的计算**（NRT / ARC-DPA / CB-DPA 解析模型）。

![工作流示意](figures/workflow.svg)

[![Python 3.7+](https://img.shields.io/badge/python-3.7%2B-blue)](https://python.org)[![LAMMPS](https://img.shields.io/badge/LAMMPS-2023-orange)](https://lammps.sandia.gov)[![License: MIT](https://img.shields.io/badge/license-MIT-green)]()

---

## 快速开始

```bash
# 1. 启动程序
python3 main.py

# 2. 配置环境路径（菜单选项 9）
#    → 设置 mpirun 路径和 LAMMPS 可执行文件路径（例如：lmp_mpi）

# 3. 清理旧结果（菜单选项 4）— 可选

# 4. 运行 TDE 计算（菜单选项 1）
#    → 输入元素 → 选择势函数 → 选择自动生成 (auto) / 人为给定 (custom) 构型
#    → 设置参数 → 确认 → LAMMPS 自动运行

# 5. 绘制 Ed 收敛图（菜单选项 2）
#    → 选择结果目录 → 保存 plot_Ed.png

# 6. 计算 DPA 对应的级联次数（菜单选项 3）
#    → 选择模型：NRT / ARC-DPA / CB-DPA

# 7. 退出（菜单选项 0）
```

---

## 目录结构

```
TDE_calculator/
├── main.py                     菜单驱动入口
├── ReadMe_zh.md                本文件
├── env.conf                    环境配置（mpirun 与 LMP 路径）
├── pot_files/                  放置 .zbl 势函数文件
├── custom_data_files/          "custom"构型来源的数据文件
├── results/                    Ed 计算结果（程序生成）
├── Ed_summary.txt              统计汇总（程序生成）
│
├── Ed_calc/                    TDE 计算（LAMMPS MD）
│   ├── run_Ed.sh               主控脚本（多构型、并行）
│   ├── alloy.conf              合金参数（由 main.py 生成）
│   ├── in.ed.equilibrate.lmp   LAMMPS 构型平衡
│   ├── in.ed.recoil.lmp        LAMMPS 撞击模拟
│   ├── check_defects.py        OVITO Wigner-Seitz 缺陷分析
│   ├── collect_results.py      元素分辨统计
│   └── generate_directions.py  随机单位矢量生成
│
├── Ed_plot/                    绘图模板
│   └── plot_Ed.py              收敛图（1×2 子图）
│
└── DPA_calc/                   DPA 解析计算器
    └── dpa_calculator.py       NRT / ARC-DPA / CB-DPA + Robinson 分配函数
```

---

## 菜单选项

```
============================================================
  TDE Calculator v0.3.0
  Radiation Damage Calculation Toolkit
============================================================
  1. TDE 计算     （Ed_calc, LAMMPS MD）
  2. 绘制 Ed 图   （Ed_plot, 收敛性）
  3. DPA 计算器   （NRT / ARC-DPA / CB-DPA）
  4. 清理         （删除所有生成文件）
  9. 环境设置     （设置 mpirun 与 LAMMPS 路径）
  0. 退出
============================================================
```

---

### 1. TDE 计算

**5 步交互设置**配置材料、势函数和模拟参数，然后自动启动 LAMMPS。

#### 第 1 步：元素
输入元素符号（空格分隔）。质量由内置周期表自动查询。

```
Elements [Hf Nb Zr Ti Ta]: W Ta Cr V
→ Auto-looked up: 183.840 180.948 51.996 50.942
```

#### 第 2 步：势函数
从 `pot_files/` 中选择势函数文件（EAM/alloy 格式），输入序号。

#### 第 3 步：构型来源
- `auto` = 程序自动生成随机固溶体（BCC 或 FCC 晶格）。
  单元素→纯金属；多元素→等原子比随机合金。
- `custom` = 用户提供已构建好的数据文件（`custom_data_files/data.1`、`.2`、…）。
  （`auto` 模式仅支持 BCC 和 FCC；其他晶系请使用 `custom`。）

还需输入晶体结构（`bcc`/`fcc`）和晶格常数（Å）。

> 💡**注意：**当前仅支持**正交晶系**，非正交晶系可通过atomsk等软件切割出正交超胞后，将data文件放入custom_data_files文件夹

#### 第 4 步：模拟参数

| 参数 | 默认值 | 含义 |
|------|--------|------|
| **NCORE** | 16 | 每个 `mpirun` 任务的 CPU 核数 |
| **NJOB** | 4 | 并发晶向数（总核数 = NJOB × NCORE） |
| **NCONFIG** | 5 | 独立随机构型数 |
| **NDIR_PER_CONFIG** | 100 | 每构型采样的晶向数 |
| **SIMTIME** | 6.0 ps | 每次能量尝试的模拟时长 |
| **EMIN** | 10 eV | 二分搜索能量下限 |
| **EMAX** | 180 eV | 二分搜索能量上限 |
| **ESTEP** | 1 eV | 能量分辨率（测试能量量化到 ESTEP 的倍数） |
| **TEMP** | 40 K | 模拟温度（低温减少热噪声） |
| **DEFECT_METHOD** | ovito | 缺陷检测算法（`displace` / `cna` / `ptm` / `ovito`） |

#### 第 5 步：确认
显示汇总信息，确认后运行 `run_Ed.sh`。

**输出**：`results/config_*/Ed_direction_*.txt`、`Ed_summary.txt`。

---

### 2. 绘制 Ed 收敛图

扫描 `results/` 中的 `config_*/Ed_direction_*.txt`，在选中目录中运行 `plot_Ed.py`。生成 1×2 子图：左为总体累积均值，右为分元素累积均值。

**输出**：`plot_Ed.png` (dpi=300)。

### 3. DPA 计算器

以交互模式（`-i`）启动 `dpa_calculator.py`。支持 NRT、ARC-DPA、CB-DPA 三种模型，以及 direct 和 Robinson 两种能量分区方式。

### 4. 清理

删除所有生成文件（`results/`、`logs/`、`dump_files/`、`Ed_summary.txt`、`rss_substitute.lmp` 等）。保留源代码、`pot_files/` 和 `custom_data_files/`。

### 9. 环境设置

设置 `mpirun` 和 LAMMPS 可执行文件路径。保存到 `env.conf`。

---

## 自定义模拟

### 修改模型尺寸

`auto` 模式下的原子数由 `Ed_calc/in.ed.equilibrate.lmp` 第 30 行的 `Nsize` 控制：

```lammps
variable Nsize equal 16    # 默认 BCC: 16³×2 = 8192 原子；FCC: 16³×4 = 16384 原子
```

要创建更大的模型，修改此行：

```lammps
variable Nsize equal 20    # 改为 BCC: 20³×2 = 16000 原子
```

### 修改势函数类型（pair_style）

当前默认的势函数类型为 **EAM/alloy**。如需切换为其他类型（如 MEAM），需同时修改两个文件：

1. `Ed_calc/in.ed.equilibrate.lmp`（第 47–48 行）：

```lammps
# 默认 eam/alloy — 更改 pair_style 并相应调整 pair_coeff：
pair_style      eam/alloy
pair_coeff      * * ${POTENTIAL} ${ELEMENTS}
```

2. `Ed_calc/in.ed.recoil.lmp`（第 34–35 行）：

```lammps
# 此处必须做相同的更改：
pair_style      eam/alloy
pair_coeff      * * ${POTENTIAL} ${ELEMENTS}
```

---

## 缺陷检测方法

通过 `DEFECT_METHOD` 参数选择四种方法之一。每种方法判定 PKA 能量是否产生了稳定的 Frenkel 对。

| 方法 | 关键参数 | 默认值 | TDE 数值 | 说明 |
|------|----------|--------|---------|------|
| **displace** | 位移阈值 | 1.0 Å | **偏低** | 过度灵敏；会将替换碰撞序列误判为缺陷 |
| **CNA** | 截断半径 | 3.5 Å | — | **不推荐**用于多元素合金 |
| **PTM** | RMSD 阈值 | 0.1 | 与 OVITO 相近 | 对 BCC/FCC 可靠；抗噪声能力强 |
| **OVITO/WS** | Wigner-Seitz方法直接计数 Frenkel pairs | — | **标准** | 最精确；需要 OVITO（`ovitos`） |

### 详细对比

- **displace**（1.0 Å 阈值，`in.ed.recoil.lmp:120`）：统计所有位移超过 1.0 Å 的原子。这是最灵敏的方法——会将**替换碰撞序列**（原子移动到相邻晶格位但未产生空位/间隙）也判定为 DEFECT，导致 TDE **系统性地偏低**。适用于快速预览。

- **CNA**（3.5 Å 截断半径，`in.ed.recoil.lmp:130`）：共近邻分析识别局部结构不匹配的原子。**不推荐**用于多元素随机固溶体，因为化学无序本身可能导致 CNA 将完好的原子误判为缺陷。

- **PTM**（0.1 RMSD 阈值，`in.ed.recoil.lmp:138`）：多面体模板匹配将理想的 BCC/FCC 模板拟合到每个原子的近邻壳层。对于 BCC 和 FCC 合金，**PTM 与 OVITO/WS 通常给出相同的 Ed**，因为两种方法检测的是同一物理事件——Frenkel 对的形成。

- **OVITO/WS**（`check_defects.py`）：Wigner-Seitz 分析将最终原子构型与参考构型对比。当 `interstitial_count > 0` 或 `vacancy_count > 0` 时判定为 DEFECT。这是**物理上最严格的方法**，直接计数空位和间隙原子数。需要 `ovitos` 命令（OVITO Python API）。

详细对比见 [`methods.md`](methods.md)。

---

## DPA 模型

DPA 计算器（`DPA_calc/dpa_calculator.py`）实现了三种解析模型：

| 模型 | 参考文献 | 效率函数 ξ(Ea) | 参数 |
|------|----------|-----------------|------|
| **NRT** | Norgett, Robinson & Torrens (1975) | $ξ = 1$ | Ea, Ed |
| **ARC-DPA** | Nordlund et al. (2018) | $ξ = (1−c)·(0.8Ea/2Ed)^b + c$ | b_arc, c_arc（MD 拟合） |
| **CB-DPA** | Chen, Bernard et al. (2020) | $ξ = (2Ed/0.8)/(2Ed/0.8+β·Ea) + β$ | β = Z/(1.5·A) |

损伤能量 Ea 由 PKA 能量通过 **direct**（Ea = E_PKA）或 **NRT 论文中使用的 Robinson 拟合（基于 Lindhard 离子阻止理论）** 方法计算。

**用法**：

```bash
# 交互模式
python3 DPA_calc/dpa_calculator.py -i

# 对比所有模型
python3 DPA_calc/dpa_calculator.py --compare --b_arc -0.568 --c_arc 0.286
```

---

## 算例展示

HfNbZrTiTa 的 TDE 计算结果（RSS，5 构型 × 100 入射方向 = 500 总计，OVITO/WS 方法）为例，势函数采用五元EAM势函数，TDE计算**采用默认设置**。500样本计算共用时约2.5小时，用 `功能2` 绘制 TDE $E_d$ 收敛图得到：

<img src="figures/plot_Ed.png" alt="plot_Ed" style="zoom:50%;" />

计算所得TDE ~ 59.4 eV，与文献报道值 62 eV十分接近 (10.1016/j.ijplas.2024.104155, 10.1016/j.matlet.2020.127832)。

---

## 依赖项

| 组件 | 要求 |
|------|------|
| **Ed_calc**（TDE） | LAMMPS、MPI （`mpirun`）、OVITO（`ovitos`，版本>=3.4.4，可选） |
| **Ed_plot**（绘图） | Python 3.7+、numpy、matplotlib、seaborn |
| **DPA_calc**（DPA） | Python 3.7+（标准库） |
| **main.py**（主程序） | Python 3.7+（标准库） |

---

## 引用

如果您在研究中使用了本工具包，请引用以下参考文献：

- **TDE Calculator（本软件）**：10.1016/j.ijplas.2026.104626, 10.48550/arXiv.2606.26019.
- **TDE 计算方法**：Byggmästar, J., Djurabekova, F., & Nordlund, K. (2024). *Phys. Rev. Materials*, 8, 115406.
- **NRT 模型**：Norgett, M.P., Robinson, M.T., & Torrens, I.M. (1975). *Nucl. Eng. Des.*, 33, 50-54.
- **ARC-DPA**：Nordlund, K., et al. (2018). *Nat. Commun.*, 9, 1084.
- **CB-DPA**：Chen, S., et al. (2020). *EPJ Web Conf.*, 239, 08003.

---

## 参考文献

| 方法 | 参考文献 |
|------|----------|
| NRT 模型 | Norgett, Robinson & Torrens (1975). *Nucl. Eng. Des.*, 33, 50-54. |
| ARC-DPA | Nordlund, Zinkle, Sand et al. (2018). *Nat. Commun.*, 9, 1084. |
| CB-DPA | Chen, Bernard, Tommasi, De Saint Jean (2020). *EPJ Web Conf.*, 239, 08003. |
| TDE（knock-on atom） | Byggmästar, Djurabekova & Nordlund (2024). *Phys. Rev. Materials*, 8, 115406. |
| PTM | Larsen, P.M., et al. (2016). *Modelling Simul. Mater. Sci. Eng.*, 24, 055007. |
| CNA | Honeycutt, J.D. & Andersen, H.C. (1987). *J. Phys. Chem.*, 91, 4950. |
