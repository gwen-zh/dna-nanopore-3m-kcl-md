# DNA 纳米孔模拟 protocol

原始构建版本日期：2026-09-11。下文保留原构建流程，不代表当前软件与作业状态。
2026-09-17 发布的 dA 轨迹实际使用 NAMD 2.14 CUDA；7AHL 经过额外膜修复和检查。
当前结果、视频及可复现的生产起点见[结果目录](results/2026-09-17/README.md)和[运行协议](results/2026-09-17/PROTOCOL.md)。
不要跳过膜几何检查，直接用下面的历史提交示例重启已修复的任务。
运行状态另见 [PORE_JOB_STATUS.md](../PORE_JOB_STATUS.md)。

## 1. 目的与体系

比较 poly(dA)40 和 poly(dT)40 在高盐条件下接近、进入纳米孔时的构象及
DNA–蛋白相互作用。四套体系均含一条 40 nt ssDNA、蛋白孔、POPC 膜、
显式水和名义浓度 3 M 的 KCl。

7AHL 是金黄色葡萄球菌 α-溶血素七聚体；3B07 是 γ-溶血素八聚体，
不是同一种蛋白的两个 PDB 版本。结构来源：[7AHL](https://www.rcsb.org/structure/7AHL)、
[3B07](https://www.rcsb.org/structure/3B07)。

本流程从仓库中的已处理蛋白/膜/水模板开始。模板不含 DNA 或离子；
继承其残基、末端、质子化状态和膜构型，不在本流程中重建晶体结构缺失片段。

| 体系目录（位于 nanopore/） | 原子总数 | DNA 原子数 | K⁺ | Cl⁻ |
|---|---:|---:|---:|---:|
| polyda40-3m-kcl-7ahl/run | 477,938 | 1,279 | 7,661 | 7,636 |
| polydt40-3m-kcl-7ahl/run | 477,964 | 1,279 | 7,662 | 7,637 |
| polyda40-3m-kcl-3b07/run | 720,355 | 1,279 | 12,091 | 12,132 |
| polydt40-3m-kcl-3b07/run | 720,372 | 1,279 | 12,092 | 12,133 |

数量来自 2026-09-11 构建输出。K⁺/Cl⁻ 数量差用于中和其余体系电荷。

## 2. DNA 构象来源与完整性

DNA 来源于对应游离 ssDNA 的连续 50 ns、3 M KCl 轨迹，而不是重新拉直的链。
本次使用原始 `prod20to50.dcd` 的以下帧，帧号从 0 开始：

| DNA | 原始 DCD 帧号 | 大致轨迹时刻 | 可用于建模的文件 |
|---|---:|---|---|
| poly(dA)40 | 2510 | 约 45.1 ns | [polyda40_donor_intact.pdb](common/polyda40_donor_intact.pdb) |
| poly(dT)40 | 2230 | 约 42.3 ns | [polydt40_donor_intact.pdb](common/polydt40_donor_intact.pdb) |

[extract_intact_donor.tcl](common/extract_intact_donor.tcl) 从原始帧读取坐标，
依据 PSF 的键连接执行 `pbc join fragment -bondlist`，再检查每条 DNA 键。
重建后 DNA 最大键长约 1.645 Å（dA）和 1.638 Å（dT）。

旧的逐原子、跨帧 unwrap 导出曾造成少量约 70 Å 的异常键。
因此，不用 `ssdna/` 中旧的 representative PDB 直接构建过孔体系。
当前沿用之前选定的时间帧，但不声称它们是修正分析后的 Rg 中位数构象。
原始轨迹未改动；旧 Rg、堆积和碱基暴露结果需另外复核。

## 3. 体系构建

由 [build_pore.tcl](common/build_pore.tcl) 和 [place_dna.py](common/place_dna.py) 执行：

1. 读取 `7ahl-template/template_noions.*` 或 `3b07-template/template_noions.*`。
   对仓库内这些清理后的模板，设置 `REMOVE_OLD_DNA=0`。
2. 将 DNA 的 1→40 残基 C1′ 端到端方向对齐 +Z；5′ 端位于孔的 cis 侧，
   朝向孔入口。只进行刚体旋转/平移，不拉伸或重塑 DNA。
3. 分别搜索绕 Z 轴的方位角和放置高度。7AHL 的起始 C1′ 高度为 88 Å，
   3B07 为 90 Å；搜索步长为 15° 和 2 Å。
4. 要求 DNA 离模拟盒边界至少 10 Å；拒绝与蛋白或膜小于 1.2 Å 的直接碰撞。
   放置程序同时统计 2.2 Å 的蛋白近接触。这些是初始几何检查，不代替平衡。
5. 使用 `top_all36_na.rtf`，构建 DNA segment `AN1`，应用 5TER/3TER、
   DEO5/DEOX 补丁，重新生成角和二面角并补齐坐标。
6. 删除距新 DNA 2.4 Å 内的整分子旧水，沿 Z 方向扩展溶剂；新增水使用
   `NW` segment 前缀，避免与模板水段重名。
7. 中和并加入 KCl：`-sc 3.0 -cation POT -anion CLA -from 5.0 -between 3.5`。

[fast_autoionize.tcl](common/fast_autoionize.tcl) 保留已安装 autoionize 的离子
计数、拓扑及替换水分子流程，只替换逐个选位的搜索。
[fast_ion_positions.py](common/fast_ion_positions.py) 使用带种子的随机候选顺序
和周期 KD-tree；离子间距至少 3.5 Å，离子与非水原子间距至少 5 Å。
随机种子和距离验证保存在各体系的 `ion_placement_validation.json` 中。
这里记录的是候选位置选择的种子；后续 autoionize 对 K⁺/Cl⁻ 的位置分配仍使用
其随机打乱步骤，未单独固定种子，因此重新构建不保证离子坐标逐位相同。

本版本 autoionize 的盐对计数为 `round(0.0187 × 3.0 × 加盐前水分子数)`，
再加中和离子。因此“3 M”是构建程序的名义设置，不是已经测得的孔内或
体相浓度；定量解释时需要核查平衡后的体相离子密度。

构建后和 GPU 启动前均运行 [validate_bonds.py](common/validate_bonds.py)：
验证 PSF/PDB 原子数、有限坐标、完整 DNA，以及周期最小镜像键长。
DNA 键必须在 0.5–2.2 Å 内；继承膜模板的少量拉长键由后续最小化松弛，
全体系检查上限为 3.5 Å。通过此检查不等同于已经平衡。

## 4. 力场、模拟盒与积分设置

软件：NAMD 3.0.2，VMD 1.9.2，Python 3（NumPy、SciPy），Slurm。

参数读取顺序见 [common.namd](common/common.namd)：

- 7AHL 先读 `par_all22_prot.prm`，提供继承 PSF 中 HB 等旧原子类型的参数。
- 随后读取 `par_all36m_prot.prm`、`par_all36_lipid.prm`、`par_all36_na.prm`。
- 水/离子读取各体系的 `water_ions_namd.prm`，其中包含 CHARMM TIP3P 水参数。
- 3B07 不启用上述 CHARMM22 补充参数。

这保留了当前模板的参数兼容性；不能把两种孔的设置笼统描述为完全统一的
CHARMM36m 蛋白模型。跨孔比较需考虑该差异。参数目录 `/opt/charmm/toppar`
和 NAMD module 名称是集群设置；在其他机器运行时需对应修改。

| 设置 | 7AHL | 3B07 |
|---|---|---|
| 盒长 X × Y × Z（Å） | 141.7 × 141.7 × 250 | 170 × 170 × 260 |
| 盒中心（Å） | (0, 0, 65) | (0, 0, 60) |
| 盒范围（Å） | X/Y ±70.85；Z −60 至 190 | X/Y ±85；Z −70 至 190 |
| PME 网格 | 144 × 144 × 250 | 180 × 180 × 270 |
| 温度 | 293 K | 293 K |

三维周期边界，`wrapAll on`、`wrapNearest on`；PME 静电；切换距离 10 Å，
截断 12 Å，pair list 14 Å；`exclude scaled1-4`、`1-4scaling 1.0`。
`rigidBonds all` 约束水及含氢共价键，而不是固定全部化学键。
其含义见 [NAMD 约束文档](https://www.ks.uiuc.edu/Research/namd/3.0.2/ug/node29.html)。

准备阶段步长 2 fs，production 为 1 fs；`nonbondedFreq 1`、
`fullElectFrequency 2`、`stepsPerCycle 20`。
模拟盒固定，没有 Langevin piston/barostat。升温后使用 Langevin 控温，
目标 293 K，阻尼 1 ps⁻¹，`langevinHydrogen off`。
不要把本流程写成 1 atm 的 NPT 平衡；游离 ssDNA 的压力耦合流程与此不同。

## 5. 最小化、升温、平衡与 production

| 顺序 | 输入 | 积分步长 | 步数 | 物理时间 | 温度/外场 |
|---|---|---|---:|---|---|
| 1 | [min.namd](common/min.namd) | 最小化，不按 MD 时间计 | 50,000 | 不计 ns | 无电场 |
| 2 | [heat.namd](common/heat.namd) | 2 fs | 250,000 | 0.5 ns | 50 → 293 K；无电场 |
| 3 | [eq1.namd](common/eq1.namd) | 2 fs | 500,000 | 1 ns | 293 K；无电场 |
| 4 | [eq2.namd](common/eq2.namd) | 2 fs | 500,000 | 1 ns | 293 K；无电场 |
| 5 | [eq3.namd](common/eq3.namd) | 2 fs | 1,000,000 | 2 ns | 293 K；无电场 |
| 6 | [prod5.namd](common/prod5.namd) | 1 fs | 5,000,000 | 5 ns | 293 K；开启电场 |

总 MD 时间为 9.5 ns，其中升温 0.5 ns、平衡 4 ns、production 5 ns。
50,000 步能量最小化单独计数，不能换算成 0.1 ns。

升温采用速度重分配：每 1,000 步（2 ps）增加 1 K，达到 293 K 后保持到
0.5 ns 结束。每阶段读取上一阶段的坐标、速度（适用时）及盒信息。

### 位置约束

[make_restraints.tcl](common/make_restraints.tcl) 将系数写入 PDB 的 B 列，
参考坐标是初始构建坐标。下表为乘过 `constraintScaling` 的系数，
单位 kcal mol⁻¹ Å⁻²，遵循 NAMD 的 `k × 位移²` 约定。
[NAMD 位置约束定义](https://www.ks.uiuc.edu/Research/namd/3.0.2/ug/node29.html)。

| 阶段 | 蛋白 Cα | DNA backbone | 脂质 P | 其余原子位置约束 |
|---|---:|---:|---:|---|
| min / heat / eq1 | 1.0 | 1.0 | 0.5 | 无 |
| eq2 | 0.2 | 0.2 | 0.1 | 无 |
| eq3 / production | 0.5 | 0 | 0.1 | 无 |

eq3 的“DNA 自由”不表示全体系无约束：蛋白和膜仍保留位置约束。
eq3 的蛋白系数由 eq2 的 0.2 调整为 0.5；这就是当前输入，不是单调减弱。
已选中 DNA backbone 317 个原子；7AHL 的蛋白 Cα / 脂质 P 为 2,051 / 497，
3B07 为 2,216 / 639。

## 6. 电场与电势差

仅 production 开启沿 +Z 的均匀电场：

```tcl
eFieldOn on
eField 0.0 0.0 $efield
```

`eFieldNormalized` 未开启（默认 no）。NAMD 此时的场强单位为
kcal mol⁻¹ Å⁻¹ e⁻¹，见 [电场文档](https://www.ks.uiuc.edu/Research/namd/3.0.2/ug/node42.html)。
使用当前输入的换算常数，盒子两端的名义电势差绝对值为：

`|ΔV| ≈ |E_z| × L_z / 23.0605`（V）

| 孔 | eField Z | Lz（Å） | 名义盒跨度电势差 |
|---|---:|---:|---:|
| 7AHL | 0.092242 | 250 | 约 1.0 V |
| 3B07 | 0.088694 | 260 | 约 1.0 V |

这不是已经从电势剖面测得的局部膜电压。DNA 带负电，电场方向的选择使其
直接电场力指向 −Z；实际迁移同时受到离子、溶剂及孔相互作用影响。
没有对 DNA 加恒速拉伸或 steered-MD 拉力。
该电场设置下的压力输出还含电场相关贡献，不用它直接解释平衡压力。

## 7. 集群执行与输出

所有 NAMD 准备及 production 作业只使用 `dept_gpu`，每套 1 GPU、
8 CPU cores、16 GB RAM、普通优先级。准备作业时限 2 天，production 3 天；
这是调度上限，不是实测完成时间，也不提供自动跨时限续跑。
构建及分析使用 `dept_cpu`；四个体系可并行，但同一体系的准备和 production
通过 `afterok` 串联，最多同时使用 4 张 GPU。

以下是从仓库根目录提交 dA40–7AHL 的示例。仅用于新运行目录；
不要在已有活跃任务或已有输出的目录重复提交：

```bash
export REPO_DIR="$PWD"
export SYSTEM_DIR="$REPO_DIR/nanopore/polyda40-3m-kcl-7ahl/run"
export PORE_PSF="$REPO_DIR/nanopore/7ahl-template/template_noions.psf"
export PORE_PDB="$REPO_DIR/nanopore/7ahl-template/template_noions.pdb"
export DONOR_PDB="$REPO_DIR/nanopore/common/polyda40_donor_intact.pdb"
export DNA_BASE=ADE REMOVE_OLD_DNA=0 TARGET_Z=88
export MIN_X=-70.85 MIN_Y=-70.85 MIN_Z=-60
export MAX_X=70.85 MAX_Y=70.85 MAX_Z=190
export SYSTEM_LABEL=polyda40_7ahl ANALYSIS_PREFIX=polyda40_7ahl

build=$(sbatch --parsable --export=ALL \
  --output="$SYSTEM_DIR/build.%j.out" --error="$SYSTEM_DIR/build.%j.err" \
  nanopore/common/build_job.sbatch)
prep=$(sbatch --parsable --export=ALL --dependency="afterok:$build" \
  --output="$SYSTEM_DIR/prep.%j.out" --error="$SYSTEM_DIR/prep.%j.err" \
  nanopore/common/prep_job.sbatch)
prod=$(sbatch --parsable --export=ALL --dependency="afterok:$prep" \
  --output="$SYSTEM_DIR/prod.%j.out" --error="$SYSTEM_DIR/prod.%j.err" \
  nanopore/common/prod_job.sbatch)
```

dT 对应 `DNA_BASE=THY` 和 dT donor；3B07 对应其模板、`TARGET_Z=90`、
X/Y 边界 ±85 Å、Z 边界 −70/190 Å。同时替换目录、体系标签和分析前缀。
各自的 `system_settings.namd` 已给出对应盒长、网格、电场及旧参数开关。
不要把 7AHL 的盒设置用于 3B07。

准备阶段的 DCD、XST 和 checkpoint 每 25,000 步输出，即 50 ps；
production 每 25,000 步输出，即 25 ps。能量/计时分别每 1,000 / 5,000 步输出。
最小化日志可能额外逐次输出能量，不能用其输出数量推算 MD 时间。
原始 DCD 不存入普通 Git；DNA-only 已有 SHA-256 清单，过孔 DCD 待生成后归档。

## 8. 质量控制和结果解释

- 开始前检查完整 DNA、PSF/PDB 对应、离子计数及初始接触，确认所加载的参数
  能覆盖 PSF 原子类型，不能只检查作业是否显示 RUNNING。
- 每阶段检查程序正常退出、没有 NaN/Inf、FATAL ERROR、异常温度或坏键。
  当前准备脚本在阶段结束时检查日志；它不是持续运行的故障监控服务。
- 进入 production 前仍应审核温度、能量、膜/孔结构和 DNA 构象的稳定性。
  当前脚本按预设时长和退出状态自动串联，并不自动判断物理上的平衡充分性。
- 计划分析 Rg、端到端距离、5′ 端/质心 Z 位置、孔内占据、DNA–蛋白接触及
  碱基 SASA、相邻碱基堆积、糖苷键二面角。它们是不同指标，不能把任意一个
  单独称为已证明的“碱基外翻”。
- 仓库中的分析脚本与已有 ssDNA 结果仍需周期镜像、孔参考坐标和时间轴复核。
  特别是现有 pore 分析对 DNA 做跨帧 unwrap，而蛋白未同步成像；碱基分析还需
  检查完整分子坐标、糖苷二面角原子顺序及 SASA 的遮挡选择。
  **本次未重算分析；暂不将其原始输出作为已经验证的相互作用结论。**
- 5 ns 是此次 production 长度，不保证能够观察到完整过孔。需要区分未捕获、
  已捕获、部分进入、停滞和完整穿越；单条轨迹不足以估计可靠的通过率或误差。

本次仓库命名调整仅作用于仓库副本。正在运行的集群目录和任务输入路径未改名，
以免破坏阶段衔接。旧日志的工作目录头已做归档路径标准化，数值记录未改动。
