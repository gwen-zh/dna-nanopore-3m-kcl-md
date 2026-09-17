# 本次 dA 运行协议与复现

本次结果来自 NAMD **2.14 CUDA** 多线程路径：`+p8 +devices 0 +idlepoll`，
每项 dept_gpu 的 1 张非 L40 GPU、8 CPU，不是旧 NAMD 3.0.2 的失败初始化分支。

## 实际准备链

两套 dA 都来自 DNA-only 45.11 ns 的完整 poly(dA)40 构象，刚体装孔后进行最小化、升温和平衡。
名义 3 M KCl 配盐口径、CHARMM 参数、293 K、固定盒和电场沿用原方案。
生产起点经历了孔体系准备，不能把它的 Rg 当作原 donor 的 Rg。

- 3B07 保留已完成的准备，在 eq1 第 400000 步检查点切换至 NAMD 2.14，完成 eq1、eq2、eq3。
  作业 57323197 的 production 从其 eq3 终点启动。
- 7AHL 的四个膜脂曾穿过蛋白芳环。局部重排后做 10000 步普通最小化和两段各 5000 步 MD/重启测试，
  逐次检查键长、重叠、穿环。作业 57323201 从通过检查的修复构象完成 eq1、eq2、eq3 后启动 production。
  修复诊断轨迹不计入这 5 ns。
- `inputs/production_start.coor/vel/xsc` 是真实生产起点，7AHL 的修复包含在二进制坐标中。
  **不要把初始 system.pdb 当作已修复、可直接启动生产的坐标。**

原准备安排是 50000 步最小化、0.5 ns 升温、1+1+2 ns 分级平衡，准备步长 2 fs。
production 为 **5000000 步 × 1 fs = 5 ns**，只有 production 开电场。
两孔继承的蛋白/参数兼容设置不同，不是只更换孔形状的完全同质对照。

## 生产参数

| 设置 | 7AHL | 3B07 |
|---|---|---|
| 盒尺寸（Å） | 141.7 × 141.7 × 250 | 170 × 170 × 260 |
| eField Z，NAMD 原生单位 | 0.092242 | 0.088694 |
| 标称盒跨度电势差 | 约 1 V | 约 1 V |
| 温度/压强耦合 | 293 K Langevin；固定盒、无 barostat | 同左 |
| 位置约束 | 蛋白 Cα 0.5、脂质 P 0.1、DNA 0 | 同左 |
| rigidBonds / timestep | all / 1 fs | 同左 |
| cutoff / switching / pairlist | 12 / 10 / 14 Å | 同左 |
| PME | 144 × 144 × 250 | 180 × 180 × 270 |
| DCD / restart 间隔 | 25000 步，即 25 ps | 同左 |

约束系数为 kcal mol⁻¹ Å⁻²，按 NAMD 约定，参考为归档的 restraints_prod.pdb。
约 1 V 是整个周期盒的名义电势差，不是测得的局部膜电压；3 M 是构建设置，不是已测量的孔内/水相浓度。

## 复现配置

每个结果目录的 production.namd 复现生产设置。inputs 是完整体系，外层 visualization.psf 仅供显示。

```bash
# 在新复制的结果目录中执行，不要覆盖集群已有模拟。
gzip -dk inputs/system.psf.gz inputs/system.pdb.gz inputs/restraints_prod.pdb.gz
gzip -dk inputs/production_start.coor.gz inputs/production_start.vel.gz
export CHARMM_TOPPAR=/path/to/charmm/toppar
namd2 +p8 +devices 0 +idlepoll production.namd > rerun.log
```

只在已有 GPU 分配的环境中运行；本集群仍须使用 dept_gpu、排除 L40、每项 1 GPU。
这里没有自动提交脚本。CHARMM 参数哈希在 forcefield_sha256.json，NAMD/参数/水离子文件需匹配。
不同硬件和并行归约会使数值轨迹逐渐分离，不承诺位级重现。

7AHL checkpoint/prod5.* 对应 5000000 步终点。
3B07 预览只提供生产起点和截至 4825000 步的数值记录/显示帧，不伪造未完成的闭合终点。
