# 分析与视频方法

作业 57323201（dA–7AHL）及 57323197（dA–3B07）的 production 步长 1 fs，
DCD 每 25000 步保存，即 25 ps。首帧是 0.025 ns，不是 0 ns。
额外读取实际 eq3 起点作为 0 ns，视频分别有 201 和 194 帧。
3B07 只读取开始分析时已完整写出的帧，不追读后续帧。

## 周期处理与指标

DNA 和各蛋白链通过 PSF 共价连接逐帧重建。DNA 根据相邻帧质心选取连续周期像。
蛋白 Cα 用正旋转 Kabsch 对齐到 production 起点，同一变换用于所有显示原子。
接触计算仍使用原始瞬时盒的周期最短距离。显示 DCD 不附加未旋转晶胞，不能作周期重启数据。

- Rg：DNA 重原子的质量加权回转半径。
- 端距：第 1、40 核苷酸 C1′ 距离；5′轴向推进为起点 z 减当前 z，正值朝 −Z。
- cis/trans 参考平面：起点全蛋白 Cα 的最大/最小 z，只是几何参照，不是孔径分析得到的狭窄口。
- `c1_below_trans_CA_plane`：位于 trans 平面以下的 C1′ 数量。本次没有达到 40 的帧。
  单独的轴向位置不等同于位于孔腔内部，DNA 质心移动也不等于完整穿孔。
- 接触：DNA–蛋白重原子周期距离 ≤3.5 Å；报告原子对数、DNA 核苷酸数及蛋白残基帧占比。
  蛋白残基必须用链+编号标识；占比的分母只包含 DCD 帧，不包括 0 ns。
- `stacked_adjacent_pairs_ring_centers`：ADE 的 9 个嘌呤环重原子几何中心；相邻中心距 ≤5.5 Å、
  法向绝对余弦 ≥cos45°、两平面高度各 2–4.5 Å、侧向偏移各 ≤3.5 Å。
  该辅助列使用环中心定义，不直接等同于旧 ssDNA 分析中不同中心定义的列，也不是外翻/SASA 指标。
- DNA 镜像间距：完整 DNA 重原子到 26 个非零邻盒镜像的最短距离；先检查每轴 DNA 跨度小于盒长。
  不混入零平移内原子距离。无短程镜像接触不排除所有长程有限尺寸效应。

没有时间平滑、按预期卷曲程度挑帧或人工改变 DNA。此分析不是离子电流、自由能、穿孔速率或独立重复统计。

## 视频验证与复现

VMD 1.9.2 / TachyonInternal 用 CPU 渲染实际帧，参考 [VMD render 文档](https://www.ks.uiuc.edu/Research/vmd/current/ug/node147.html)。
前侧蛋白透明，水/盐隐藏，膜只画 P 原子；接触着色是显示辅助，数量来自周期距离计算。
VMD 独立计算首末帧 Rg、端距及接触数并与 Python 比较，记录在 `render_checks_*.csv`。
MP4 编码后完整解码核对帧数，并检查分子画面非空、没有冻结为同一张图。
10 fps 是播放设置，不能当作真实时间尺度；没有坐标插值。

安装 NumPy、SciPy、pandas、matplotlib、Pillow、imageio-ffmpeg 和 VMD，然后：

```bash
gzip -dk visualization.dcd.gz
VMD_BIN=/path/to/vmd python /path/to/nanopore/analysis/movie.py /path/to/result-directory
```

默认 VMD 路径为 `/opt/bin/vmd`，字体为 Linux DejaVu Sans。完整数值分析需原始全体系 PSF、DCD 和生产起点；
来源/哈希见 `summary.json`。轻量轨迹不能用于溶剂/离子分析或 MD 重启。
