# matchfamily — 匹配与加权（jamovi 模块）

当所有混杂因素均可观测时估计二值处理的效应：倾向得分与重叠、协变量平衡、回归调整、加权、双重稳健估计、匹配、熵平衡和对未观测混杂的敏感性。 英文完整说明见 [README.md](README.md)。

## 分析（jamovi 菜单：Matching）

| 分析 | 内容 |
|---|---|
| **倾向得分与重叠** | logit、probit 或 CBPS，共同支撑与修剪 |
| **协变量平衡** | 标准化均值差、方差比、KS 统计量，Love 图 |
| **加权与双重稳健估计** | RA、IPW、AIPW、IPWRA，ATE/ATT/ATC/ATO |
| **匹配估计量** | k 近邻、多种距离、精确匹配、卡尺、偏差调整 |
| **熵平衡** | 精确匹配均值（可选方差与偏度） |
| **比较估计量** | 各方法并列与森林图 |
| **对未观测混杂的敏感性** | E 值与 Rosenbaum 界 |

每个分析的输出部分都有"方法说明"和"数据与设置指南"。

引擎为基础 R 编写，与 Stata（`teffects`、`ebalance`）及 R 包（`WeightIt`、`cobalt`、`Matching`、`EValue`、`rbounds`、`PSweight`、`CBPS`）对照，估计量与标准误一致到 1e-8 量级（IPWRA 在 LaLonde 数据上标准误相差 1–2%，Stata 自身报告未收敛）。示例数据：`lalonde`、`sim_confounded.csv`、`sim_overlap.csv`、`sim_hidden.csv`。

## 安装
- **jamovi 模块库**：在 jamovi 右上角打开模块按钮，进入 *jamovi library*，搜索并安装（模块上架后）。
- **侧载**：下载 `.jmo` 文件，在同一菜单选 *Install from file*。`.jmo` 只适用于构建它的操作系统、处理器类型和 jamovi 版本。
- **从源码**：需要 R 和 `jmvtools`：`jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")`（路径按本机 jamovi 修改）。

## 引用与许可
引用：Sun, M. *matchfamily* (jamovi module). https://github.com/mingchao-sun/matchfamily 。请同时引用所用方法的原始文献（输出中列出）。

许可：见 `DESCRIPTION`（GPL）；第三方代码与数据见 `NOTICE`。问题与建议：https://github.com/mingchao-sun/matchfamily/issues
