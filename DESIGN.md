# matchfamily 设计文档：匹配与加权家族

状态：**已实现，v1.0.0（2026-11-16）**。以下为开工前的设计，与最终实现的差别见本文件末尾“实现后的差别”。模块名暂定 `matchfamily`，菜单组 `Match`，界面与帮助文档用英文。
沿用 ivfamily、rdfamily 的做法：一个家族一个模块，每个分析配 Method notes、Data and setup guide、示例数据、NEWS、NOTICE，验证从设计阶段就带着对照。

## 1. 目标

在"选择可观测"（所有混杂变量都被观测到）的假设下，让用户在一个模块里完成观察数据因果推断的标准流程：
估计倾向得分并检查重叠，用多种估计量（回归调整、加权、双重稳健、匹配、熵平衡）估计处理效应，检查协变量平衡，
最后评估结论对未观测混杂的敏感程度。

这是 Angrist-Pischke 之外使用最广的识别思路，Stata 的 `teffects` 是内置命令，所以数值对照可以做得很严。

## 2. 分析清单

### 第一批（v1.0，7 个）

| 名称 | 菜单标题 | 内容 | 主要输出 |
|---|---|---|---|
| `psmodel` | Propensity Score and Overlap | logit 或 probit 倾向得分；重叠检查和修剪 | 系数表、c 统计量（AUC）、两组得分分布图、共同支撑范围、按规则（min-max、Crump 0.1 到 0.9、用户阈值）修剪后的样本量 |
| `balance` | Covariate Balance | 处理组与对照组的协变量平衡，处理前后对比 | 每个协变量的标准化均值差、方差比、KS 统计量；加权或匹配前后；Love plot；超过阈值（0.1）的个数 |
| `wtest` | Weighting and Doubly Robust Estimators | 回归调整（RA）、IPW、AIPW、IPWRA；ATE、ATT、ATC | 效应估计、标准误（M-估计，计入倾向得分估计误差）、置信区间、权重诊断（极端权重、有效样本量）、可选修剪 |
| `nnmatch` | Matching Estimators | 最近邻匹配（k 个近邻、有放回、卡尺、对部分变量精确匹配）、马氏距离或倾向得分距离、偏差校正（Abadie-Imbens）；ATE、ATT | 效应估计、Abadie-Imbens 标准误、匹配后样本的平衡表、每个处理单位的匹配数分布 |
| `ebal` | Entropy Balancing | 熵平衡（Hainmueller 2012）：精确平衡均值（可选方差、偏度）的权重 | 权重摘要、有效样本量、平衡表、加权后的效应估计 |
| `matchcmp` | Compare Matching and Weighting Estimators | 同一估计目标下并排比较朴素差值、RA、IPW、AIPW、IPWRA、最近邻匹配、熵平衡 | 估计量表、森林图、各方法的平衡摘要 |
| `sensobs` | Sensitivity to Unobserved Confounding | E-value；匹配样本的 Rosenbaum 界 | E-value 点估计和置信限；Gamma 与 p 值（符号检验、Wilcoxon 符号秩）；解释说明 |

### 第二批（v1.x）
- 倾向得分分层（子分类）、重叠加权（ATO）、CBPS、协变量精确匹配（CEM）、最优匹配与全匹配、无放回匹配。
- Cinelli-Hazlett 稳健性值与 Oster δ（回归调整的遗漏变量敏感性）。
- 连续处理（广义倾向得分）、多值处理、聚类与面板结构。
- 双重机器学习（DML）、因果森林：计算量大，要先解决 jamovi 里的异步计算和运行时间。

## 3. 选项设计

jamovi 里一个变量只能放一个框（ivover 测试时踩过），所以设计成"额外框"而不是"子集框"。

| 区域 | 选项 |
|---|---|
| 变量 | Outcome（`balance`、`psmodel` 不需要）、Treatment（二值 0/1，两水平因子按第二水平为 1）、Covariates（倾向得分模型的协变量，数值或名义）、Additional outcome-model covariates（只进入结果模型的变量，`wtest` 和 `nnmatch` 的偏差校正用）、Exact-match variables（`nnmatch`，精确匹配的变量，额外于 Covariates）、Weights（抽样权重）、Cluster |
| 倾向得分 | 链接（logit 或 probit）、交互项与平方项（默认无，用户在数据里生成）、修剪规则与阈值 |
| 估计目标 | ATE、ATT、ATC |
| 估计 | `wtest`：估计量（RA、IPW、AIPW、IPWRA）、权重归一化（Hajek 或 Horvitz-Thompson）；`nnmatch`：近邻数 k、距离、卡尺、偏差校正、是否纳入所有并列近邻（默认与 teffects 一致）；`ebal`：平衡的矩阶数 |
| 推断 | 标准误类型（解析、聚类）、可选自助法（重复次数、种子）、置信水平 |
| Output（末尾折叠） | Method notes、Data and setup guide |

要求：每个选项组默认折叠；错误提示指出具体原因（处理变量不是 0/1、某一组没有观测、完全分离、重叠不足、极端权重）；
缺失值删除并计数；每个分析在重叠严重不足时给出警告而不是静默估计。

## 4. 实现方案

**引擎自写，只依赖基础 R（`stats::glm`、矩阵运算），不包装 MatchIt、WeightIt。** 理由：
1. 版本不一致吃过亏：jamovi 打包的包版本与本机不同（rdhte、rdrobust），自写引擎行为稳定。
2. 标准误的定义是这类方法最容易出问题的地方（是否计入倾向得分估计误差、Abadie-Imbens 方差、并列近邻的处理），
   必须与 Stata `teffects` 对齐，包装的包给不出这些中间量。
3. 方法本身不复杂：IPW、AIPW、RA 是 M-估计；匹配是近邻搜索加 Abadie-Imbens 方差；熵平衡是对偶问题的牛顿迭代；Rosenbaum 界和 E-value 是闭式。

结构：
```
matchfamily/
├── DESCRIPTION  LICENSE  NEWS.md  NOTICE.md  README.md
├── jamovi/  0000.yaml  00refs.yaml  <name>.{a,r,u}.yaml
├── R/  engine_ps.R       # 倾向得分（glm）、重叠、修剪
│       engine_weight.R   # RA、IPW、AIPW、IPWRA 及影响函数标准误
│       engine_match.R    # 最近邻搜索、Abadie-Imbens 方差、偏差校正
│       engine_ebal.R     # 熵平衡
│       engine_balance.R  # 标准化均值差、方差比、KS
│       engine_sens.R     # E-value、Rosenbaum 界
│       common.R  guide.R  method_notes.R  <name>.b.R
├── data/  示例数据
├── tests/  smoke.R  edge.R  reference.R  ui_checklist.md
└── dev/  gen_*.py（批量生成 yaml）、make_data.R、icon
```

性能：最近邻搜索对处理组和对照组做分块计算，倾向得分距离（一维）用排序加二分把复杂度降到 O(N log N)；马氏距离按块做，样本很大时给出提示。

## 5. 示例数据

| 数据 | 来源 | 用途 |
|---|---|---|
| LaLonde 职业培训 | `MatchIt::lalonde`（GPL-2，Dehejia-Wahba 的样本，处理组与 PSID 对照） | 所有分析；已知实验基准约 1,800 美元（NSW 随机实验） |
| 模拟：强混杂，真实 ATE 与 ATT 已知 | `dev/make_data.R` | 估计量恢复真值；重叠充分 |
| 模拟：重叠不足 | `dev/make_data.R` | 修剪、极端权重、偏差校正的演示 |
| 模拟：隐藏混杂 | `dev/make_data.R` | `sensobs`：E-value 与 Rosenbaum 界的演示 |

Cattaneo 的 `cattaneo2`（Stata 手册数据）不能随模块再分发，只在对照脚本里使用。

## 6. 验证计划（设计阶段就确定）

| 分析 | R 参考 | Stata 参考 | 备注 |
|---|---|---|---|
| `psmodel`：系数、c 统计量 | `stats::glm`、`pROC` 不需要，自行用秩和公式验证 | `logit`、`lroc` | 简单 |
| `wtest`：RA、IPW、AIPW、IPWRA 的 ATE 和 ATT 及标准误 | `WeightIt` 加 `sandwich`（需安装） | **`teffects ra / ipw / aipw / ipwra`（内置，M-估计标准误）** | 最重要的一组，数值应到 1e-6 |
| `nnmatch`：ATE、ATT、Abadie-Imbens 标准误、偏差校正 | `Matching::Match`（需安装）、`MatchIt` | **`teffects nnmatch`、`teffects psmatch`、`psmatch2`（已装）** | 并列近邻和方差估计的细节要逐项对齐 |
| `ebal`：权重、平衡 | `WeightIt` 或 `ebal` 包 | `ebalance`（SSC，需安装） | 权重应到 1e-6 |
| `balance`：SMD、方差比、KS | `cobalt`（需安装） | `pstest`（已装）、`tebalance summarize`（teffects 后） | 注意 SMD 分母的定义（合并方差还是处理组方差） |
| `sensobs`：E-value | `EValue` 包（需安装） | `evalue`（SSC，需安装） | 闭式 |
| `sensobs`：Rosenbaum 界 | `rbounds` 包（需安装） | `rbounds`（SSC，需安装） | 符号检验的界是二项分布，可手算 |

验证规则（沿用本项目的教训）：
1. 所有测试在 **jamovi 自带 R 与随模块打包的库**下再跑一遍；
2. 每个分析在界面里按清单点一遍，尤其是变量框的行为；
3. 对照脚本固定数据和随机数种子，容差写进脚本，结果写入 `PLAN.md`；
4. Stata 脚本注释里不要写 `/*`，批处理前先 `capture file close _all`；
5. 选项必须真的起作用：每个选项都要有"改变后结果变化"的检验，防止像 rdhte 的核函数那样被悄悄忽略。

## 7. 已知难点与风险

| 风险 | 说明 | 对策 |
|---|---|---|
| IPW 标准误 | 是否计入倾向得分估计误差会差一个量级；teffects 用联立 M-估计 | 一开始就按联立 M-估计写影响函数，对照 teffects |
| 并列近邻 | 离散协变量下并列很多，teffects 包含所有并列近邻，MatchIt 随机取一个 | 默认与 teffects 一致，随机取一个作为选项，文档写明 |
| Abadie-Imbens 方差 | 有放回匹配下用匹配次数 K_M(i) 的方差估计，异方差版本需要近邻估计条件方差 | 对照 `teffects nnmatch, vce(robust, nn(2))`，不一致时以 Stata 为准并写明 |
| 重叠不足 | 极端权重使估计不稳定，用户容易忽略 | 每个分析都显示最大权重、有效样本量，重叠差时给明确警告 |
| 倾向得分完全分离 | `glm` 不收敛或系数发散 | 检测并给出提示：哪个协变量几乎完美预测处理 |
| 平衡指标定义 | SMD 的分母各软件不同 | 文档写明用的是哪种（默认处理组方差，与 cobalt 的 ATT 默认一致），并提供选项 |
| 解释责任 | 用户容易把匹配当成"控制了混杂" | Method notes 明确"只对可观测混杂成立"，`sensobs` 紧跟在结果后 |
| 因果估计目标 | ATE、ATT 混淆 | 每个表的标题写明估计目标 |

## 8. 里程碑

| 阶段 | 内容 | 交付物 |
|---|---|---|
| M0 | 脚手架：DESCRIPTION、0000.yaml、`gen_*.py`、示例数据、图标 | 可安装的空模块 |
| M1 | `engine_ps.R`、`engine_weight.R`；`psmodel`、`wtest` | RA、IPW、AIPW、IPWRA 与 `teffects` 对上 |
| M2 | `engine_balance.R`、`engine_ebal.R`；`balance`、`ebal` | 对上 `cobalt` 或 `pstest`、`ebalance` |
| M3 | `engine_match.R`；`nnmatch` | 对上 `teffects nnmatch`、`psmatch` |
| M4 | `sensobs`、`matchcmp` | E-value、Rosenbaum 界对上 |
| M5 | 帮助文档、界面实测、打包库下重跑测试、图标 | v1.0 候选版 |

## 9. 需要你确认的决策

1. **自写引擎还是包装 MatchIt 和 WeightIt？**（推荐自写，理由见第 4 节。）
2. **v1 估计目标**：ATE、ATT、ATC；重叠加权（ATO）放 v1.x？（推荐如此。）
3. **匹配范围**：v1 只做有放回最近邻（与 `teffects nnmatch` 一致），无放回、最优、全匹配放 v1.x？（推荐如此。）
4. **敏感性分析**：放进本家族的 `sensobs`，v1 只做 E-value 和 Rosenbaum 界，Cinelli-Hazlett 和 Oster δ 放 v1.x？（推荐如此。）
5. **模块名与菜单组**：`matchfamily` 和 `Match`，还是别的？
6. **示例数据**：LaLonde 数据通过安装 CRAN 包 `MatchIt` 提取（同时它也是一个对照包）。可以安装吗？

## 10. 需要在本机安装的包（确认后再装）

- R（CRAN）：`MatchIt`、`WeightIt`、`cobalt`、`Matching`、`EValue`、`rbounds`；已有 `sandwich`。
- Stata（SSC）：`ebalance`、`evalue`、`rbounds`；`teffects` 内置，`psmatch2` 和 `pstest` 已装。

## 实现后的差别（2026-11-16）

- 7 个分析全部实现并装入 jamovi；验证结果见 `README.md` 的 Validation 表和根目录 `PLAN.md`。
- `nnmatch` 没有“倾向得分距离”：点估计与 `teffects psmatch` 一致，但 Abadie-Imbens (2016) 的估计倾向得分方差修正无法复现（Stata 的调整量随 k 变化，手册公式与论文公式均对不上）。引擎里保留实验性的 `mf_psmatch`，未接入界面，列为 v1.x。
- `wtest`、`ebal`、`nnmatch`、`matchcmp` 的 Estimand 默认显示 ATE、ATT、ATC 全部；`balance`、`sensobs` 仍为单选。
- `nnmatch` 在有卡尺或精确匹配变量时剔除无匹配的单位并报告（Stata 报错）。
- 示例数据：LaLonde 用 `.omv`（jamovi 导入 CSV 会把 `age`、`educ` 猜成 Nominal，导致共线）。
- `sensobs` 的 Rosenbaum 界只对 ATT 配对，数据模式下的 E-value 用 AIPW 估计；未做 Stata `evalue`、`rbounds` 对照（对 R 包 `EValue`、`rbounds` 已对）。
- 熵平衡对 Stata `ebalance` 的对照受其默认容差 0.015 影响（已在 `15_ebalance.do` 改为 1e-6，可重跑）。
