---
name: single-cell-metabolism
description: 根据单细胞 RNA 数据和口语化要求，完成代谢通路评分、模块通量推断、样本层比较及结果交付。当用户要求分析单细胞代谢、比较细胞群或条件的代谢状态，或指定 GSVA、scMetabolism、VISION、scFEA 处理单细胞表达矩阵时使用；提供输入适配、完整流程和可修改的代码示例。
---

# 单细胞代谢推断

## 1. 规则：何时触发、如何调用

### 触发规则

- 用户提供单细胞 RNA 数据，要求计算代谢活性、代谢通路、模块通量或代谢物 balance。
- 用户要求比较不同细胞类型、条件或供体的代谢状态。
- 用户指定 GSVA、scMetabolism、VISION、scFEA，或明确调用本 Skill。

输入为 H5AD、Seurat/SingleCellExperiment RDS、10x MatrixMarket 目录或 CSV/TSV。质谱峰表对应实验代谢组工作流；本 Skill 的输入是 RNA 表达矩阵。

### 调用规则

每次开始一项分析，先明确询问以下三类信息。用户在当前请求或此前对话中已经给出的答案直接沿用，只补问尚未明确的部分。

1. **使用什么数据？**“请提供要分析的文件或目录路径。数据是什么物种、什么格式？有细胞类型、实验分组和样本/供体信息吗？”用户不清楚文件结构时，先读取结构，再向用户列出相关列和表达层供选择。
2. **使用哪个算法，或者具体有什么需求？**“你想使用 GSVA、scMetabolism/AUCell、VISION、scFEA，还是多种方法一起运行？如果还没选算法，请说明想看通路活性、推断通量，还是比较某些细胞群或实验条件。”根据需求给出推荐及简短理由，由用户选择；用户明确交由 Codex 选择时，说明采用的方法并继续。
3. **关键参数怎样设置？**根据选定的数据、算法和比较需求，提前询问下表中适用的项目。给出数据中的实际候选值和建议值，允许用户一次回复“按建议设置”。

| 适用情况 | 运行前的具体问题 | 对应设置 |
| --- | --- | --- |
| 所有分析 | “分析全部细胞，还是指定细胞类型、样本或子集？使用哪一层表达数据，其是否已归一化？” | select、layer、输入预处理 |
| 通路评分 | “使用 KEGG、REACTOME，还是你自己的 GMT？每条通路至少匹配多少基因、整体覆盖率下限是否采用建议值 10 和 15%？” | database、gene_sets_gmt、min_genes_per_set、min_gene_overlap |
| scFEA | “本次用于流程试跑还是正式分析？训练轮次如何设置？示例从 100 epochs 起步，正式分析结合损失与结果稳定性调整；是否需要限制细胞范围？” | scfea_epochs、细胞选择与训练检查 |
| 组间比较 | “比较哪两组，哪组作为对照？哪一列标识独立样本/供体？两组是否配对，是否按细胞类型分别比较？” | case、control、sample、condition、design、strata |
| 样本层统计 | “每个样本分组至少保留多少细胞、至少多少独立样本或配对才做检验？是否采用示例的 20 个细胞、3 个样本/配对？” | min_cells、min_samples |

提问采用简短、可直接回答的表述；算法确定后再展开该算法的参数问题。关键答案尚未明确时，先做输入与环境的只读检查。得到回答或用户明确授权采用建议后，将数据、算法和关键参数汇总写入配置并开始计算。示例阈值是可调整的起点，样本量是否足以回答问题按研究设计判断。

随后按以下约定执行：

- 阅读对应示例，复用最接近任务的配置和代码，按确认的列名、数据层和实验设计修改；将本次配置和分析代码保存在项目目录。
- 计算细节沿用原包 API 和已有适配代码；新需求在示例基础上扩展。随机种子、日志等常规执行参数使用已记录默认值；用户提出资源或复现要求时再调整。
- 本机环境使用 `D:\WorkBuddy\env`，路径中的 `_` 和 `@` 原样保留。首次运行或环境缺失时读 [environment.md](references/environment.md)。
- 组间推断以独立生物学样本为单位；多因素、纵向重复测量按具体研究设计调整示例。

| 用户意图 | 方法 | 输出 |
| --- | --- | --- |
| 人类代谢通路活性，希望推荐算法 | 推荐 scMetabolism / AUCell，按上述提问规则确定 | 通路 × 细胞评分 |
| 指定 GSVA | GSVA 单细胞稀疏算法 | 基因集 × 细胞评分 |
| 指定 VISION | VISION | 签名 × 细胞评分 |
| 通量、代谢物平衡或 scFEA | scFEA | 细胞 × 模块通量及 balance |
| 四种方法都算 | 依次运行四项 | 各方法矩阵、汇总和图表 |

人类通路评分可建议 KEGG，并询问是否采用 REACTOME 或自定义库。小鼠使用 scFEA，或提供匹配小鼠的 GMT 后用 GSVA/VISION。其他物种需要先适配基因集与输入检查。GSVA/VISION 支持自定义 GMT；scMetabolism 使用安装包内置人类基因集。

## 2. 完整流程：从数据到输出

**询问数据、算法/需求和关键参数 → 输入整理 → 配置与检查 → 算法计算 → 样本汇总与比较 → 图表 → 结果报告。**

| 步骤 | 执行动作 | 留下的结果 |
| --- | --- | --- |
| ① 提问并明确任务 | 按调用规则询问数据、算法或需求、关键参数；结合文件结构明确细胞范围、比较组和独立单位 | 用户回答、选定参数与数据概览 |
| ② 整理输入 | 选目标细胞；使用与方法相容的非负表达层；整理 gene symbol；核对 cell ID 与元信息；需要映射/合并基因时记录规则和映射表 | 分析输入、细胞数量与转换记录 |
| ③ 整理配置 | 复制配置示例，填真实路径和列名；选择方法、数据库和分组；比较设计单独保存 | config.json；需要比较时加 comparison.json |
| ④ 检查并计算 | 执行 check，读取 ready 和覆盖率；修正实际问题后 run；依次调用原包实现 | preflight、coverage、评分/通量矩阵、运行日志与清单 |
| ⑤ 汇总与比较 | 校验矩阵的 ID、维度和有限值；按样本 × 条件 × 细胞类型取评分中位数；按设计比较并校正多重检验 | 样本汇总表；比较表；排除记录 |
| ⑥ 绘图 | 生成方法各自的代谢热图；按问题增加目标通路的供体点图、配对连线或已有 UMAP 的评分着色 | 带单位、分组和样本数的图 |
| ⑦ 交付 | 根据图表撰写主要发现，说明输入、筛选、方法、覆盖率、统计单位、效果方向和完成状态，链接到结果 | report.md、图、表、配置与复用后的代码 |

输入整理方式和命令见 [examples.md](references/examples.md)，参数含义见 [config.md](references/config.md)。默认关闭 ALRA/MAGIC 填补；大数据先选择与问题对应的细胞子集。scFEA 示例设置为 100 epochs，按运行前确认的轮次和训练检查方案执行。

基础运行器的 `summary_*.csv` 是细胞分组中位数，`sample_column` 用于统计样本数。需要组间推断时，接着运行第 3 部分的样本分析示例：按独立样本聚合，配对设计用 Wilcoxon，非配对设计用 Mann–Whitney U，BH 校正后报告效应与 FDR。缺少独立样本信息时输出描述性表和图。

GSVA、AUCell、VISION 表示 RNA 驱动的通路评分；scFEA 表示模型推断通量与 balance。各方法按各自量纲解释，代谢物浓度和绝对通量的结论结合实验测量。报告中注明稀疏 GSVA、预处理和 scFEA 训练设置。

交付完成的判断：所选方法的完成/失败状态明确；成功矩阵已检查；对应的汇总、图和结果报告已生成；用户能够从报告定位配置、代码与输出。部分方法失败时，报告已完成结果并说明具体失败步骤。

## 3. 示例代码：调用时的参考起点

先阅读对应示例，再修改路径、列名、细胞选择与统计设计；保留原算法调用和输出含义。示例参数可按数据和问题调整。

| 场景 | 参考文件 | 如何复用 |
| --- | --- | --- |
| 人类数据，选择细胞并跑四算法 | [human-all-methods.json](examples/human-all-methods.json) | 修改 input、output、layer、gene_column、select 和元信息列；删减 methods 可运行单方法 |
| 从检查到报告的一次完整调用 | [run_workflow.py](examples/run_workflow.py) | 复用 check → run → 汇总/报告顺序；配置相对于自身所在目录解析 |
| 配对样本比较 | [paired-comparison.json](examples/paired-comparison.json) | 填 sample、condition、case、control；strata 可填 cell_type |
| 描述分析、非配对比较、热图和报告 | [summarize_results.py](examples/summarize_results.py) | 复用样本中位数、统计与 BH、绘图函数；对复杂设计替换比较函数 |
| H5AD、RDS、CSV/10x 的输入适配 | [examples.md](references/examples.md) | 参考读取、子集与配置写法；按输入实际结构调整 |
| 原始算法 API 与上游关系 | [score.R](scripts/score.R)、[run.py](scripts/run.py)、[methods.md](references/methods.md) | 直接参考已使用的包调用、兼容处理与矩阵方向 |
| scFEA 上游示例数据 | [scfea-source-demo.json](examples/scfea-source-demo.json) | 用于算法调用练习；更换输出目录后运行 |

本机完整调用：

```powershell
$py = 'D:\WorkBuddy\env\single_cell_metabolism_py312\Scripts\python.exe'
$skill = 'D:\codex\omic-skill\single-cell-metabolism'
& $py "$skill\examples\run_workflow.py" --config 'D:\project\config.json' --analysis 'D:\project\comparison.json'
```

省略 `--analysis` 时按配置中的样本列和分组生成描述性汇总、热图和报告。有已完成矩阵时，直接调用 `summarize_results.py`，复用原始计算结果。
