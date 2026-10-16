# sc-metabolism-skill

从单细胞 RNA 表达数据出发，使用 **GSVA、scMetabolism/AUCell、VISION 和 scFEA** 分析代谢通路评分、模块通量及代谢物 balance。

Skill 名称为 `single-cell-metabolism`，入口是 [SKILL.md](SKILL.md)。它包含三部分：

1. **规则**：先询问使用的数据、算法或研究需求、关键参数。
2. **完整流程**：输入整理 → 检查与计算 → 样本汇总与比较 → 图表 → 报告。
3. **示例代码**：复用已有配置和调用代码，按数据列名与实验设计修改。

## 能做什么

| 方法 | 主要输出 |
| --- | --- |
| GSVA | 代谢基因集 × 细胞评分，使用原包单细胞稀疏算法 |
| scMetabolism/AUCell | 代谢通路 × 细胞富集评分 |
| VISION | 代谢签名 × 细胞评分 |
| scFEA | 细胞 × 模块通量、代谢物 balance |

输入支持 H5AD、Seurat/SCE RDS、10x MatrixMarket 目录及 CSV/TSV。人类通路评分可选 KEGG/REACTOME；GSVA/VISION 可用自定义 GMT；scFEA 提供人类和小鼠模块。

结果分析示例提供样本中位数汇总、配对/非配对比较、BH 校正、热图和 Markdown 报告。RNA 评分与模型通量按各自量纲解释。

## 安装与使用

将仓库克隆到 Codex 的 Skill 目录，并保留文件夹名 `single-cell-metabolism`：

```powershell
git clone https://github.com/Pengzhi-Gao/sc-metabolism-skill.git "$env:USERPROFILE\.codex\skills\single-cell-metabolism"
```

已有同名 Skill 时，将仓库克隆到新的目录供核对。仓库包含运行所需的 `upstream/` 源码和资源；Python、R 及软件包的安装见 [运行环境](references/environment.md)。示例中的 `D:\WorkBuddy\env` 是开发机环境路径，使用时改为本机路径。

在 Codex 中调用：

> 使用 $single-cell-metabolism 分析我的单细胞数据，先询问数据、算法需求和关键参数。

直接执行完整示例：

```powershell
python examples/run_workflow.py --config /path/to/config.json --analysis /path/to/comparison.json
```

从 [四算法配置](examples/human-all-methods.json) 和 [配对比较配置](examples/paired-comparison.json) 开始修改。省略 `--analysis` 时生成描述性汇总、热图和报告。输入适配、单方法调用及已有结果再分析见 [示例说明](references/examples.md)。

## 目录

| 路径 | 用途 |
| --- | --- |
| [SKILL.md](SKILL.md) | 触发/调用规则、完整流程、示例索引 |
| [examples/](examples/) | 可修改的配置、完整调用、统计、绘图和报告示例 |
| [scripts/](scripts/) | 输入适配、检查与原始算法调用 |
| [references/](references/) | 参数、环境、上游来源与使用说明 |
| [tests/](tests/) | 输入、输出与样本统计测试 |
| [upstream/](upstream/) | 原始算法源码、资源与许可证 |
| [upstream/SOURCES.json](upstream/SOURCES.json) | 随仓库发布的上游源码版本和内容摘要 |

## 验证

2026-10-16 的人类 H5AD 实测使用 400 个细胞、25,712 个基因及 4 位配对患者，四种方法完成计算，并生成样本比较表、5 张热图和报告。GSVA/scMetabolism/VISION 分别输出 68/73/73 条通路；scFEA 输出 168 个模块和 70 个 balance 条目。该数据中的配对比较未检出 FDR < 0.05 的条目。

```powershell
python -m unittest discover -s tests -v
```

完整实测覆盖上述人类 H5AD 场景；其他输入、物种和自定义基因集按各自数据验证。

## 来源与许可

上游仓库和调用关系见 [methods.md](references/methods.md)。第三方源码保留各自的许可证和说明。scFEA 的 LICENSE 规定学术非商业用途及营利用途的单独授权条件，见 [scFEA LICENSE](upstream/scFEA-master/LICENSE)。
