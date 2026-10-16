# 示例的选择与修改

这些示例提供数据处理、算法调用和结果交付的写法。调用时将所需配置复制到分析项目目录，再按实际输入修改；`scripts/` 保留通用适配，项目代码记录这次研究的选择。

直接使用示例时按下方命令从 Skill 目录调用。若把运行示例复制到项目，将 `run_workflow.py` 与 `summarize_results.py` 放在同一目录，并加 `--skill-dir 'D:\codex\omic-skill\single-cell-metabolism'` 指定通用适配代码的位置。结果分析示例可单独复制运行。

## A. 人类 H5AD：从口语到配置

用户说：“只看 CD8 T，四种方法都算，比较同一病人治疗前后，重点看糖酵解。”

这句话已明确细胞范围、算法和配对需求。继续询问数据路径与物种，以及表达层、通路库、匹配阈值、scFEA 训练轮次和样本统计阈值；给出建议值供用户选择。文件中的实际列名和取值可先检查后展示：

```python
from pathlib import Path
import anndata as ad

source = Path('D:/project/data/cells.h5ad')
a = ad.read_h5ad(source, backed='r')
print(a.shape, list(a.layers.keys()))
print(a.obs.columns.tolist(), a.var.columns.tolist())
print(a.obs.groupby(['donor_id', 'condition', 'cell_type'], observed=True).size())
print(a.var_names[:5])
a.file.close()
```

取得关键参数的回答或用户采用建议的授权后，确认 counts 层、gene symbol 列、分组值和供体关系，再修改 [human-all-methods.json](../examples/human-all-methods.json)：

- `input` / `output`：真实文件和新的结果目录。
- `layer`：counts，或已经确认适合该方法的 X；记录其归一化状态。
- `gene_column`：符号列名；var_names 已是符号时填 null。
- `select`：`{"column":"cell_type","values":["CD8 T"]}`。
- `group_by`：`["donor_id","condition"]`，保留供体与条件分层。
- `sample_column`：donor_id。

修改 [paired-comparison.json](../examples/paired-comparison.json) 中的真实条件值，例如 case=After、control=Before。`paired` 使用同一独立样本的两种条件；每位供体同一条件的技术重复先按研究设计合并。

```powershell
$py = 'D:\WorkBuddy\env\single_cell_metabolism_py312\Scripts\python.exe'
$skill = 'D:\codex\omic-skill\single-cell-metabolism'
& $py "$skill\examples\run_workflow.py" --config 'D:\project\config.json' --analysis 'D:\project\comparison.json'
```

输出在 config 的 output 目录，结果分析在其 `analysis/` 下。重点通路在全库计算和多重校正之后筛选，完整矩阵与比较表保留。

## B. 较大 H5AD 或复杂选择：先保存研究子集

```python
from pathlib import Path
import anndata as ad

source = Path('D:/project/data/cells.h5ad')
target = Path('D:/project/data/cd8_selected.h5ad')
if target.exists():
    raise FileExistsError(target)
a = ad.read_h5ad(source, backed='r')
keep = a.obs['cell_type'].eq('CD8 T') & a.obs['condition'].isin(['Before', 'After'])
if not keep.any():
    raise ValueError('筛选结果为空，请核对元信息取值')
subset = a[keep].to_memory()
a.file.close()
target.parent.mkdir(parents=True, exist_ok=True)
subset.write_h5ad(target)
subset.obs.groupby(['donor_id', 'condition'], observed=True).size().to_csv(target.with_suffix('.cell_counts.csv'))
```

根据记录的数量判断样本是否仍具备比较条件。基因是 Ensembl ID 时，使用该物种、该注释版本的映射表填 gene symbol；多个条目映射到同一符号时，按数据层明确合并规则。原始 counts 可求和，归一化值需按原预处理方案处理；保存映射表和未匹配数量。

## C. 其他输入格式

### Seurat RDS

在 R 中检查 assay、layer 和元信息；如需子集，在导出前选择：

```r
library(SeuratObject)
obj <- readRDS("D:/project/data/cells.rds")
print(names(obj@assays))
print(Layers(obj[["RNA"]]))
print(colnames(obj@meta.data))
selected <- rownames(obj@meta.data)[obj@meta.data$cell_type == "CD8 T"]
if (!length(selected)) stop("no selected cells")
obj <- subset(obj, cells = selected)
saveRDS(obj, "D:/project/data/cd8_selected.rds")
```

配置 input 指向新 RDS、layer 指向一个精确层名。分层的 counts.* 按实验设计先合并，或选定具体层。转换器读取 Seurat 的 RNA assay（没有 RNA 时用 DefaultAssay），也接受 SCE/SummarizedExperiment 的指定 assay。基因标识取矩阵行名。

### CSV/TSV

```python
import pandas as pd

# 原始表头须保留唯一 cell ID；生成矩阵时明确 genes × cells。
counts = pd.DataFrame({'cell_1@A': [5, 0], 'cell_2@B': [0, 8]}, index=['HK1', 'GPI'])
counts.to_csv('D:/project/data/counts.csv')
pd.DataFrame({'donor_id': ['A', 'B'], 'condition': ['Control', 'Treatment']},
             index=counts.columns).to_csv('D:/project/data/metadata.csv')
```

这里的两基因矩阵仅演示文件格式；实际计算使用具有足够代谢基因覆盖的表达矩阵。配置 input 为 counts.csv，metadata 为 metadata.csv，orientation 为 genes_by_cells。细胞 × 基因表使用 cells_by_genes。

### 10x MatrixMarket

input 指向包含 matrix.mtx(.gz)、features.tsv/genes.tsv(.gz)、barcodes.tsv(.gz) 的目录，metadata 指向细胞元信息 CSV。Gene Expression 特征会被读取；select 可依据外部元信息先选细胞。

## D. 只用一种方法或更换物种

在四算法配置基础上修改：

| 任务 | 配置变化 |
| --- | --- |
| 默认人类通路评分 | methods=["scmetabolism"] |
| 人类 REACTOME | database="REACTOME"，methods 选需要的 R 方法 |
| 自定义 GMT 的 GSVA/VISION | methods=["gsva","vision"]，gene_sets_gmt 指向文件 |
| 小鼠模块通量 | species="mouse"，methods=["scfea"] |
| 小鼠通路评分 | species="mouse"，methods=["gsva","vision"]，gene_sets_gmt 使用小鼠符号 |

所有方法的 API 示例在 [score.R](../scripts/score.R) 与 [run.py](../scripts/run.py)。GSVA 使用 `gsvaParam`；scMetabolism 使用 `method="AUCell"`；VISION 使用 `Vision → analyze → getSignatureScores`；scFEA 使用仓库 CLI 加兼容启动器。

## E. 复用已有矩阵做统计和图表

```powershell
$py = 'D:\WorkBuddy\env\single_cell_metabolism_py312\Scripts\python.exe'
$example = 'D:\codex\omic-skill\single-cell-metabolism\examples\summarize_results.py'
& $py $example --run-dir 'D:\project\results\metabolism' --output 'D:\project\results\analysis_v2' --analysis 'D:\project\comparison.json'
```

`comparison.json` 的字段：

| 字段 | 示例与含义 |
| --- | --- |
| sample | donor_id：独立生物学样本列 |
| condition | condition：实验组别列 |
| case / control | Treatment / Control：效应为 case − control |
| design | paired、unpaired 或 descriptive |
| strata | [] 或 ["cell_type"]；每个层内分别比较 |
| min_cells | 20：每个样本 × 条件 × 分层的最少细胞数；按数据调整 |
| min_samples | 3：开始检验的最少配对数或每组样本数；这是执行阈值，研究功效由样本设计决定 |
| top_n | 12：热图展示的高变异特征数 |

配对分析仅使用两种条件均存在的供体，并在比较表记录未配对数量。非配对分析要求两组供体独立。样本数不足时保留描述性效应，p/FDR 留空。BH 校正覆盖每个输出矩阵内全部已检验特征与分层；跨方法确认性分析需要另外定义检验家族。

样本聚合是单细胞评分的中位数。需要原始 counts pseudobulk、协变量调整、混合模型或纵向设计时，复用输入、聚合与输出代码并替换统计步骤，记录模型公式与比较方向。

## F. 结果交付

示例自动生成 sample_scores/group_scores、comparison、excluded_groups、热图、analysis_config、output_checks 和 report.md。Codex 阅读这些结果后，将报告补成以下内容：

1. 本次研究问题、细胞与样本数、选择条件。
2. 使用的方法、基因覆盖、归一化和模型训练设置。
3. 主要通路/模块的变化方向、效应、供体一致性和对应图表。
4. 独立单位、检验设计、FDR 与排除记录；没有比较设计时报告描述结果。
5. 各方法完成状态、结果路径以及回答当前问题所需的下一步实验或分析。

热图每行进行了标准化；比较表中的效应保留原评分量纲。阅读报告时先核对对应方法和样本范围，再解释生物学意义。
