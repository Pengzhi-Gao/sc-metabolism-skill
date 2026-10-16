# 参数与口语映射

计算配置供 `scripts/run.py --config` 使用，示例见 [human-all-methods.json](../examples/human-all-methods.json)。JSON 使用 UTF-8；数据、资源和输出路径相对于配置文件所在目录解析。未知字段会报错。

调用时先按 [SKILL.md 的调用规则](../SKILL.md#调用规则) 询问数据、算法或具体需求、关键参数。下表默认值用于提出建议；关键设置按用户回答填写，或在用户明确同意采用建议后使用。已给出的答案直接沿用。

## 计算参数

| 参数 | 默认值 | 含义 |
| --- | --- | --- |
| input / output / species | 必填 | 输入、新结果目录、human 或 mouse |
| methods | ["scmetabolism"] | gsva、scmetabolism、vision、scfea，可选一项或多项 |
| database | KEGG | R 通路方法可选 KEGG/REACTOME；scFEA 使用自己的模块资源 |
| layer | counts | H5AD 的 X/层名；RDS 的精确 layer/assay 名 |
| gene_column | null | H5AD var 中的符号列；null 使用 var_names |
| orientation | genes_by_cells | 表格也支持 cells_by_genes |
| metadata | null | 外部 CSV，首列为与表达矩阵一致的 cell ID；合并 H5AD 元信息时使用不同列名 |
| select | null | {"column":"cell_type","values":["CD8 T"]}；多列组合先生成子集 |
| group_by | [] | 基础运行器按这些列取细胞评分中位数 |
| sample_column | null | 生物学样本列；基础汇总统计样本数，后续示例用于样本聚合 |
| gene_sets_gmt | null | GSVA/VISION 的自定义 GMT；小鼠通路评分需提供 |
| min_genes_per_set | 10 | R 通路输出的最少匹配基因数；GSVA/VISION 在建集时筛选，scMetabolism 按原包 GMT 的精确匹配数过滤返回通路 |
| min_gene_overlap | 0.15 | 运行前参考基因覆盖率下限 |
| pathway_contains | null | 仅筛选基础汇总表名称；原始矩阵保留全库 |
| cores | 1 | 传给 AUCell/VISION；GSVA 和 scFEA 未映射此线程参数 |
| r_seed | 16 | R 评分种子；scFEA 使用原脚本内部种子 |
| scfea_epochs | 100 | scFEA 训练轮次；按训练损失与结果稳定性调整 |
| max_dense_elements | 50000000 | 输入适配的稠密矩阵元素上限；大数据先选子集 |
| rscript / r_libs_user | Rscript / null | R 可执行文件和包库 |
| python | 当前解释器 | scFEA 子进程解释器 |

H5AD/10x 可在稠密转换前选择细胞；CSV/RDS 建议先生成所需子集。VISION 至少 30 个细胞，scFEA 单次最多 10,000 个细胞。CSV/TSV 中列 ID、基因与细胞标识保持唯一。

## 口语要求如何落到代码

| 用户要求 | 配置或后续步骤 |
| --- | --- |
| “只看 CD8 T” | select 对应 cell_type 及其真实取值 |
| “四种都跑” | methods 填四项，依次运行 |
| “按病人看治疗前后” | group_by=[donor_id,condition]；sample_column=donor_id；比较示例使用 paired |
| “两组不同病人” | 比较示例使用 unpaired，并核对独立样本 |
| “各细胞类型分别比较” | 比较配置 strata=[cell_type] |
| “看看分布，没有样本信息” | 描述性聚合和绘图 |
| “重点看糖酵解” | 完整计算后筛选通路名称；scFEA 另按模块资源关联名称 |
| “保留原始零值” | 沿用默认关闭 ALRA/MAGIC 的调用 |

两组比较的参数保存在单独的 comparison.json，字段见 [examples.md](examples.md#e-复用已有矩阵做统计和图表)。它们供结果分析示例使用，与算法计算参数分开记录。

## 输出

基础运行器保存 preflight.json、coverage.csv、cell_metadata.csv（有元信息时）、run_manifest.json、run.log，以及 scores_<method>.csv 或 scfea_flux.csv/scfea_balance.csv。R 矩阵为 pathway × cell，scFEA 为 cell × module/compound。

`summary_*.csv` 为细胞分组中位数。结果分析示例进一步生成 sample_scores/group_scores、comparison、excluded_groups、heatmap、analysis_config、output_checks 和 report.md。比较表的效应为 case − control；配对时先计算供体内差值。

```powershell
python scripts/run.py check --config config.json
python scripts/run.py run --config config.json
```

读取 check 输出的 ready，确认其为 true 后运行。已有 run_manifest.json 的目录会在写入前被拒绝；继续分析既有矩阵时直接调用 summarize_results.py，并指定新的分析输出目录。
