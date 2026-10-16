# 来源与实现关系

推送 [原文](https://mp.weixin.qq.com/s/XfDnfgwrB57ocuKhItLVXw) 截图显示 `GSVA`、`RunMetabolism`、`VISION`、`scFEA` 四条 R 脚本。正文点名 scFEA、scMetabolism；`RunMetabolism` 对应 `scMetabolism`，`VISION` 和 `GSVA` 也可作为 scMetabolism 的评分后端。

| 方法 | 原仓库 | 本 Skill 调用 | 产出 |
| --- | --- | --- | --- |
| GSVA | https://github.com/rcastelo/GSVA | `GSVA::gsva(gsvaParam(...))` | 代谢基因集评分 |
| scMetabolism | https://github.com/wu-yc/scMetabolism | `scMetabolism::sc.metabolism(..., method="AUCell")` | KEGG/REACTOME 代谢通路评分 |
| VISION | https://github.com/YosefLab/VISION | `VISION::Vision`, `analyze`, `getSignatureScores` | 代谢签名评分 |
| scFEA | https://github.com/changwn/scFEA | `scfea_no_impute.py` 启动未改动的 `src/scFEA.py` | 模块通量与 metabolite balance |
| scop 参考封装 | https://github.com/mengxu98/scop | 阅读 `RunGSVA.R`、`RunMetabolism.R`、`RunscFEA.R`，运行时不依赖 scop | 参数与输出行为参考 |

本 Skill 读取 `upstream/scMetabolism/data/*.gmt` 作为 GSVA/VISION 的代谢基因集；原包 scMetabolism 用其安装版自带基因集。GSVA 使用 MatrixMarket 稀疏矩阵和 `kcdf="auto"`，调用原包面向单细胞数据的稀疏计算；首次真实数据测试中，稠密 CSV/强制 Poisson 在 25,712 基因×400 细胞上耗时约 16 分钟。

默认关闭 ALRA/MAGIC 填补。scFEA 的兼容启动器用空的 `magic` 模块满足原脚本的导入要求，并把已删除的 `DataFrame.append` 映射到等价的 `pandas.concat`，随后执行仓库中的原始计算脚本。scFEA 的 GPU 输出路径直接调用 `.numpy()`，因此包装器固定使用 CPU。单次 scFEA 运行最多 10,000 个细胞；处理更大数据时，按研究问题分层或抽样，并报告筛选范围。

R 方法每次运行调用 `set.seed(r_seed)`。AUCell 对表达值相同的基因做随机排序，固定种子使同一输入和软件版本下的评分可复现。配置和运行清单记录所用种子。

`min_genes_per_set` 在 GSVA/VISION 中筛选输入基因集，在 scMetabolism 中按安装包 GMT 与输入基因的精确交集过滤返回通路。各方法最终保留的通路数以结果为准。包装器校验输出的细胞 ID、唯一特征和有限值；scFEA 还与当前模块、代谢物资源逐项核对。

后续结果分析参考 `examples/summarize_results.py`：按样本、条件及指定分层取细胞评分中位数，再依据配对或非配对设计比较。该步骤使用评分层的样本汇总。

仓库源码保存在 Skill 内的 `upstream`。GSVA、scMetabolism 为浅克隆；VISION、scFEA、scop 为 GitHub 源码压缩包。Skill 调用这些原始实现，软件许可见各仓库文件。

scFEA 仓库的 LICENSE 允许学术非商业用途，营利用途要求单独获得授权。使用者应在实际部署场景中核对许可。
