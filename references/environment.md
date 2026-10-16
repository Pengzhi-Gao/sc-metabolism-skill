# 运行环境

检查与 CSV 输入需要 Python 3.10+、`numpy`、`pandas`。10x 需要 `scipy`；H5AD 需要 `anndata`、`h5py`。结果分析示例还使用 `scipy` 和 `matplotlib`。RDS 的 Seurat 转换使用 SeuratObject 5 的 `Layers`/`LayerData` API；SCE/SummarizedExperiment 转换需要 SummarizedExperiment。

R 方法需要 `GSVA`、`scMetabolism`、`VISION`、`AUCell` 及其 R 依赖。在独立 R 包库中安装这些包。包装器使用 `gsvaParam` 参数对象 API，选择支持该 API 的 GSVA 版本。

```r
install.packages(c("pak", "BiocManager"))
BiocManager::install(c("GSVA", "AUCell", "GSEABase"))
pak::pak(c("wu-yc/scMetabolism", "YosefLab/VISION@v3.0.2"))
```

scFEA 使用独立 Python 环境，安装 CPU 版 `torch`、`matplotlib`、`tqdm`、`pandas`、`numpy`。包装器使用 CPU、关闭 MAGIC 填补，并处理原脚本的顶层 `magic` 导入。scFEA 每次对输入重新训练；少量 epoch 适合流程测试，科学分析需检查训练收敛与结果稳定性。

本机使用 `D:\WorkBuddy\env\single_cell_metabolism_py312\Scripts\python.exe` 执行 `run.py`，并将配置的 `python` 指向同一路径。该隔离环境使用 Python 3.12.14 和 PyTorch 2.14.1+cpu，从已通过真实数据测试的工作区环境复制包建立。其他机器可先建虚拟环境，再安装 [requirements-scfea.txt](../requirements-scfea.txt)。H5AD 与 10x 输入分别还需 `anndata`/`h5py` 与 `scipy`，按需要安装。

本机 2026-10-15 的真实数据测试已将 R 4.6.1、Rtools45 和 R 包库放在 `D:\WorkBuddy\env`：`R-4.6.1\bin\Rscript.exe`、`rtools45`、`Rlib-4.6`。已验证 GSVA 2.6.12、AUCell 1.34.0、VISION 3.0.2、scMetabolism 0.2.2。配置中填写 `rscript` 与 `r_libs_user` 后即可独立运行 R 方法；scFEA 使用同目录下的隔离 Python 3.12 环境。GSVA 运行需要 `scipy`，包装器会将表达矩阵转为 MatrixMarket 稀疏格式，使用 GSVA 单细胞稀疏算法。

安装后执行 `check`，核验依赖与输入数据的基因覆盖率。检查通过后运行所选方法，核对输出矩阵维度、细胞 ID 和计算日志。

本机基础 Python 解释器目前位于 Codex 的运行时缓存，虚拟环境中的 pyvenv.cfg 记录了该路径。迁移机器时，在 `D:\WorkBuddy\env` 下安装基础 Python 并重建虚拟环境，随后按所需输入格式安装依赖。

运行器在本 Skill 目录内定位 `upstream`，其中保存 GSVA、scMetabolism、VISION-master、scFEA-master 和 scop-main。移动 Skill 时保留内部目录结构。源码与调用关系见 [methods.md](methods.md)。
