#!/usr/bin/env python3
"""Validated input adapter and dispatcher for the four upstream methods."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


SKILL = Path(__file__).resolve().parents[1]
ROOT = SKILL
UPSTREAM = ROOT / "upstream"
METHODS = {"gsva", "scmetabolism", "vision", "scfea"}
DEFAULTS = {
    "methods": ["scmetabolism"], "database": "KEGG", "layer": "counts",
    "orientation": "genes_by_cells", "min_genes_per_set": 10,
    "min_gene_overlap": 0.15, "cores": 1, "scfea_epochs": 100,
    "max_dense_elements": 50_000_000, "rscript": "Rscript", "python": None,
    "r_libs_user": None,
    "r_seed": 16,
    "metadata": None, "select": None, "group_by": [], "sample_column": None,
    "gene_sets_gmt": None, "gene_column": None, "pathway_contains": None,
}


def absolute_path(value, base):
    p = Path(value).expanduser()
    return (p if p.is_absolute() else base / p).resolve()


def config_load(path):
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("配置必须是 JSON 对象")
    unknown = set(raw) - (set(DEFAULTS) | {"input", "output", "species"})
    if unknown:
        raise ValueError("未知配置字段: " + ", ".join(sorted(unknown)))
    c = {**DEFAULTS, **raw}
    for key in ("input", "output", "species"):
        if not c.get(key):
            raise ValueError(f"缺少必填字段: {key}")
    c["species"] = {"人": "human", "human": "human", "homo_sapiens": "human",
                    "鼠": "mouse", "小鼠": "mouse", "mouse": "mouse", "mus_musculus": "mouse"}.get(str(c["species"]).lower(), str(c["species"]).lower())
    if c["species"] not in {"human", "mouse"}:
        raise ValueError("当前只支持 human 或 mouse；其他物种需要单独适配")
    c["methods"] = [str(x).lower() for x in c["methods"]]
    if not c["methods"] or set(c["methods"]) - METHODS or len(set(c["methods"])) != len(c["methods"]):
        raise ValueError("methods 必须是不重复的 gsva/scmetabolism/vision/scfea 列表")
    c["database"] = str(c["database"]).upper()
    if c["database"] not in {"KEGG", "REACTOME"}:
        raise ValueError("database 必须是 KEGG 或 REACTOME")
    if c["orientation"] not in {"genes_by_cells", "cells_by_genes"}:
        raise ValueError("orientation 必须是 genes_by_cells 或 cells_by_genes")
    if isinstance(c["group_by"], str):
        c["group_by"] = [c["group_by"]]
    if not isinstance(c["group_by"], list):
        raise ValueError("group_by 必须是列名数组")
    for key in ("input", "output", "metadata", "gene_sets_gmt", "r_libs_user"):
        if c.get(key):
            c[key] = str(absolute_path(c[key], path.parent))
    if c["species"] != "human" and "scmetabolism" in c["methods"]:
        raise ValueError("scMetabolism 原生代谢基因集是人类符号；小鼠请选择 scFEA 或提供匹配 GMT 后使用 GSVA/VISION")
    if "scmetabolism" in c["methods"] and c["gene_sets_gmt"]:
        raise ValueError("scMetabolism 使用原包内置基因集，不能接收 gene_sets_gmt；自定义 GMT 请选择 GSVA/VISION")
    if c["species"] != "human" and any(x in c["methods"] for x in ("gsva", "vision")) and not c["gene_sets_gmt"]:
        raise ValueError("非人类 GSVA/VISION 需要 gene_sets_gmt")
    for key in ("cores", "scfea_epochs", "min_genes_per_set", "max_dense_elements", "r_seed"):
        c[key] = int(c[key])
        if c[key] < 1:
            raise ValueError(f"{key} 必须大于零")
    c["min_gene_overlap"] = float(c["min_gene_overlap"])
    if not 0 <= c["min_gene_overlap"] <= 1:
        raise ValueError("min_gene_overlap 应在 0 和 1 之间")
    return c


def _read_table(path, orientation):
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    # Check original headers before pandas can rename duplicates with '.1'.
    with path.open(encoding="utf-8-sig", newline="") as handle:
        header = next(csv.reader(handle, delimiter=sep), [])[1:]
    if not header or any(not value.strip() for value in header):
        raise ValueError("矩阵列 ID 为空")
    if len(header) != len(set(header)):
        raise ValueError("矩阵列 ID 重复；请保留唯一的基因/细胞标识")
    frame = pd.read_csv(path, sep=sep, index_col=0)
    if orientation == "cells_by_genes":
        frame = frame.T
    return frame


def _read_rds(path, c, tmp):
    rscript = shutil.which(c["rscript"]) or c["rscript"]
    out = tmp / "from_rds.csv"
    meta = tmp / "from_rds_metadata.csv"
    cmd = [rscript, str(SKILL / "scripts" / "export_rds.R"), str(path), str(out), str(meta),
           c["layer"], str(c["max_dense_elements"])]
    proc = subprocess.run(cmd, text=True, capture_output=True, check=False, env=_r_env(c))
    if proc.returncode:
        raise RuntimeError(f"RDS 转换失败: {proc.stderr[-2000:]}")
    return _read_table(out, "genes_by_cells"), pd.read_csv(meta, index_col=0)


def _read_10x(path):
    try:
        from scipy.io import mmread
    except ImportError as exc:
        raise RuntimeError("10x 输入需要 scipy") from exc
    def find(name):
        for candidate in (name, name + ".gz"):
            p = path / candidate
            if p.exists():
                return p
        return None
    matrix = find("matrix.mtx")
    genes = find("features.tsv") or find("genes.tsv")
    barcodes = find("barcodes.tsv")
    if not all((matrix, genes, barcodes)):
        raise ValueError("10x 目录缺少 matrix.mtx、features.tsv/genes.tsv 或 barcodes.tsv")
    opener = lambda p: gzip.open(p, "rt", encoding="utf-8") if p.suffix == ".gz" else open(p, encoding="utf-8")
    with opener(genes) as f:
        features = list(csv.reader(f, delimiter="\t"))
    with opener(barcodes) as f:
        cells = [row[0] for row in csv.reader(f, delimiter="\t")]
    with opener(matrix) as f:
        m = mmread(f).tocsr()
    if m.shape != (len(features), len(cells)):
        raise ValueError("10x 矩阵与 features/barcodes 数量不匹配")
    keep = [i for i, row in enumerate(features) if len(row) < 3 or row[2] == "Gene Expression"]
    names = [features[i][1] if len(features[i]) > 1 else features[i][0] for i in keep]
    return m[keep, :], names, cells


def load_input(c, tmp):
    path = Path(c["input"])
    if not path.exists():
        raise FileNotFoundError(path)
    meta = None
    external_meta = None
    if c["metadata"]:
        external_meta = pd.read_csv(c["metadata"], index_col=0)
        external_meta.index = external_meta.index.astype(str)
        if external_meta.index.duplicated().any():
            raise ValueError("元信息中的细胞 ID 重复")
    if path.is_dir():
        data, genes, cells = _read_10x(path)
        if c["select"] and external_meta is not None and c["select"]["column"] in external_meta:
            selector = c["select"]
            values = selector["values"]
            if isinstance(values, str):
                values = [values]
            matched = external_meta.index[external_meta[selector["column"]].astype(str).isin([str(v) for v in values])]
            keep = [i for i, cell in enumerate(cells) if cell in matched]
            if not keep:
                raise ValueError("select 未选中任何细胞")
            data = data[:, keep]
            cells = [cells[i] for i in keep]
        if data.shape[0] * data.shape[1] > c["max_dense_elements"]:
            raise ValueError("10x 矩阵转密集格式超过 max_dense_elements；请先按细胞类型/样本拆分")
        frame = pd.DataFrame(data.toarray(), index=genes, columns=cells)
    elif path.suffix.lower() in {".csv", ".tsv", ".txt"}:
        frame = _read_table(path, c["orientation"])
    elif path.suffix.lower() == ".h5ad":
        try:
            import anndata as ad
        except ImportError as exc:
            raise RuntimeError("H5AD 输入需要 anndata 和 h5py") from exc
        a = ad.read_h5ad(path)
        if c["select"]:
            values = c["select"]["values"]
            if isinstance(values, str):
                values = [values]
            column = c["select"]["column"]
            if column in a.obs:
                selected = a.obs[column].astype(str).isin([str(v) for v in values])
                if not selected.any():
                    raise ValueError("select 未选中任何细胞")
                a = a[selected].copy()
            elif external_meta is not None and column in external_meta:
                selected_cells = set(external_meta.index[external_meta[column].astype(str).isin([str(v) for v in values])])
                selected = [cell in selected_cells for cell in a.obs_names.astype(str)]
                if not any(selected):
                    raise ValueError("select 未选中任何细胞")
                a = a[selected].copy()
        x = a.X if c["layer"] in {"X", "x", None} else a.layers[c["layer"]]
        if x.shape[0] * x.shape[1] > c["max_dense_elements"]:
            raise ValueError("H5AD 矩阵转密集格式超过 max_dense_elements；请先创建目标细胞子集")
        genes = a.var[c["gene_column"]].astype(str) if c["gene_column"] else a.var_names.astype(str)
        frame = pd.DataFrame(x.T.toarray() if hasattr(x, "toarray") else np.asarray(x.T),
                             index=genes, columns=a.obs_names.astype(str))
        meta = a.obs.copy()
        meta.index = meta.index.astype(str)
    elif path.suffix.lower() == ".rds":
        frame, meta = _read_rds(path, c, tmp)
    else:
        raise ValueError("输入应为 h5ad、rds、基因×细胞 CSV/TSV 或 10x 目录")
    if frame.empty or frame.shape[0] * frame.shape[1] > c["max_dense_elements"]:
        raise ValueError("矩阵为空或超过 max_dense_elements")
    if frame.index.isna().any() or frame.columns.isna().any():
        raise ValueError("基因或细胞 ID 存在缺失")
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    if frame.columns.duplicated().any():
        raise ValueError("细胞 ID 重复，无法安全匹配元信息")
    if frame.index.duplicated().any():
        raise ValueError("基因 ID 重复；请先明确合并或去重规则")
    if set(c["methods"]) & {"gsva", "vision"} and frame.index.str.upper().duplicated().any():
        raise ValueError("GSVA/VISION 输入含大小写归一后重复的基因 ID；请先明确映射规则")
    frame = frame.apply(pd.to_numeric, errors="raise")
    values = frame.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("表达矩阵必须是非负有限数值")
    if external_meta is not None:
        meta = external_meta if meta is None else meta.join(external_meta, how="left", validate="one_to_one")
    if meta is not None:
        missing = frame.columns.difference(meta.index)
        if len(missing):
            raise ValueError(f"元信息缺少 {len(missing)} 个细胞 ID")
        meta = meta.loc[frame.columns].copy()
    if c["select"]:
        if meta is None:
            raise ValueError("select 需要细胞元信息")
        selector = c["select"]
        col = selector["column"]
        if col not in meta:
            raise ValueError(f"select 列不存在: {col}")
        vals = selector["values"]
        if isinstance(vals, str):
            vals = [vals]
        chosen = meta.index[meta[col].astype(str).isin([str(v) for v in vals])]
        if not len(chosen):
            raise ValueError("select 未选中任何细胞")
        frame = frame.loc[:, chosen]
        meta = meta.loc[chosen]
    for col in c["group_by"] + ([c["sample_column"]] if c["sample_column"] else []):
        if meta is None or col not in meta:
            raise ValueError(f"元信息缺少列: {col}")
        if meta[col].isna().any():
            raise ValueError(f"元信息列存在缺失: {col}")
    return frame, meta


def gmt_path(c):
    if c["gene_sets_gmt"]:
        return Path(c["gene_sets_gmt"])
    filename = "KEGG_metabolism_nc.gmt" if c["database"] == "KEGG" else "REACTOME_metabolism.gmt"
    return UPSTREAM / "scMetabolism" / "data" / filename


def read_gmt(path):
    if not path.is_file():
        raise FileNotFoundError(f"代谢基因集不存在: {path}")
    sets = {}
    with path.open(encoding="utf-8") as f:
        for line in f:
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) >= 3:
                sets[fields[0]] = [x for x in fields[2:] if x]
    if not sets:
        raise ValueError("GMT 没有有效基因集")
    return sets


def coverage_table(c, frame):
    sets = read_gmt(gmt_path(c)) if set(c["methods"]) & {"gsva", "vision", "scmetabolism"} else {}
    reference = {g.upper() for genes in sets.values() for g in genes}
    observed = {g.upper() for g in frame.index}
    observed_exact = set(frame.index)
    rows = [{"feature_type": "pathway", "feature": name, "reference_genes": len(set(genes)),
             "matched_genes": len({g.upper() for g in genes} & observed)} for name, genes in sets.items()]
    if "scfea" in c["methods"]:
        fname = "module_gene_m168.csv" if c["species"] == "human" else "module_gene_complete_mouse_m168.csv"
        module = pd.read_csv(UPSTREAM / "scFEA-master" / "data" / fname, index_col=0)
        scfea_reference = {str(g) for g in module.to_numpy().ravel() if pd.notna(g)}
        for name, values in module.iterrows():
            genes = {str(g) for g in values if pd.notna(g)}
            rows.append({"feature_type": "scfea_module", "feature": name,
                         "reference_genes": len(genes), "matched_genes": len(genes & observed_exact)})
    else:
        scfea_reference = set()
    result = {
        "pathway_genes": len(reference),
        "pathway_genes_matched": len(reference & observed),
        "pathway_overlap_fraction": round(len(reference & observed) / len(reference), 4) if reference else None,
        "scmetabolism_exact_overlap_fraction": (
            round(len({g for genes in sets.values() for g in genes} & observed_exact) /
                  len({g for genes in sets.values() for g in genes}), 4)
            if reference and "scmetabolism" in c["methods"] else None),
        "scfea_genes": len(scfea_reference),
        "scfea_genes_matched": len(scfea_reference & observed_exact),
        "scfea_overlap_fraction": round(len(scfea_reference & observed_exact) / len(scfea_reference), 4) if scfea_reference else None,
    }
    return pd.DataFrame(rows), result


def dependency_report(c):
    report = {"python": sys.executable, "numpy": bool(importlib.util.find_spec("numpy")),
              "pandas": bool(importlib.util.find_spec("pandas"))}
    path = Path(c["input"])
    if path.is_dir():
        report["scipy"] = bool(importlib.util.find_spec("scipy"))
    if "gsva" in c["methods"]:
        report["scipy"] = bool(importlib.util.find_spec("scipy"))
    if path.suffix.lower() == ".h5ad":
        report["anndata"] = bool(importlib.util.find_spec("anndata"))
    if path.suffix.lower() == ".rds" or set(c["methods"]) & {"gsva", "vision", "scmetabolism"}:
        rscript = shutil.which(c["rscript"]) or (c["rscript"] if Path(c["rscript"]).is_file() else None)
        report["Rscript"] = bool(rscript)
        if rscript:
            packages = []
            if "gsva" in c["methods"]:
                packages.append("GSVA")
            if "vision" in c["methods"]:
                packages.append("VISION")
            if "scmetabolism" in c["methods"]:
                packages.extend(["scMetabolism", "AUCell"])
            if packages:
                names = ",".join(f"'{x}'" for x in sorted(set(packages)))
                code = f"cat(paste(vapply(c({names}), requireNamespace, logical(1), quietly=TRUE), collapse=','))"
                proc = subprocess.run([rscript, "-e", code], capture_output=True, text=True,
                                      check=False, env=_r_env(c))
                vals = proc.stdout.strip().split(",") if proc.returncode == 0 else []
                report.update({f"R:{p}": vals[i] == "TRUE" if i < len(vals) else False
                               for i, p in enumerate(sorted(set(packages)))})
                for package in sorted(set(packages)):
                    if report[f"R:{package}"]:
                        version = subprocess.run([rscript, "-e", f"cat(as.character(utils::packageVersion('{package}')))"],
                                                 capture_output=True, text=True, check=False, env=_r_env(c))
                        if version.returncode == 0:
                            report[f"R:{package}_version"] = version.stdout.strip()
    if "scfea" in c["methods"]:
        py = c["python"] or sys.executable
        proc = subprocess.run([py, "-c", "import torch, matplotlib, tqdm, pandas, numpy"],
                              capture_output=True, text=True, check=False)
        report["scfea_python"] = py
        report["scfea_dependencies"] = proc.returncode == 0
        if proc.returncode == 0:
            v = subprocess.run([py, "-c", "import torch; print(torch.__version__)"], capture_output=True, text=True, check=False)
            if v.returncode == 0:
                report["torch_version"] = v.stdout.strip()
        if proc.returncode:
            report["scfea_dependency_error"] = proc.stderr[-500:]
    return report


def source_version(name):
    folder = UPSTREAM / name
    if (folder / ".git").is_dir():
        proc = subprocess.run(["git", "-C", str(folder), "rev-parse", "HEAD"], capture_output=True, text=True, check=False)
        if proc.returncode == 0:
            return proc.stdout.strip()
    archive = UPSTREAM / {"VISION-master": "VISION-master.zip", "scFEA-master": "scFEA.zip", "scop-main": "scop.zip"}.get(name, "")
    if archive.is_file():
        return "zip-sha256:" + hashlib.sha256(archive.read_bytes()).hexdigest()
    snapshots = UPSTREAM / "SOURCES.json"
    if snapshots.is_file():
        record = json.loads(snapshots.read_text(encoding="utf-8")).get(name, {})
        if record.get("tree_sha256"):
            return "snapshot-tree-sha256:" + record["tree_sha256"]
    return "unavailable"


def _r_env(c):
    env = dict(os.environ)
    if c.get("r_libs_user"):
        env["R_LIBS_USER"] = c["r_libs_user"]
    if os.name == "nt":
        for key in ("LC_ALL", "LC_CTYPE"):
            if env.get(key) == "C.UTF-8":
                env.pop(key)
    return env


def preflight(c, frame, meta):
    coverage, overlap = coverage_table(c, frame)
    deps = dependency_report(c)
    warnings = []
    too_many_scfea = "scfea" in c["methods"] and frame.shape[1] > 10000
    if too_many_scfea:
        warnings.append("scFEA 原仓库提示超过 10000 细胞耗时显著；运行阶段会拒绝该规模")
    too_few_vision = "vision" in c["methods"] and frame.shape[1] < 30
    if too_few_vision:
        warnings.append("VISION 至少需要 30 个细胞用于邻域分析")
    if not c["sample_column"] and c["group_by"]:
        warnings.append("没有独立样本/供体列，只输出描述性组汇总，不做细胞级显著性检验")
    if c["species"] == "human" and frame.index.to_series().str.match(r"^ENS(G|MUSG)\d+").mean() > 0.5:
        warnings.append("多数基因是 Ensembl ID；请用 gene_column 或预先映射为物种对应 gene symbol")
    low = [k for k in ("pathway_overlap_fraction", "scmetabolism_exact_overlap_fraction", "scfea_overlap_fraction")
           if overlap[k] is not None and overlap[k] < c["min_gene_overlap"]]
    if low:
        warnings.append("基因集覆盖率低于 min_gene_overlap: " + ", ".join(low))
    pathway_coverage = coverage.loc[coverage["feature_type"] == "pathway"] if len(coverage) else coverage
    insufficient_sets = bool(len(pathway_coverage)) and not (pathway_coverage["matched_genes"] >= c["min_genes_per_set"]).any()
    if insufficient_sets:
        warnings.append("没有任何代谢基因集达到 min_genes_per_set")
    ready = (not low and not too_many_scfea and not too_few_vision and not insufficient_sets
             and not any(v is False for v in deps.values())
             and not any("Ensembl ID" in w for w in warnings))
    out = {"ready": ready, "input": c["input"], "species": c["species"], "methods": c["methods"],
           "n_genes": int(frame.shape[0]), "n_cells": int(frame.shape[1]), "metadata_columns": list(meta.columns) if meta is not None else [],
           "overlap": overlap, "dependencies": deps, "warnings": warnings,
           "gene_sets_gmt": str(gmt_path(c)) if not set(c["methods"]) <= {"scfea"} else None}
    return out, coverage


def _run_logged(cmd, cwd, logfile, env=None):
    with logfile.open("a", encoding="utf-8") as f:
        f.write("\nCOMMAND: " + json.dumps(cmd, ensure_ascii=False) + "\n")
        f.flush()
        proc = subprocess.run(cmd, cwd=cwd, stdout=f, stderr=subprocess.STDOUT, text=True, env=env, check=False)
    if proc.returncode:
        raise RuntimeError(f"后端返回 {proc.returncode}；查看 {logfile}")


def write_sparse_gsva_input(frame, tmp):
    """Preserve zero entries as sparse data for GSVA's single-cell algorithm."""
    try:
        from scipy import sparse
        from scipy.io import mmwrite
    except ImportError as exc:
        raise RuntimeError("GSVA 稀疏输入需要 scipy") from exc
    target = tmp / "expression_gsva.mtx"
    matrix = sparse.csc_matrix(frame.to_numpy(dtype=float, copy=False))
    mmwrite(target, matrix)
    pd.Series(frame.index, dtype=str).to_csv(str(target) + ".genes.csv", index=False, header=False)
    pd.Series(frame.columns, dtype=str).to_csv(str(target) + ".cells.csv", index=False, header=False)
    return target


def summarize(scores, meta, c, outpath, feature_name="pathway"):
    if meta is None or not c["group_by"]:
        return
    if c["pathway_contains"]:
        scores = scores.loc[[x for x in scores.index if c["pathway_contains"].lower() in str(x).lower()]]
    rows = []
    groups = meta.groupby(c["group_by"], dropna=False, observed=True, sort=True)
    for key, group in groups:
        key = key if isinstance(key, tuple) else (key,)
        cells = group.index.intersection(scores.columns)
        if not len(cells):
            continue
        med = scores.loc[:, cells].median(axis=1)
        n_samples = group[c["sample_column"]].nunique() if c["sample_column"] else None
        for feature, value in med.items():
            row = dict(zip(c["group_by"], key))
            row.update({feature_name: feature, "median_score": float(value), "n_cells": len(cells), "n_samples": n_samples})
            rows.append(row)
    pd.DataFrame(rows).to_csv(outpath, index=False)


def validate_scores(scores, cells, features=None):
    """Validate the shared feature-by-cell output contract."""
    if scores.empty or scores.index.has_duplicates or scores.columns.has_duplicates:
        raise ValueError("输出矩阵为空或含重复 ID")
    if scores.index.isna().any() or scores.columns.isna().any():
        raise ValueError("输出矩阵含缺失 ID")
    if not np.isfinite(scores.to_numpy(dtype=float)).all():
        raise ValueError("输出矩阵包含非有限数值")
    if list(scores.columns) != list(cells):
        raise ValueError("输出细胞 ID/顺序不匹配")
    if features is not None and list(scores.index) != list(features):
        raise ValueError("输出特征 ID/顺序与模块资源不匹配")


def ensure_new_output(output):
    if (output / "run_manifest.json").exists():
        raise FileExistsError("输出目录已有 run_manifest.json；请选择新目录，保留既往分析")


def run(c, frame, meta, check, coverage, output):
    if not check["ready"]:
        raise RuntimeError("前置检查未通过；查看 preflight.json")
    ensure_new_output(output)
    output.mkdir(parents=True, exist_ok=True)
    coverage.to_csv(output / "coverage.csv", index=False)
    tmp = output / "work"
    tmp.mkdir(exist_ok=True)
    matrix = tmp / "expression.csv"
    if any(method != "gsva" for method in c["methods"]):
        frame.to_csv(matrix)
    if meta is not None:
        meta.to_csv(output / "cell_metadata.csv", index_label="cell_id")
    logfile = output / "run.log"
    manifest = {"started_at_utc": datetime.now(timezone.utc).isoformat(), "input": c["input"],
                "config": {k: v for k, v in c.items() if k not in {"python", "rscript"}},
                "n_genes": int(frame.shape[0]), "n_cells": int(frame.shape[1]), "overlap": check["overlap"],
                "runtime": check["dependencies"],
                "sources": {x: source_version(x) for x in ("GSVA", "scMetabolism", "VISION-master", "scFEA-master", "scop-main")},
                "methods": {}}
    for method in c["methods"]:
        try:
            if method == "scfea":
                if frame.shape[1] > 10000:
                    raise ValueError("scFEA 单次最多 10000 细胞；请基于研究问题选择子集")
                src = UPSTREAM / "scFEA-master"
                data_dir = src / "data"
                mouse = c["species"] == "mouse"
                suffix = "complete_mouse_c70_m168.csv" if mouse else "c70_m168.csv"
                # Upstream prefixes res_dir with './' when saving its loss plot;
                # an absolute Windows path would become invalid ('./D:\\...').
                relative_result = os.path.relpath(output, src).replace("\\", "/")
                cmd = [c["python"] or sys.executable, str(SKILL / "scripts" / "scfea_no_impute.py"),
                       "--data_dir", str(data_dir), "--input_dir", str(tmp), "--res_dir", relative_result,
                       "--test_file", matrix.name,
                       "--moduleGene_file", "module_gene_complete_mouse_m168.csv" if mouse else "module_gene_m168.csv",
                       "--stoichiometry_matrix", "cmMat_" + suffix,
                       "--cName_file", "cName_" + suffix,
                       "--sc_imputation", "False", "--train_epoch", str(c["scfea_epochs"]),
                       "--output_flux_file", str(output / "scfea_flux.csv"),
                       "--output_balance_file", str(output / "scfea_balance.csv")]
                # Upstream scFEA converts inference tensors directly to NumPy.
                # Force CPU even on hosts with CUDA, where that conversion fails.
                env = dict(os.environ, CUDA_VISIBLE_DEVICES="-1")
                _run_logged(cmd, src, logfile, env)
                outputs = [output / "scfea_flux.csv", output / "scfea_balance.csv"]
                module_file = data_dir / ("module_gene_complete_mouse_m168.csv" if mouse else "module_gene_m168.csv")
                module_ids = pd.read_csv(module_file, index_col=0).index.tolist()
                compound_names = pd.read_csv(data_dir / ("cName_" + suffix), nrows=0).columns.tolist()
                for p, expected in zip(outputs, (module_ids, compound_names)):
                    x = _read_table(p, "genes_by_cells")
                    validate_scores(x.T, frame.columns, expected)
                    summarize(x.T, meta, c, output / ("summary_" + p.stem + ".csv"), "feature")
            else:
                target = output / f"scores_{method}.csv"
                rscript = shutil.which(c["rscript"]) or c["rscript"]
                method_input = write_sparse_gsva_input(frame, tmp) if method == "gsva" else matrix
                cmd = [rscript, str(SKILL / "scripts" / "score.R"), str(method_input), str(target),
                       method, str(gmt_path(c)), c["database"], str(c["cores"]),
                       str(c["min_genes_per_set"]), str(c["r_seed"])]
                _run_logged(cmd, ROOT, logfile, _r_env(c))
                scores = _read_table(target, "genes_by_cells")
                validate_scores(scores, frame.columns)
                summarize(scores, meta, c, output / f"summary_{method}.csv")
                outputs = [target]
            manifest["methods"][method] = {
                "status": "complete", "outputs": [str(p) for p in outputs],
                "backend_input_format": "matrix_market_sparse" if method == "gsva" else "csv_dense",
            }
        except Exception as exc:
            manifest["methods"][method] = {"status": "failed", "error": str(exc)}
    manifest["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
    (output / "run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    if any(x["status"] != "complete" for x in manifest["methods"].values()):
        raise RuntimeError("至少一个方法失败；查看 run_manifest.json 和 run.log")
    return manifest


def main():
    parser = argparse.ArgumentParser(description="单细胞 RNA 代谢推断：检查输入并调用原始算法")
    parser.add_argument("action", choices=["check", "run"])
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    path = Path(args.config).resolve()
    c = config_load(path)
    output = Path(c["output"])
    ensure_new_output(output)
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="scmet_input_", dir=output) as d:
        frame, meta = load_input(c, Path(d))
        check, coverage = preflight(c, frame, meta)
        (output / "preflight.json").write_text(json.dumps(check, ensure_ascii=False, indent=2), encoding="utf-8")
        print(json.dumps(check, ensure_ascii=False, indent=2))
        if args.action == "run":
            run(c, frame, meta, check, coverage, output)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
