"""Reference: checked score matrices -> sample summaries -> comparisons -> figures/report.

Copy and adapt this example to the study design. The comparison unit is the
biological sample named by `sample`; scores are aggregated by sample medians.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def bh_adjust(values):
    values = np.asarray(values, dtype=float)
    result = np.full(values.shape, np.nan)
    valid = np.flatnonzero(np.isfinite(values))
    order = valid[np.argsort(values[valid])]
    if len(order):
        q = values[order] * len(order) / np.arange(1, len(order) + 1)
        result[order] = np.minimum(1, np.minimum.accumulate(q[::-1])[::-1])
    return result


def sample_medians(scores, meta, keys, min_cells):
    """Return feature-by-unit values and metadata, excluding small units explicitly."""
    rows, columns, excluded = [], [], []
    groups = meta.groupby(keys, observed=True, sort=True, dropna=False) if keys else [((), meta)]
    for key, group in groups:
        key = key if isinstance(key, tuple) else (key,)
        record = dict(zip(keys, key)) | {"n_cells": len(group)}
        if len(group) < min_cells:
            excluded.append(record)
            continue
        unit = f"u{len(rows) + 1}"
        rows.append({"unit": unit, **record})
        columns.append(scores.loc[:, group.index].median(axis=1).rename(unit))
    if not rows:
        raise ValueError("没有分组达到 min_cells，请检查细胞选择与阈值")
    return pd.concat(columns, axis=1), pd.DataFrame(rows).set_index("unit"), pd.DataFrame(excluded, columns=keys + ["n_cells"])


def compare_samples(values, units, spec):
    """One two-group contrast, with BH correction across all tested strata/features."""
    sample, condition = spec["sample"], spec["condition"]
    case, control, design = spec["case"], spec["control"], spec["design"]
    strata = spec["strata"]
    groups = units.groupby(strata, observed=True, sort=True) if strata else [((), units)]
    rows = []
    for key, group in groups:
        key = key if isinstance(key, tuple) else (key,)
        group = group[group[condition].isin([case, control])]
        case_units = group[group[condition] == case]
        control_units = group[group[condition] == control]
        case_ids, control_ids = set(case_units[sample]), set(control_units[sample])
        paired_ids = sorted(case_ids & control_ids)
        if design == "unpaired" and case_ids & control_ids:
            raise ValueError("unpaired 设计中出现跨组相同样本；请核对独立单位或选择 paired")
        for feature in values.index:
            a = values.loc[feature, case_units.index].to_numpy(dtype=float)
            b = values.loc[feature, control_units.index].to_numpy(dtype=float)
            n_pairs, missing_pairs = 0, 0
            if design == "paired":
                n_pairs, missing_pairs = len(paired_ids), len(case_ids ^ control_ids)
                case_lookup = dict(zip(case_units[sample], case_units.index))
                control_lookup = dict(zip(control_units[sample], control_units.index))
                a = values.loc[feature, [case_lookup[s] for s in paired_ids]].to_numpy(dtype=float)
                b = values.loc[feature, [control_lookup[s] for s in paired_ids]].to_numpy(dtype=float)
                delta = a - b
                effect = float(np.median(delta)) if len(delta) else np.nan
                enough = n_pairs >= spec["min_samples"]
                p = (1.0 if np.all(delta == 0) else float(stats.wilcoxon(delta, alternative="two-sided", method="auto").pvalue)) if enough else np.nan
                test = "paired_wilcoxon"
            else:
                effect = float(np.median(a) - np.median(b)) if len(a) and len(b) else np.nan
                enough = min(len(a), len(b)) >= spec["min_samples"]
                p = float(stats.mannwhitneyu(a, b, alternative="two-sided", method="auto").pvalue) if enough else np.nan
                test = "mann_whitney_u"
            rows.append(dict(zip(strata, key)) | {
                "feature": feature, "case": case, "control": control,
                "n_case": len(a), "n_control": len(b), "n_pairs": n_pairs,
                "n_unpaired_excluded": missing_pairs,
                "case_median": float(np.median(a)) if len(a) else np.nan,
                "control_median": float(np.median(b)) if len(b) else np.nan,
                "effect_case_minus_control": effect, "test": test,
                "p_value": p, "status": "tested" if enough else "descriptive_low_n",
            })
    result = pd.DataFrame(rows)
    result["fdr_bh"] = bh_adjust(result["p_value"])
    return result


def heatmap(values, units, output, top_n):
    # Rank features by variation across aggregate units; no p-value filtering.
    selected = values.var(axis=1, ddof=0).nlargest(top_n).index
    plot = values.loc[selected]
    std = plot.std(axis=1, ddof=0).replace(0, 1)
    scaled = plot.sub(plot.mean(axis=1), axis=0).div(std, axis=0)
    labels = [" | ".join(str(v) for k, v in row.items() if k != "n_cells") or unit
              for unit, row in units.iterrows()]
    fig, ax = plt.subplots(figsize=(min(22, max(7, len(labels) * .6)), max(3.5, len(selected) * .32)))
    im = ax.imshow(scaled, cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
    ax.set_xticks(range(len(labels)), labels, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(len(selected)), selected, fontsize=8)
    ax.set_title("Metabolic scores: aggregate medians (row z-score)")
    fig.colorbar(im, ax=ax, label="Row z-score")
    fig.tight_layout()
    fig.savefig(output, dpi=160, bbox_inches="tight")
    plt.close(fig)


def load_analysis(path, config, meta):
    raw = json.loads(Path(path).read_text(encoding="utf-8")) if path else {}
    defaults = {"sample": config.get("sample_column"), "condition": None,
                "case": None, "control": None, "design": "descriptive",
                "strata": [], "min_cells": 1, "min_samples": 3, "top_n": 12}
    unknown = set(raw) - set(defaults)
    if unknown:
        raise ValueError(f"未知结果分析字段: {sorted(unknown)}")
    spec = defaults | raw
    if spec["design"] not in {"descriptive", "paired", "unpaired"}:
        raise ValueError("design 应为 descriptive、paired 或 unpaired")
    if not isinstance(spec["strata"], list) or len(spec["strata"]) != len(set(spec["strata"])):
        raise ValueError("strata 应为不重复的元信息列名列表")
    if spec["sample"] in spec["strata"] or spec["condition"] in spec["strata"]:
        raise ValueError("sample、condition 和 strata 使用不同列")
    for name in ("min_cells", "min_samples", "top_n"):
        if not isinstance(spec[name], int) or spec[name] < 1:
            raise ValueError(f"{name} 应为正整数")
    if spec["design"] != "descriptive":
        if not all(spec[k] for k in ("sample", "condition", "case", "control")):
            raise ValueError("组间比较需要 sample、condition、case 和 control")
        if spec["sample"] == spec["condition"] or spec["case"] == spec["control"]:
            raise ValueError("样本列、条件列及两组取值需要分别指定")
    keys = list(dict.fromkeys(spec["strata"] + [x for x in (spec["sample"], spec["condition"]) if x]))
    if not path:
        keys = list(dict.fromkeys(keys + config.get("group_by", [])))
    for key in keys:
        if key not in meta or meta[key].isna().any():
            raise ValueError(f"元信息列缺失或含空值: {key}")
        meta[key] = meta[key].astype(str)
    if spec["design"] != "descriptive":
        if not {spec["case"], spec["control"]} <= set(meta[spec["condition"]]):
            raise ValueError("case/control 与元信息条件取值不匹配")
    return spec, keys


def build_report(run_dir, output, analysis=None):
    run_dir, output = Path(run_dir).resolve(), Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("结果分析目录已有文件，请使用新目录")
    manifest = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    config = manifest["config"]
    meta_path = run_dir / "cell_metadata.csv"
    meta = pd.read_csv(meta_path, index_col=0, dtype=str) if meta_path.exists() else None
    if meta is not None and meta.index.has_duplicates:
        raise ValueError("元信息细胞 ID 重复")
    output.mkdir(parents=True, exist_ok=True)
    report = ["# 单细胞代谢结果", "", f"输入：`{manifest['input']}`", "",
              f"计算规模：{manifest['n_cells']} 细胞 × {manifest['n_genes']} 基因。", "",
              f"表达层：{config.get('layer')}；通路库：{config.get('database')}；R 种子：{config.get('r_seed')}；scFEA epochs：{config.get('scfea_epochs')}。", "",
              f"[运行清单]({(run_dir / 'run_manifest.json').as_posix()}) · [基因覆盖]({(run_dir / 'coverage.csv').as_posix()}) · [日志]({(run_dir / 'run.log').as_posix()})", "",
              "评分矩阵保持各方法原始量纲；热图按行标准化，仅用于展示相对变化。", "",
              "## 方法与输出", "", "| 方法 | 状态 | 特征数 |", "| --- | --- | --- |"]
    notes, results, evidence = [], [], []
    effective_spec = None
    for method, state in manifest["methods"].items():
        if state["status"] != "complete":
            report.append(f"| {method} | failed: {state.get('error', '')} | — |")
            continue
        for stored_path in state["outputs"]:
            source = run_dir / Path(stored_path).name
            raw = pd.read_csv(source, index_col=0)
            scores = raw.T if method == "scfea" else raw
            scores.columns = scores.columns.astype(str)
            if scores.empty or not np.isfinite(scores.to_numpy(dtype=float)).all():
                raise ValueError(f"输出为空或包含非有限值: {source}")
            if scores.index.has_duplicates or scores.columns.has_duplicates:
                raise ValueError(f"输出 ID 重复: {source}")
            cells = meta.copy() if meta is not None else pd.DataFrame(index=scores.columns)
            if list(cells.index) != list(scores.columns):
                raise ValueError(f"细胞 ID/顺序不匹配: {source}")
            spec, keys = load_analysis(analysis, config, cells)
            effective_spec = spec | {"aggregation_columns": keys}
            values, units, excluded = sample_medians(scores, cells, keys, spec["min_cells"])
            stem = source.stem.removeprefix("scores_")
            table = values.T.rename_axis(index="unit", columns="feature").reset_index().melt(
                id_vars="unit", var_name="feature", value_name="median_score")
            table = table.merge(units, left_on="unit", right_index=True, validate="many_to_one")
            table_name = f"{'sample' if spec['sample'] else 'group'}_scores_{stem}.csv"
            table.to_csv(output / table_name, index=False)
            excluded.to_csv(output / f"excluded_groups_{stem}.csv", index=False)
            comparison = None
            if spec["design"] != "descriptive":
                comparison = compare_samples(values, units, spec)
                comparison.to_csv(output / f"comparison_{stem}.csv", index=False)
            figure = f"heatmap_{stem}.png"
            heatmap(values, units, output / figure, spec["top_n"])
            report.append(f"| {stem} | complete | {len(scores)} |")
            section = [f"### {stem}", "", f"聚合单位 {len(units)} 个；因细胞数阈值排除 {len(excluded)} 个分组。", "",
                       f"[聚合表]({table_name}) · [排除记录](excluded_groups_{stem}.csv)", "", f"![{stem}]({figure})", ""]
            if comparison is not None:
                tested = comparison[comparison.status == "tested"]
                significant = tested[tested.fdr_bh < .05]
                section += [f"[{spec['design']} 比较结果](comparison_{stem}.csv)：检验 {len(tested)} 项，BH FDR < 0.05 为 {len(significant)} 项。", "",
                            "BH 校正在本结果矩阵的全部已检验特征与分层中完成；图中条目按聚合值变异度选取。", ""]
                if len(significant):
                    top = significant.sort_values("fdr_bh").iloc[0]
                    section += [f"最小 FDR 条目：{top['feature']}；效应（case − control）={top['effect_case_minus_control']:.4g}，FDR={top['fdr_bh']:.4g}。", ""]
                elif len(tested):
                    section += ["本次比较未检出 FDR < 0.05 的条目；结合效应方向和样本数解读。", ""]
                largest = comparison.loc[comparison.effect_case_minus_control.abs().nlargest(3).index]
                section += ["效应绝对值较大的条目（描述性排序）：", "", "| 条目 | case − control | FDR |", "| --- | --- | --- |"]
                for _, row in largest.iterrows():
                    context = " / ".join(str(row[k]) for k in spec["strata"])
                    label = (context + ": " if context else "") + str(row["feature"])
                    q = f"{row.fdr_bh:.4g}" if np.isfinite(row.fdr_bh) else "未检验"
                    section.append(f"| {label} | {row.effect_case_minus_control:.4g} | {q} |")
                section.append("")
            results.extend(section)
            evidence.append({"matrix": str(source), "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                             "features": len(scores), "cells": len(scores.columns), "aggregate_units": len(units)})
    if effective_spec is not None:
        (output / "analysis_config.json").write_text(json.dumps(effective_spec, ensure_ascii=False, indent=2), encoding="utf-8")
        notes = ["## 分析设置", "", f"设计：{effective_spec['design']}；样本列：{effective_spec['sample']}；条件列：{effective_spec['condition']}。", "",
                 f"每个样本分组至少 {effective_spec['min_cells']} 个细胞；每组或配对数少于 {effective_spec['min_samples']} 时输出描述值，p/FDR 留空。", "",
                 "配对效应为供体内差值的中位数；非配对效应为两组供体中位数之差。样本聚合量是单细胞评分的中位数。", "",
                 "## 图表与结果", ""]
    report += [""] + notes + results + ["## 解释", "", "GSVA、AUCell、VISION 表示 RNA 驱动的通路评分；scFEA 表示模型推断的模块通量与 balance。", "",
        "结合基因覆盖、样本设计和 scFEA 训练日志解释结果；涉及代谢物浓度或绝对通量的结论需结合相应实验测量。", ""]
    (output / "report.md").write_text("\n".join(report), encoding="utf-8")
    (output / "output_checks.json").write_text(json.dumps(evidence, ensure_ascii=False, indent=2), encoding="utf-8")
    return output / "report.md"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--analysis")
    args = parser.parse_args()
    try:
        print(build_report(args.run_dir, args.output, args.analysis))
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
