"""Reference entrypoint: inspect config -> check -> score -> summarize/report.

Usage: python run_workflow.py --config config.json [--analysis comparison.json]
When copying both example .py files into a project, pass --skill-dir explicitly.
Reuse this sequence; adapt input preparation and interpretation to the study.
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

from summarize_results import build_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--analysis")
    parser.add_argument("--skill-dir", default=str(Path(__file__).resolve().parents[1]))
    args = parser.parse_args()
    config_path = Path(args.config).resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    output = Path(config["output"])
    output = output if output.is_absolute() else config_path.parent / output
    runner = Path(args.skill_dir).resolve() / "scripts" / "run.py"
    if not runner.is_file():
        raise FileNotFoundError("请用 --skill-dir 指向包含 scripts/run.py 的 Skill 目录")
    command = [sys.executable, str(runner)]
    subprocess.run(command + ["check", "--config", str(config_path)], check=True)
    preflight = json.loads((output / "preflight.json").read_text(encoding="utf-8"))
    if not preflight["ready"]:
        raise RuntimeError("前置检查未通过；按 preflight.json 修正后重新运行")
    computation = subprocess.run(command + ["run", "--config", str(config_path)], check=False)
    if (output / "run_manifest.json").exists():
        report = build_report(output, output / "analysis", args.analysis)
        print(f"REPORT: {report}")
    if computation.returncode:
        raise RuntimeError("部分方法运行失败，查看运行清单及报告中的方法状态")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
