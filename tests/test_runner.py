import importlib.util
import json
import tempfile
import unittest
from unittest import mock
from pathlib import Path

import numpy as np
import pandas as pd


MODULE = Path(__file__).resolve().parents[1] / "scripts" / "run.py"
TEST_TMP = Path(__file__).resolve().parents[1] / ".test_tmp"
TEST_TMP.mkdir(exist_ok=True)
spec = importlib.util.spec_from_file_location("metabolism_runner", MODULE)
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class RunnerContracts(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=TEST_TMP)
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)

    def config(self, **overrides):
        data = {"input": "counts.csv", "output": "out", "species": "human", "methods": ["scfea"], **overrides}
        path = self.base / "config.json"
        path.write_text(json.dumps(data), encoding="utf-8")
        return runner.config_load(path)

    def test_literal_ids_and_metadata_selection(self):
        pd.DataFrame({"cell_1@a": [1, 2], "cell_2@b": [3, 4]}, index=["GENE_A", "GENE@B"]).to_csv(self.base / "counts.csv")
        pd.DataFrame({"cell_type": ["T", "B"], "donor_id": ["p_1", "p_2"]},
                     index=["cell_1@a", "cell_2@b"]).to_csv(self.base / "meta.csv")
        c = self.config(metadata="meta.csv", select={"column": "cell_type", "values": ["T"]},
                        group_by=["cell_type"], sample_column="donor_id")
        with tempfile.TemporaryDirectory(dir=TEST_TMP) as d:
            matrix, metadata = runner.load_input(c, Path(d))
        self.assertEqual(list(matrix.index), ["GENE_A", "GENE@B"])
        self.assertEqual(list(matrix.columns), ["cell_1@a"])
        self.assertEqual(metadata.loc["cell_1@a", "donor_id"], "p_1")

    def test_duplicate_genes_are_rejected(self):
        (self.base / "counts.csv").write_text("gene,c1,c2\nGENE_A,1,2\nGENE_A,3,4\n", encoding="utf-8")
        with tempfile.TemporaryDirectory(dir=TEST_TMP) as d:
            with self.assertRaisesRegex(ValueError, "基因 ID 重复"):
                runner.load_input(self.config(), Path(d))

    def test_duplicate_csv_headers_are_rejected_before_renaming(self):
        (self.base / "counts.csv").write_text('gene,"cell_1@a","cell_1@a"\nHK1,1,2\n', encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "列 ID 重复"):
            runner.load_input(self.config(), self.base)

    def test_unknown_configuration_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "scfea_epohs"):
            self.config(scfea_epohs=7)

    def test_output_values_and_feature_contract(self):
        for bad_value in (np.nan, np.inf):
            with self.assertRaisesRegex(ValueError, "非有限"):
                runner.validate_scores(pd.DataFrame({"c1": [bad_value]}, index=["M_1"]), ["c1"], ["M_1"])
        with self.assertRaisesRegex(ValueError, "特征"):
            runner.validate_scores(pd.DataFrame({"c1": [1.]}, index=["wrong"]), ["c1"], ["M_1"])

    def test_existing_run_is_preserved_before_preflight(self):
        out = self.base / "out"
        out.mkdir()
        (out / "run_manifest.json").write_text('{"original":true}', encoding="utf-8")
        (out / "preflight.json").write_text('{"original":true}', encoding="utf-8")
        self.config()
        with mock.patch("sys.argv", ["run.py", "run", "--config", str(self.base / "config.json")]):
            with self.assertRaises(FileExistsError):
                runner.main()
        self.assertEqual(json.loads((out / "preflight.json").read_text()), {"original": True})

    def test_summary_uses_descriptive_sample_counts(self):
        scores = pd.DataFrame({"c1": [1.0], "c2": [3.0], "c3": [8.0]}, index=["Glycolysis"])
        meta = pd.DataFrame({"condition": ["A", "A", "B"], "donor": ["d_1", "d_2", "d_3"]},
                            index=["c1", "c2", "c3"])
        out = self.base / "summary.csv"
        runner.summarize(scores, meta, {"group_by": ["condition"], "sample_column": "donor",
                                         "pathway_contains": None}, out)
        result = pd.read_csv(out)
        self.assertEqual(result.loc[result.condition == "A", "median_score"].iloc[0], 2.0)
        self.assertEqual(result.loc[result.condition == "A", "n_samples"].iloc[0], 2)
        self.assertNotIn("p_value", result.columns)

    def test_low_gene_coverage_is_flagged(self):
        pd.DataFrame({"c1": [1], "c2": [2]}, index=["UNMATCHED_GENE"]).to_csv(self.base / "counts.csv")
        c = self.config(methods=["gsva"], min_gene_overlap=0.1)
        with tempfile.TemporaryDirectory(dir=TEST_TMP) as d:
            frame, metadata = runner.load_input(c, Path(d))
        check, _ = runner.preflight(c, frame, metadata)
        self.assertFalse(check["ready"])
        self.assertEqual(check["overlap"]["pathway_genes_matched"], 0)
        self.assertTrue(any("覆盖率" in warning for warning in check["warnings"]))

    def test_h5ad_uses_gene_symbols_and_selects_before_dense_limit(self):
        try:
            import anndata as ad
        except ImportError:
            self.skipTest("anndata unavailable")
        a = ad.AnnData(X=np.zeros((4, 3)), obs=pd.DataFrame({"cell_type": ["T", "B", "T", "B"]},
                                                             index=["cell_1@a", "cell_2@b", "cell_3@c", "cell_4@d"]),
                       var=pd.DataFrame({"symbol": ["GENE_A", "GENE@B", "GENE_C"]},
                                        index=["ENSG1", "ENSG2", "ENSG3"]))
        a.layers["counts"] = np.array([[1, 2, 3], [0, 0, 0], [4, 5, 6], [0, 0, 0]])
        a.write_h5ad(self.base / "cells.h5ad")
        c = self.config(input="cells.h5ad", gene_column="symbol", max_dense_elements=6,
                        select={"column": "cell_type", "values": ["T"]})
        with tempfile.TemporaryDirectory(dir=TEST_TMP) as d:
            frame, meta = runner.load_input(c, Path(d))
        self.assertEqual(frame.shape, (3, 2))
        self.assertEqual(list(frame.index), ["GENE_A", "GENE@B", "GENE_C"])
        self.assertEqual(list(frame.columns), ["cell_1@a", "cell_3@c"])
        self.assertEqual(frame.loc["GENE_C", "cell_3@c"], 6)
        self.assertEqual(meta.loc["cell_3@c", "cell_type"], "T")

    def test_10x_matrix_and_external_selection(self):
        try:
            from scipy import sparse
            from scipy.io import mmwrite
        except ImportError:
            self.skipTest("scipy unavailable")
        folder = self.base / "tenx"
        folder.mkdir()
        mmwrite(folder / "matrix.mtx", sparse.csr_matrix([[1, 0, 3], [0, 2, 4]]))
        (folder / "features.tsv").write_text("id1\tGENE_A\tGene Expression\nid2\tGENE@B\tGene Expression\n", encoding="utf-8")
        (folder / "barcodes.tsv").write_text("cell_1@a\ncell_2@b\ncell_3@c\n", encoding="utf-8")
        pd.DataFrame({"cell_type": ["T", "B", "T"]}, index=["cell_1@a", "cell_2@b", "cell_3@c"]).to_csv(self.base / "meta.csv")
        c = self.config(input="tenx", metadata="meta.csv", max_dense_elements=4,
                        select={"column": "cell_type", "values": ["T"]})
        with tempfile.TemporaryDirectory(dir=TEST_TMP) as d:
            frame, _ = runner.load_input(c, Path(d))
        self.assertEqual(frame.shape, (2, 2))
        self.assertEqual(list(frame.columns), ["cell_1@a", "cell_3@c"])
        self.assertEqual(frame.loc["GENE@B", "cell_3@c"], 4)


if __name__ == "__main__":
    unittest.main()
