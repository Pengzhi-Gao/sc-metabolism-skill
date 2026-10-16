import importlib.util
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

module_path = Path(__file__).resolve().parents[1] / "examples" / "summarize_results.py"
module_spec = importlib.util.spec_from_file_location("metabolism_example", module_path)
example = importlib.util.module_from_spec(module_spec)
module_spec.loader.exec_module(example)


class ResultExamples(unittest.TestCase):
    def test_donors_receive_one_aggregate_each(self):
        scores = pd.DataFrame([[100., 100., 100., 0., 0.]], index=["path"], columns=list("abcde"))
        meta = pd.DataFrame({"donor": ["A", "A", "A", "B", "C"]}, index=scores.columns)
        values, units, excluded = example.sample_medians(scores, meta, ["donor"], 1)
        self.assertEqual(values.loc["path"].median(), 0.)
        self.assertEqual(list(units.n_cells), [3, 1, 1])
        self.assertTrue(excluded.empty)

    def test_paired_comparison_aligns_donors_and_reports_exclusions(self):
        units = pd.DataFrame({"donor": ["B", "A", "A", "B", "C", "C", "D"],
                              "group": ["T", "C", "T", "C", "C", "T", "T"]},
                             index=[f"u{x}" for x in range(7)])
        values = pd.DataFrame([[6., 1., 3., 2., 3., 9., 999.], [1.] * 7],
                              index=["change", "constant"], columns=units.index)
        spec = {"sample": "donor", "condition": "group", "case": "T", "control": "C",
                "design": "paired", "strata": [], "min_samples": 3}
        result = example.compare_samples(values, units, spec).set_index("feature")
        self.assertEqual(result.loc["change", "effect_case_minus_control"], 4.)
        self.assertEqual(result.loc["change", "n_pairs"], 3)
        self.assertEqual(result.loc["change", "n_unpaired_excluded"], 1)
        self.assertEqual(result.loc["constant", "p_value"], 1.)
        spec["min_samples"] = 4
        low = example.compare_samples(values, units, spec)
        self.assertTrue(low.p_value.isna().all())
        spec["design"] = "unpaired"
        with self.assertRaisesRegex(ValueError, "跨组相同样本"):
            example.compare_samples(values, units, spec)

    def test_bh_uses_tested_features_and_preserves_missing(self):
        actual = example.bh_adjust([.01, .04, .03, np.nan])
        np.testing.assert_allclose(actual[:3], [.03, .04, .04])
        self.assertTrue(np.isnan(actual[3]))

    def test_unpaired_effect_uses_independent_sample_scores(self):
        units = pd.DataFrame({"donor": list("ABCDEF"), "group": ["T"] * 3 + ["C"] * 3}, index=list("abcdef"))
        values = pd.DataFrame([[7, 8, 9, 1, 2, 3]], index=["path"], columns=units.index)
        spec = {"sample": "donor", "condition": "group", "case": "T", "control": "C",
                "design": "unpaired", "strata": [], "min_samples": 3}
        row = example.compare_samples(values, units, spec).iloc[0]
        self.assertEqual(row.effect_case_minus_control, 6.)
        self.assertEqual(row.status, "tested")
        self.assertTrue(np.isfinite(row.p_value))


if __name__ == "__main__":
    unittest.main()
