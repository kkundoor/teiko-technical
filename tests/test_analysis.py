import tempfile
import unittest
from pathlib import Path

from analysis import (
    baseline_melanoma_miraclib_cohort,
    baseline_samples_by_project,
    baseline_subjects_by_response,
    baseline_subjects_by_sex,
    benjamini_hochberg,
    melanoma_male_responder_baseline_b_cell_summary,
    part3_subject_frequencies,
    relative_frequencies,
    responder_statistics,
)
from data_validation import load_validated_rows
from load_data import create_database


ROOT = Path(__file__).resolve().parent.parent
DATA_PATH = ROOT / "cell-count.csv"


def sample_row(**overrides):
    row = {
        "project": "prj1",
        "subject": "subject1",
        "condition": "melanoma",
        "age": 57,
        "sex": "M",
        "treatment": "miraclib",
        "response": "yes",
        "sample": "sample1",
        "sample_type": "PBMC",
        "time_from_treatment_start": 0,
        "b_cell": 10,
        "cd8_t_cell": 20,
        "cd4_t_cell": 30,
        "nk_cell": 20,
        "monocyte": 20,
    }
    row.update(overrides)
    return row


class AnalysisTestCase(unittest.TestCase):
    def create_test_database(self, rows):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)

        db_path = Path(temp.name) / "test.db"
        create_database(rows, db_path)

        return db_path


class RelativeFrequencyTests(AnalysisTestCase):
    def test_relative_frequency_output_shape_and_values(self):
        db_path = self.create_test_database([
            sample_row()
        ])

        rows = relative_frequencies(db_path)

        self.assertEqual(len(rows), 5)

        self.assertEqual(
            list(rows[0]),
            [
                "sample",
                "total_count",
                "population",
                "count",
                "percentage",
            ],
        )

        expected = {
            "b_cell": (10, 10.0),
            "cd8_t_cell": (20, 20.0),
            "cd4_t_cell": (30, 30.0),
            "nk_cell": (20, 20.0),
            "monocyte": (20, 20.0),
        }

        for row in rows:
            self.assertEqual(row["sample"], "sample1")
            self.assertEqual(row["total_count"], 100)

            count, percentage = expected[row["population"]]

            self.assertEqual(row["count"], count)
            self.assertAlmostEqual(
                row["percentage"],
                percentage,
            )

    def test_each_sample_has_five_population_rows(self):
        db_path = self.create_test_database([
            sample_row(),
            sample_row(
                subject="subject2",
                sample="sample2",
                b_cell=5,
                cd8_t_cell=5,
                cd4_t_cell=5,
                nk_cell=5,
                monocyte=5,
            ),
        ])

        rows = relative_frequencies(db_path)

        populations_by_sample = {}

        for row in rows:
            populations_by_sample.setdefault(
                row["sample"], []
            ).append(row)

        self.assertEqual(
            set(populations_by_sample),
            {"sample1", "sample2"},
        )

        for sample_rows in populations_by_sample.values():
            self.assertEqual(len(sample_rows), 5)

    def test_percentages_sum_to_one_hundred_per_sample(self):
        db_path = self.create_test_database([
            sample_row(),
            sample_row(
                subject="subject2",
                sample="sample2",
                b_cell=1,
                cd8_t_cell=1,
                cd4_t_cell=1,
                nk_cell=1,
                monocyte=2,
            ),
        ])

        rows = relative_frequencies(db_path)

        totals = {}

        for row in rows:
            totals.setdefault(row["sample"], 0.0)
            totals[row["sample"]] += row["percentage"]

        for total in totals.values():
            self.assertAlmostEqual(
                total,
                100.0,
                places=10,
            )

    def test_total_count_is_repeated_for_each_population(self):
        db_path = self.create_test_database([
            sample_row(
                b_cell=3,
                cd8_t_cell=7,
                cd4_t_cell=11,
                nk_cell=13,
                monocyte=17,
            )
        ])

        rows = relative_frequencies(db_path)

        self.assertEqual(
            {row["total_count"] for row in rows},
            {51},
        )


class BaselineCohortTests(AnalysisTestCase):
    def test_baseline_cohort_applies_all_four_filters(self):
        rows = [
            sample_row(
                sample="include1",
                subject="subject1",
                project="prj1",
                response="yes",
                sex="M",
            ),
            sample_row(
                sample="include2",
                subject="subject2",
                project="prj3",
                response="no",
                sex="F",
            ),
            sample_row(
                sample="wrong_time",
                subject="subject3",
                time_from_treatment_start=7,
            ),
            sample_row(
                sample="wrong_type",
                subject="subject4",
                sample_type="WB",
            ),
            sample_row(
                sample="wrong_treatment",
                subject="subject5",
                treatment="phauximab",
            ),
            sample_row(
                sample="wrong_condition",
                subject="subject6",
                condition="carcinoma",
            ),
        ]

        db_path = self.create_test_database(rows)

        cohort = baseline_melanoma_miraclib_cohort(db_path)

        self.assertEqual(
            {row["sample"] for row in cohort},
            {"include1", "include2"},
        )

    def test_project_summary_counts_samples(self):
        rows = [
            sample_row(
                sample="sample1",
                subject="subject1",
                project="prj1",
            ),
            sample_row(
                sample="sample1b",
                subject="subject1",
                project="prj1",
            ),
            sample_row(
                sample="sample2",
                subject="subject2",
                project="prj3",
            ),
        ]

        db_path = self.create_test_database(rows)

        summary = baseline_samples_by_project(db_path)

        self.assertEqual(
            summary,
            [
                {"project": "prj1", "sample_count": 2},
                {"project": "prj3", "sample_count": 1},
            ],
        )

    def test_response_summary_counts_distinct_subjects(self):
        rows = [
            sample_row(
                sample="sample1",
                subject="subject1",
                response="yes",
            ),
            sample_row(
                sample="sample1b",
                subject="subject1",
                response="yes",
            ),
            sample_row(
                sample="sample2",
                subject="subject2",
                response="no",
            ),
            sample_row(
                sample="sample3",
                subject="subject3",
                response="yes",
            ),
        ]

        db_path = self.create_test_database(rows)

        summary = baseline_subjects_by_response(db_path)

        self.assertEqual(
            summary,
            [
                {"response": "no", "subject_count": 1},
                {"response": "yes", "subject_count": 2},
            ],
        )

    def test_sex_summary_counts_distinct_subjects(self):
        rows = [
            sample_row(
                sample="sample1",
                subject="subject1",
                sex="M",
            ),
            sample_row(
                sample="sample1b",
                subject="subject1",
                sex="M",
            ),
            sample_row(
                sample="sample2",
                subject="subject2",
                sex="F",
            ),
            sample_row(
                sample="sample3",
                subject="subject3",
                sex="M",
            ),
        ]

        db_path = self.create_test_database(rows)

        summary = baseline_subjects_by_sex(db_path)

        self.assertEqual(
            summary,
            [
                {"sex": "F", "subject_count": 1},
                {"sex": "M", "subject_count": 2},
            ],
        )



class Part4SubsetTests(AnalysisTestCase):
    def test_b_cell_summary_uses_all_treatments_and_sample_types(self):
        rows = [
            sample_row(
                sample="included_pbmc",
                subject="included_pbmc",
                treatment="miraclib",
                sample_type="PBMC",
                b_cell=100,
            ),
            sample_row(
                sample="included_wb",
                subject="included_wb",
                treatment="phauximab",
                sample_type="WB",
                b_cell=200,
            ),
            sample_row(
                sample="wrong_sex",
                subject="wrong_sex",
                sex="F",
                b_cell=900,
            ),
            sample_row(
                sample="wrong_response",
                subject="wrong_response",
                response="no",
                b_cell=900,
            ),
            sample_row(
                sample="wrong_time",
                subject="wrong_time",
                time_from_treatment_start=7,
                b_cell=900,
            ),
            sample_row(
                sample="wrong_condition",
                subject="wrong_condition",
                condition="carcinoma",
                b_cell=900,
            ),
        ]

        db_path = self.create_test_database(rows)

        result = (
            melanoma_male_responder_baseline_b_cell_summary(
                db_path
            )
        )

        self.assertEqual(
            result["sample_count"],
            2,
        )

        self.assertAlmostEqual(
            result["average_b_cells"],
            150.0,
        )

    def test_b_cell_summary_handles_empty_subset(self):
        db_path = self.create_test_database([
            sample_row(
                sample="nonmatching",
                subject="nonmatching",
                condition="carcinoma",
            )
        ])

        result = (
            melanoma_male_responder_baseline_b_cell_summary(
                db_path
            )
        )

        self.assertEqual(
            result["sample_count"],
            0,
        )

        self.assertIsNone(
            result["average_b_cells"]
        )

    def test_b_cell_summary_matches_provided_dataset(self):
        rows = load_validated_rows(
            DATA_PATH
        )

        matching = [
            row
            for row in rows
            if (
                row["condition"] == "melanoma"
                and row["sex"] == "M"
                and row["response"] == "yes"
                and row["time_from_treatment_start"] == 0
            )
        ]

        self.assertGreater(
            len(matching),
            0,
        )

        expected_average = (
            sum(
                row["b_cell"]
                for row in matching
            )
            / len(matching)
        )

        db_path = self.create_test_database(
            rows
        )

        result = (
            melanoma_male_responder_baseline_b_cell_summary(
                db_path
            )
        )

        self.assertEqual(
            result["sample_count"],
            len(matching),
        )

        self.assertAlmostEqual(
            result["average_b_cells"],
            expected_average,
        )


class StatisticalAnalysisTests(AnalysisTestCase):
    def test_subject_frequencies_average_repeated_samples(self):
        rows = [
            sample_row(
                sample="yes_day0",
                subject="subject1",
                response="yes",
                time_from_treatment_start=0,
                b_cell=60,
                cd8_t_cell=10,
                cd4_t_cell=10,
                nk_cell=10,
                monocyte=10,
            ),
            sample_row(
                sample="yes_day7",
                subject="subject1",
                response="yes",
                time_from_treatment_start=7,
                b_cell=40,
                cd8_t_cell=15,
                cd4_t_cell=15,
                nk_cell=15,
                monocyte=15,
            ),
        ]

        db_path = self.create_test_database(rows)

        values = part3_subject_frequencies(db_path)

        b_cell = next(
            row
            for row in values
            if row["population"] == "b_cell"
        )

        self.assertEqual(b_cell["sample_count"], 2)
        self.assertAlmostEqual(
            b_cell["mean_percentage"],
            50.0,
        )


    def test_part3_cohort_includes_repeated_timepoints_and_filters_scope(self):
        rows = [
            sample_row(
                sample="included_day0",
                subject="included",
                response="yes",
                time_from_treatment_start=0,
            ),
            sample_row(
                sample="included_day7",
                subject="included",
                response="yes",
                time_from_treatment_start=7,
            ),
            sample_row(
                sample="included_day14",
                subject="included",
                response="yes",
                time_from_treatment_start=14,
            ),
            sample_row(
                sample="wrong_condition",
                subject="wrong_condition",
                condition="carcinoma",
            ),
            sample_row(
                sample="wrong_treatment",
                subject="wrong_treatment",
                treatment="phauximab",
            ),
            sample_row(
                sample="wrong_type",
                subject="wrong_type",
                sample_type="WB",
            ),
            sample_row(
                sample="missing_response",
                subject="missing_response",
                response=None,
            ),
        ]

        db_path = self.create_test_database(rows)

        values = part3_subject_frequencies(db_path)

        subjects = {
            row["subject"]
            for row in values
        }

        self.assertEqual(
            subjects,
            {"included"},
        )

        self.assertEqual(
            len(values),
            5,
        )

        self.assertEqual(
            {
                row["sample_count"]
                for row in values
            },
            {3},
        )
    def test_benjamini_hochberg_adjustment(self):
        adjusted = benjamini_hochberg([
            0.01,
            0.04,
            0.03,
            0.002,
        ])

        expected = [
            0.02,
            0.04,
            0.04,
            0.008,
        ]

        for actual, expected_value in zip(
            adjusted,
            expected,
        ):
            self.assertAlmostEqual(
                actual,
                expected_value,
            )

    def test_responder_statistics_count_subjects_not_samples(self):
        rows = [
            sample_row(
                sample="yes1_day0",
                subject="yes1",
                response="yes",
                b_cell=60,
                cd8_t_cell=10,
                cd4_t_cell=10,
                nk_cell=10,
                monocyte=10,
            ),
            sample_row(
                sample="yes1_day7",
                subject="yes1",
                response="yes",
                time_from_treatment_start=7,
                b_cell=70,
                cd8_t_cell=8,
                cd4_t_cell=8,
                nk_cell=7,
                monocyte=7,
            ),
            sample_row(
                sample="yes2",
                subject="yes2",
                response="yes",
                b_cell=65,
                cd8_t_cell=9,
                cd4_t_cell=9,
                nk_cell=9,
                monocyte=8,
            ),
            sample_row(
                sample="no1_day0",
                subject="no1",
                response="no",
                b_cell=20,
                cd8_t_cell=20,
                cd4_t_cell=20,
                nk_cell=20,
                monocyte=20,
            ),
            sample_row(
                sample="no1_day7",
                subject="no1",
                response="no",
                time_from_treatment_start=7,
                b_cell=10,
                cd8_t_cell=23,
                cd4_t_cell=23,
                nk_cell=22,
                monocyte=22,
            ),
            sample_row(
                sample="no2",
                subject="no2",
                response="no",
                b_cell=15,
                cd8_t_cell=22,
                cd4_t_cell=21,
                nk_cell=21,
                monocyte=21,
            ),
        ]

        db_path = self.create_test_database(rows)

        results = responder_statistics(db_path)

        b_cell = next(
            row
            for row in results
            if row["population"] == "b_cell"
        )

        self.assertEqual(
            b_cell["responder_n"],
            2,
        )
        self.assertEqual(
            b_cell["nonresponder_n"],
            2,
        )
        self.assertGreater(
            b_cell["rank_biserial"],
            0,
        )

if __name__ == "__main__":
    unittest.main()