import tempfile
import unittest
from pathlib import Path

from analysis import (
    baseline_melanoma_miraclib_cohort,
    baseline_samples_by_project,
    baseline_subjects_by_response,
    baseline_subjects_by_sex,
    relative_frequencies,
)
from load_data import create_database


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


if __name__ == "__main__":
    unittest.main()