import csv
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
DATA_PATH = ROOT / "cell-count.csv"


from data_validation import (
    COLUMNS,
    DataValidationError,
    load_validated_rows,
)


def base_row(**overrides):
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


def write_csv(path, rows, columns=None):
    columns = columns or COLUMNS

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


class SourceValidationTests(unittest.TestCase):
    def validate_rows(self, rows, columns=None):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "input.csv"
            write_csv(path, rows, columns)
            return load_validated_rows(path)

    def test_valid_rows_are_parsed_and_blank_response_is_null(self):
        rows = [
            base_row(
                response="",
                time_from_treatment_start=-7,
            ),
            base_row(
                sample="sample2",
                sample_type="WB",
                response="",
                time_from_treatment_start=0,
            ),
        ]

        parsed = self.validate_rows(rows)

        self.assertEqual(len(parsed), 2)
        self.assertIsNone(parsed[0]["response"])
        self.assertEqual(parsed[0]["age"], 57)
        self.assertEqual(parsed[0]["b_cell"], 10)
        self.assertEqual(
            parsed[0]["time_from_treatment_start"], -7
        )

    def test_duplicate_sample_id_is_rejected(self):
        rows = [
            base_row(),
            base_row(time_from_treatment_start=7),
        ]

        with self.assertRaisesRegex(
            DataValidationError, "duplicate sample id"
        ):
            self.validate_rows(rows)

    def test_inconsistent_subject_metadata_is_rejected(self):
        rows = [
            base_row(),
            base_row(
                sample="sample2",
                age=58,
                time_from_treatment_start=7,
            ),
        ]

        with self.assertRaisesRegex(
            DataValidationError, "inconsistent metadata"
        ):
            self.validate_rows(rows)

    def test_negative_cell_count_is_rejected(self):
        with self.assertRaisesRegex(
            DataValidationError, "cannot be negative"
        ):
            self.validate_rows([
                base_row(b_cell=-1)
            ])

    def test_zero_total_count_is_rejected(self):
        with self.assertRaisesRegex(
            DataValidationError,
            "cell-count total must be greater than zero",
        ):
            self.validate_rows([
                base_row(
                    b_cell=0,
                    cd8_t_cell=0,
                    cd4_t_cell=0,
                    nk_cell=0,
                    monocyte=0,
                )
            ])

    def test_schema_drift_is_rejected(self):
        columns = COLUMNS[:-1]

        row = base_row()
        row.pop("monocyte")

        with self.assertRaisesRegex(
            DataValidationError, "missing columns: monocyte"
        ):
            self.validate_rows([row], columns)

    def test_provided_dataset_passes_validation(self):
        rows = load_validated_rows(
            DATA_PATH
        )

        self.assertGreater(len(rows), 0)
        self.assertEqual(
            len(rows),
            len({row["sample"] for row in rows}),
        )


if __name__ == "__main__":
    unittest.main()
