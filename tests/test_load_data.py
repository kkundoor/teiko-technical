import sqlite3
import tempfile
from contextlib import closing
import unittest
from pathlib import Path

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


class DatabaseLoadTests(unittest.TestCase):
    def create_test_database(self, rows):
        temp = tempfile.TemporaryDirectory()
        db_path = Path(temp.name) / "test.db"

        create_database(rows, db_path)

        self.addCleanup(temp.cleanup)

        return db_path

    def test_rows_are_normalized_into_relational_tables(self):
        rows = [
            sample_row(),
            sample_row(
                sample="sample2",
                time_from_treatment_start=7,
            ),
            sample_row(
                subject="subject2",
                sample="sample3",
                project="prj3",
                age=64,
                sex="F",
                response=None,
                sample_type="WB",
            ),
        ]

        db_path = self.create_test_database(rows)

        with closing(sqlite3.connect(db_path)) as connection:
            subjects = connection.execute(
                "SELECT COUNT(*) FROM subjects"
            ).fetchone()[0]

            samples = connection.execute(
                "SELECT COUNT(*) FROM samples"
            ).fetchone()[0]

            counts = connection.execute(
                "SELECT COUNT(*) FROM cell_counts"
            ).fetchone()[0]

        self.assertEqual(subjects, 2)
        self.assertEqual(samples, 3)
        self.assertEqual(counts, 15)

    def test_cell_counts_preserve_population_and_value(self):
        db_path = self.create_test_database([
            sample_row()
        ])

        with closing(sqlite3.connect(db_path)) as connection:
            measurements = connection.execute(
                """
                SELECT population, count
                FROM cell_counts
                WHERE sample = ?
                ORDER BY population
                """,
                ("sample1",),
            ).fetchall()

        self.assertEqual(
            measurements,
            [
                ("b_cell", 10),
                ("cd4_t_cell", 30),
                ("cd8_t_cell", 20),
                ("monocyte", 20),
                ("nk_cell", 20),
            ],
        )

    def test_null_response_is_preserved(self):
        db_path = self.create_test_database([
            sample_row(response=None)
        ])

        with closing(sqlite3.connect(db_path)) as connection:
            response = connection.execute(
                """
                SELECT response
                FROM subjects
                WHERE subject = ?
                """,
                ("subject1",),
            ).fetchone()[0]

        self.assertIsNone(response)

    def test_foreign_keys_are_valid(self):
        db_path = self.create_test_database([
            sample_row()
        ])

        with closing(sqlite3.connect(db_path)) as connection:
            connection.execute("PRAGMA foreign_keys = ON")

            errors = connection.execute(
                "PRAGMA foreign_key_check"
            ).fetchall()

        self.assertEqual(errors, [])

    def test_rebuilding_database_is_idempotent(self):
        rows = [sample_row()]

        with tempfile.TemporaryDirectory() as tmp:
            db_path = Path(tmp) / "test.db"

            create_database(rows, db_path)
            create_database(rows, db_path)

            with closing(sqlite3.connect(db_path)) as connection:
                samples = connection.execute(
                    "SELECT COUNT(*) FROM samples"
                ).fetchone()[0]

                counts = connection.execute(
                    "SELECT COUNT(*) FROM cell_counts"
                ).fetchone()[0]

        self.assertEqual(samples, 1)
        self.assertEqual(counts, 5)


if __name__ == "__main__":
    unittest.main()
