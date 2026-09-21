from pathlib import Path
import sqlite3

from data_validation import CELL_COLUMNS, load_validated_rows


ROOT = Path(__file__).resolve().parent
DATA_PATH = ROOT / "cell-count.csv"
SCHEMA_PATH = ROOT / "schema.sql"
DB_PATH = ROOT / "cell-count.db"


def create_database(rows, db_path=DB_PATH):
    db_path = Path(db_path)

    # The database is a generated artifact. Rebuilding from scratch
    # makes repeated pipeline runs deterministic.
    if db_path.exists():
        db_path.unlink()

    schema = SCHEMA_PATH.read_text(encoding="utf-8-sig")

    connection = sqlite3.connect(db_path)

    try:
        connection.execute("PRAGMA foreign_keys = ON")
        connection.executescript(schema)

        subjects = {}

        for row in rows:
            subjects[row["subject"]] = (
                row["subject"],
                row["project"],
                row["condition"],
                row["age"],
                row["sex"],
                row["treatment"],
                row["response"],
            )

        with connection:
            connection.executemany(
                """
                INSERT INTO subjects (
                    subject,
                    project,
                    condition,
                    age,
                    sex,
                    treatment,
                    response
                )
                VALUES (?, ?, ?, ?, ?, ?, ?)
                """,
                subjects.values(),
            )

            connection.executemany(
                """
                INSERT INTO samples (
                    sample,
                    subject,
                    sample_type,
                    time_from_treatment_start
                )
                VALUES (?, ?, ?, ?)
                """,
                (
                    (
                        row["sample"],
                        row["subject"],
                        row["sample_type"],
                        row["time_from_treatment_start"],
                    )
                    for row in rows
                ),
            )

            connection.executemany(
                """
                INSERT INTO cell_counts (
                    sample,
                    population,
                    count
                )
                VALUES (?, ?, ?)
                """,
                (
                    (
                        row["sample"],
                        population,
                        row[population],
                    )
                    for row in rows
                    for population in CELL_COLUMNS
                ),
            )

        foreign_key_errors = connection.execute(
            "PRAGMA foreign_key_check"
        ).fetchall()

        if foreign_key_errors:
            raise RuntimeError(
                f"foreign key check failed: {foreign_key_errors}"
            )

        integrity = connection.execute(
            "PRAGMA integrity_check"
        ).fetchone()[0]

        if integrity != "ok":
            raise RuntimeError(
                f"database integrity check failed: {integrity}"
            )

    except Exception:
        connection.close()

        if db_path.exists():
            db_path.unlink()

        raise

    else:
        connection.close()


def build_database(db_path=DB_PATH):
    rows = load_validated_rows(DATA_PATH)

    create_database(rows, db_path)


def table_count(db_path, table):
    connection = sqlite3.connect(db_path)

    try:
        return connection.execute(
            f"SELECT COUNT(*) FROM {table}"
        ).fetchone()[0]
    finally:
        connection.close()


def main():
    build_database(DB_PATH)

    print(f"Created {DB_PATH.name}")
    print(f"Subjects:    {table_count(DB_PATH, 'subjects'):,}")
    print(f"Samples:     {table_count(DB_PATH, 'samples'):,}")
    print(
        f"Cell counts: "
        f"{table_count(DB_PATH, 'cell_counts'):,}"
    )


if __name__ == "__main__":
    main()
