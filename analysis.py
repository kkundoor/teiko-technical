from pathlib import Path
import sqlite3


ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "cell-count.db"

POPULATION_ORDER = [
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
]


def relative_frequencies(db_path=DB_PATH):
    connection = sqlite3.connect(db_path)
    connection.row_factory = sqlite3.Row

    try:
        rows = connection.execute(
            """
            WITH sample_totals AS (
                SELECT
                    sample,
                    SUM(count) AS total_count
                FROM cell_counts
                GROUP BY sample
            )
            SELECT
                c.sample AS sample,
                t.total_count AS total_count,
                c.population AS population,
                c.count AS count,
                100.0 * c.count / t.total_count AS percentage
            FROM cell_counts AS c
            JOIN sample_totals AS t
                ON t.sample = c.sample
            ORDER BY
                c.sample,
                CASE c.population
                    WHEN 'b_cell' THEN 1
                    WHEN 'cd8_t_cell' THEN 2
                    WHEN 'cd4_t_cell' THEN 3
                    WHEN 'nk_cell' THEN 4
                    WHEN 'monocyte' THEN 5
                END
            """
        ).fetchall()

        return [dict(row) for row in rows]

    finally:
        connection.close()


def baseline_melanoma_miraclib_cohort(db_path=DB_PATH):
    connection = sqlite3.connect(db_path)
    connection.row_factory = sqlite3.Row

    try:
        rows = connection.execute(
            """
            SELECT
                sa.sample,
                sa.subject,
                s.project,
                s.response,
                s.sex
            FROM samples AS sa
            JOIN subjects AS s
                ON s.subject = sa.subject
            WHERE
                s.condition = 'melanoma'
                AND s.treatment = 'miraclib'
                AND sa.sample_type = 'PBMC'
                AND sa.time_from_treatment_start = 0
            ORDER BY sa.sample
            """
        ).fetchall()

        return [dict(row) for row in rows]

    finally:
        connection.close()


def baseline_samples_by_project(db_path=DB_PATH):
    connection = sqlite3.connect(db_path)
    connection.row_factory = sqlite3.Row

    try:
        rows = connection.execute(
            """
            SELECT
                s.project,
                COUNT(sa.sample) AS sample_count
            FROM samples AS sa
            JOIN subjects AS s
                ON s.subject = sa.subject
            WHERE
                s.condition = 'melanoma'
                AND s.treatment = 'miraclib'
                AND sa.sample_type = 'PBMC'
                AND sa.time_from_treatment_start = 0
            GROUP BY s.project
            ORDER BY s.project
            """
        ).fetchall()

        return [dict(row) for row in rows]

    finally:
        connection.close()


def baseline_subjects_by_response(db_path=DB_PATH):
    connection = sqlite3.connect(db_path)
    connection.row_factory = sqlite3.Row

    try:
        rows = connection.execute(
            """
            SELECT
                s.response,
                COUNT(DISTINCT s.subject) AS subject_count
            FROM samples AS sa
            JOIN subjects AS s
                ON s.subject = sa.subject
            WHERE
                s.condition = 'melanoma'
                AND s.treatment = 'miraclib'
                AND sa.sample_type = 'PBMC'
                AND sa.time_from_treatment_start = 0
            GROUP BY s.response
            ORDER BY s.response
            """
        ).fetchall()

        return [dict(row) for row in rows]

    finally:
        connection.close()


def baseline_subjects_by_sex(db_path=DB_PATH):
    connection = sqlite3.connect(db_path)
    connection.row_factory = sqlite3.Row

    try:
        rows = connection.execute(
            """
            SELECT
                s.sex,
                COUNT(DISTINCT s.subject) AS subject_count
            FROM samples AS sa
            JOIN subjects AS s
                ON s.subject = sa.subject
            WHERE
                s.condition = 'melanoma'
                AND s.treatment = 'miraclib'
                AND sa.sample_type = 'PBMC'
                AND sa.time_from_treatment_start = 0
            GROUP BY s.sex
            ORDER BY s.sex
            """
        ).fetchall()

        return [dict(row) for row in rows]

    finally:
        connection.close()