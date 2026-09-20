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
