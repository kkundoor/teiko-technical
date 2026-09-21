from pathlib import Path
import sqlite3
from collections import defaultdict
from statistics import median

from scipy.stats import mannwhitneyu


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


def melanoma_male_responder_baseline_b_cell_summary(
    db_path=DB_PATH,
):
    connection = sqlite3.connect(db_path)
    connection.row_factory = sqlite3.Row

    try:
        row = connection.execute(
            """
            SELECT
                COUNT(*) AS sample_count,
                AVG(c.count) AS average_b_cells
            FROM cell_counts AS c
            JOIN samples AS sa
                ON sa.sample = c.sample
            JOIN subjects AS s
                ON s.subject = sa.subject
            WHERE
                s.condition = 'melanoma'
                AND s.sex = 'M'
                AND s.response = 'yes'
                AND sa.time_from_treatment_start = 0
                AND c.population = 'b_cell'
            """
        ).fetchone()

        return {
            "sample_count": row["sample_count"],
            "average_b_cells": (
                float(row["average_b_cells"])
                if row["average_b_cells"] is not None
                else None
            ),
        }

    finally:
        connection.close()


def part3_subject_frequencies(db_path=DB_PATH):
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
            ),
            sample_frequencies AS (
                SELECT
                    s.subject,
                    s.response,
                    sa.sample,
                    sa.time_from_treatment_start,
                    c.population,
                    100.0 * c.count / t.total_count AS percentage
                FROM cell_counts AS c
                JOIN sample_totals AS t
                    ON t.sample = c.sample
                JOIN samples AS sa
                    ON sa.sample = c.sample
                JOIN subjects AS s
                    ON s.subject = sa.subject
                WHERE
                    s.condition = 'melanoma'
                    AND s.treatment = 'miraclib'
                    AND sa.sample_type = 'PBMC'
                    AND s.response IN ('yes', 'no')
            )
            SELECT
                subject,
                response,
                population,
                COUNT(*) AS sample_count,
                AVG(percentage) AS mean_percentage
            FROM sample_frequencies
            GROUP BY
                subject,
                response,
                population
            ORDER BY
                subject,
                CASE population
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


def benjamini_hochberg(p_values):
    count = len(p_values)

    if count == 0:
        return []

    order = sorted(
        range(count),
        key=lambda index: p_values[index],
    )

    adjusted = [0.0] * count
    running_min = 1.0

    for rank in range(count, 0, -1):
        index = order[rank - 1]

        candidate = (
            p_values[index]
            * count
            / rank
        )

        running_min = min(
            running_min,
            candidate,
            1.0,
        )

        adjusted[index] = running_min

    return adjusted


def responder_statistics(db_path=DB_PATH, alpha=0.05):
    rows = part3_subject_frequencies(db_path)

    grouped = defaultdict(
        lambda: {"yes": [], "no": []}
    )

    for row in rows:
        grouped[row["population"]][
            row["response"]
        ].append(row["mean_percentage"])

    results = []
    p_values = []

    for population in POPULATION_ORDER:
        responder = grouped[population]["yes"]
        nonresponder = grouped[population]["no"]

        if not responder or not nonresponder:
            raise ValueError(
                f"missing comparison group for {population}"
            )

        u_statistic, p_value = mannwhitneyu(
            responder,
            nonresponder,
            alternative="two-sided",
            method="asymptotic",
        )

        rank_biserial = (
            2.0 * u_statistic
            / (len(responder) * len(nonresponder))
            - 1.0
        )

        result = {
            "population": population,
            "responder_n": len(responder),
            "nonresponder_n": len(nonresponder),
            "responder_median": median(responder),
            "nonresponder_median": median(nonresponder),
            "median_difference": (
                median(responder)
                - median(nonresponder)
            ),
            "u_statistic": float(u_statistic),
            "p_value": float(p_value),
            "rank_biserial": float(rank_biserial),
        }

        results.append(result)
        p_values.append(float(p_value))

    q_values = benjamini_hochberg(p_values)

    for result, q_value in zip(results, q_values):
        result["q_value"] = q_value
        result["significant"] = q_value < alpha

    return results