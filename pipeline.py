from analysis import (
    baseline_melanoma_miraclib_cohort,
    baseline_samples_by_project,
    baseline_subjects_by_response,
    baseline_subjects_by_sex,
    melanoma_male_responder_baseline_b_cell_summary,
    relative_frequencies,
    responder_statistics,
)


def main():
    relative = relative_frequencies()
    statistics = responder_statistics()
    baseline = baseline_melanoma_miraclib_cohort()
    projects = baseline_samples_by_project()
    responses = baseline_subjects_by_response()
    sexes = baseline_subjects_by_sex()
    b_cell_summary = melanoma_male_responder_baseline_b_cell_summary()

    print("Part 2")
    print(f"  Relative-frequency rows: {len(relative):,}")

    print("\nPart 3")
    print("  Responder vs. non-responder statistics:")

    for result in statistics:
        print(
            f"    {result['population']}: "
            f"p={result['p_value']:.6f}, "
            f"adjusted_p={result['q_value']:.6f}, "
            f"significant={result['significant']}"
        )

    print("\nPart 4")
    print(f"  Baseline samples: {len(baseline):,}")

    print(
        "  Samples by project: "
        + ", ".join(
            f"{row['project']}={row['sample_count']}"
            for row in projects
        )
    )

    print(
        "  Subjects by response: "
        + ", ".join(
            f"{row['response']}={row['subject_count']}"
            for row in responses
        )
    )

    print(
        "  Subjects by sex: "
        + ", ".join(
            f"{row['sex']}={row['subject_count']}"
            for row in sexes
        )
    )


    print(
        "\n  Melanoma male responders at time 0 "
        "(all treatments and sample types):"
    )

    print(
        "    Matching samples: "
        f"{b_cell_summary['sample_count']:,}"
    )

    average_b_cells = b_cell_summary[
        "average_b_cells"
    ]

    print(
        "    Average B cells: "
        + (
            f"{average_b_cells:.2f}"
            if average_b_cells is not None
            else "N/A"
        )
    )


if __name__ == "__main__":
    main()