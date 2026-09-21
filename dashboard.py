from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

from analysis import (
    POPULATION_ORDER,
    baseline_melanoma_miraclib_cohort,
    baseline_samples_by_project,
    baseline_subjects_by_response,
    baseline_subjects_by_sex,
    melanoma_male_responder_baseline_b_cell_summary,
    part3_subject_frequencies,
    relative_frequencies,
    responder_statistics,
)
from load_data import build_database


ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "cell-count.db"

POPULATION_LABELS = {
    "b_cell": "B cells",
    "cd8_t_cell": "CD8 T cells",
    "cd4_t_cell": "CD4 T cells",
    "nk_cell": "NK cells",
    "monocyte": "Monocytes",
}

POPULATION_LABEL_ORDER = [
    POPULATION_LABELS[population]
    for population in POPULATION_ORDER
]


st.set_page_config(
    page_title="Immune Cell Population Analysis",
    layout="wide",
)


if not DB_PATH.exists():
    try:
        with st.spinner(
            "Building SQLite database from cell-count.csv..."
        ):
            build_database(DB_PATH)
    except Exception as error:
        st.error(
            "Could not build the SQLite database from "
            "cell-count.csv."
        )
        st.exception(error)
        st.stop()


db_version = DB_PATH.stat().st_mtime_ns


@st.cache_data
def load_relative_frequencies(db_version):
    return pd.DataFrame(relative_frequencies())


@st.cache_data
def load_subject_frequencies(db_version):
    return pd.DataFrame(part3_subject_frequencies())


@st.cache_data
def load_responder_statistics(db_version):
    return pd.DataFrame(responder_statistics())


@st.cache_data
def load_baseline_cohort(db_version):
    return pd.DataFrame(
        baseline_melanoma_miraclib_cohort()
    )


@st.cache_data
def load_project_summary(db_version):
    return pd.DataFrame(
        baseline_samples_by_project()
    )


@st.cache_data
def load_response_summary(db_version):
    return pd.DataFrame(
        baseline_subjects_by_response()
    )


@st.cache_data
def load_sex_summary(db_version):
    return pd.DataFrame(
        baseline_subjects_by_sex()
    )


@st.cache_data
def load_b_cell_summary(db_version):
    return melanoma_male_responder_baseline_b_cell_summary()


st.title("Immune Cell Population Analysis")

st.caption(
    "Results are queried from `cell-count.db`, generated from "
    "the supplied CSV by `load_data.py`."
)

part2_tab, part3_tab, part4_tab = st.tabs(
    [
        "Part 2: Relative Cell Frequencies",
        "Part 3: Responder Comparison",
        "Part 4: Baseline and Subset Analysis",
    ]
)


# ============================================================
# PART 2
# ============================================================

with part2_tab:
    st.header("Relative Cell Frequencies")

    st.write(
        "For each sample, the five measured cell populations "
        "are summed to obtain the total count. Each population "
        "is shown as a percentage of that total."
    )

    relative = load_relative_frequencies(db_version)

    sample_ids = relative["sample"].drop_duplicates().tolist()

    selected_sample = st.selectbox(
        "Sample",
        sample_ids,
        index=0,
    )

    sample_data = relative[
        relative["sample"] == selected_sample
    ].copy()

    sample_data["population_label"] = (
        sample_data["population"].map(
            POPULATION_LABELS
        )
    )

    total_count = int(
        sample_data["total_count"].iloc[0]
    )

    metric_left, metric_right = st.columns(2)

    metric_left.metric(
        "Selected sample",
        selected_sample,
    )

    metric_right.metric(
        "Total measured cells",
        f"{total_count:,}",
    )

    figure = px.bar(
        sample_data,
        x="population_label",
        y="percentage",
        hover_data={
            "count": ":,",
            "percentage": ":.3f",
            "population_label": False,
        },
        labels={
            "population_label": "Population",
            "percentage": "Relative frequency (%)",
            "count": "Cell count",
        },
        category_orders={
            "population_label": POPULATION_LABEL_ORDER
        },
        color_discrete_sequence=["#6FA3D2"],
    )

    figure.update_layout(
        showlegend=False,
        xaxis_title=None,
    )

    st.plotly_chart(
        figure,
        width="stretch",
    )

    display_sample = sample_data[
        [
            "population_label",
            "count",
            "percentage",
        ]
    ].copy()

    display_sample.columns = [
        "Population",
        "Count",
        "Percentage",
    ]

    display_sample["Percentage"] = (
        display_sample["Percentage"].round(3)
    )

    st.dataframe(
        display_sample,
        hide_index=True,
        width="stretch",
    )

    with st.expander("View Complete Relative-Frequency Summary"):
        st.caption(
            f"{len(relative):,} population-level rows across "
            f"{relative['sample'].nunique():,} samples."
        )

        full_summary = relative[
            [
                "sample",
                "total_count",
                "population",
                "count",
                "percentage",
            ]
        ].copy()

        st.dataframe(
            full_summary,
            hide_index=True,
            width="stretch",
            column_config={
                "percentage": st.column_config.NumberColumn(
                    "percentage",
                    format="%.4f",
                ),
            },
        )


# ============================================================
# PART 3
# ============================================================

with part3_tab:
    st.header(
        "Responder vs. Non-Responder Comparison"
    )

    st.write(
        "Cohort: melanoma patients receiving miraclib, using "
        "PBMC samples only. Each subject contributes one mean "
        "relative frequency per population across days 0, 7, "
        "and 14."
    )

    subject_frequencies = load_subject_frequencies(
        db_version
    )

    statistics = load_responder_statistics(
        db_version
    )

    subject_frequencies["Population"] = (
        subject_frequencies["population"].map(
            POPULATION_LABELS
        )
    )

    subject_frequencies["Response"] = (
        subject_frequencies["response"].map(
            {
                "yes": "Responder",
                "no": "Non-responder",
            }
        )
    )

    independent_subjects = (
        subject_frequencies["subject"].nunique()
    )

    responder_count = (
        subject_frequencies.loc[
            subject_frequencies["response"] == "yes",
            "subject",
        ].nunique()
    )

    nonresponder_count = (
        subject_frequencies.loc[
            subject_frequencies["response"] == "no",
            "subject",
        ].nunique()
    )

    significant_count = int(
        statistics["significant"].sum()
    )

    col1, col2, col3, col4 = st.columns(4)

    col1.metric("Subjects", f"{independent_subjects:,}")
    col2.metric("Responders", f"{responder_count:,}")
    col3.metric(
        "Non-responders",
        f"{nonresponder_count:,}",
    )
    col4.metric(
        "BH-Adjusted p < 0.05",
        significant_count,
    )

    if significant_count == 0:
        st.info(
            "No population remains statistically significant "
            "after Benjamini-Hochberg correction at FDR 0.05."
        )

    boxplot = px.box(
        subject_frequencies,
        x="Population",
        y="mean_percentage",
        color="Response",
        color_discrete_map={
            "Responder": "#5B8E7D",
            "Non-responder": "#7C8796",
        },
        category_orders={
            "Population": POPULATION_LABEL_ORDER,
            "Response": [
                "Responder",
                "Non-responder",
            ],
        },
        labels={
            "mean_percentage":
                "Mean relative frequency per subject (%)",
        },
        points=False,
    )

    boxplot.update_layout(
        boxmode="group",
        xaxis_title=None,
        legend_title_text="",
    )

    st.plotly_chart(
        boxplot,
        width="stretch",
    )

    display_statistics = statistics.copy()

    display_statistics["Population"] = (
        display_statistics["population"].map(
            POPULATION_LABELS
        )
    )

    display_statistics = display_statistics[
        [
            "Population",
            "responder_n",
            "nonresponder_n",
            "responder_median",
            "nonresponder_median",
            "median_difference",
            "p_value",
            "q_value",
            "rank_biserial",
            "significant",
        ]
    ]

    display_statistics.columns = [
        "Population",
        "Responder n",
        "Non-responder n",
        "Responder median (%)",
        "Non-responder median (%)",
        "Difference (pp)",
        "Raw p",
        "BH-adjusted p",
        "Rank-biserial correlation",
        "Adjusted p < 0.05",
    ]

    numeric_columns = [
        "Responder median (%)",
        "Non-responder median (%)",
        "Difference (pp)",
        "Raw p",
        "BH-adjusted p",
        "Rank-biserial correlation",
    ]

    display_statistics[numeric_columns] = (
        display_statistics[numeric_columns].round(4)
    )

    st.dataframe(
        display_statistics,
        hide_index=True,
        width="stretch",
    )

    strongest = statistics.loc[
        statistics["q_value"].idxmin()
    ]

    st.caption(
        "Smallest BH-adjusted p-value: "
        f"{POPULATION_LABELS[strongest['population']]} "
        f"(raw p={strongest['p_value']:.4f}, "
        f"adjusted p={strongest['q_value']:.4f}, "
        f"median difference="
        f"{strongest['median_difference']:+.3f} pp)."
    )

    with st.expander("Statistical Method"):
        st.write(
            "Each subject contributes one value per population: "
            "the mean relative frequency across that subject's "
            "days 0, 7, and 14 samples. Responders and "
            "non-responders are compared with a two-sided "
            "Mann-Whitney U test. The five population-level "
            "p-values are adjusted using Benjamini-Hochberg FDR "
            "correction at 0.05. Rank-biserial correlation is "
            "reported as an effect-size measure. This evaluates "
            "an overall responder/non-responder difference across "
            "the observed period; it does not estimate a "
            "time-specific treatment trajectory or out-of-sample "
            "predictive performance."
        )


# ============================================================
# PART 4
# ============================================================

with part4_tab:
    st.header("Baseline and Subset Analysis")

    st.subheader(
        "Baseline Miraclib PBMC Cohort"
    )

    st.write(
        "Melanoma patients receiving miraclib, using PBMC "
        "samples at time from treatment start = 0."
    )

    cohort = load_baseline_cohort(db_version)
    projects = load_project_summary(db_version)
    responses = load_response_summary(db_version)
    sexes = load_sex_summary(db_version)
    b_cell_summary = load_b_cell_summary(db_version)

    st.metric(
        "Baseline cohort samples",
        f"{len(cohort):,}",
    )

    project_col, response_col, sex_col = st.columns(3)

    with project_col:
        st.subheader("Samples by Project")

        project_chart = px.bar(
            projects,
            x="project",
            y="sample_count",
            labels={
                "project": "Project",
                "sample_count": "Samples",
            },
            text_auto=True,
            color_discrete_sequence=["#6FA3D2"],
        )

        project_chart.update_layout(
            showlegend=False,
            xaxis_title=None,
        )

        st.plotly_chart(
            project_chart,
            width="stretch",
        )

    with response_col:
        st.subheader("Subjects by Response")

        response_display = responses.copy()

        response_display["response"] = (
            response_display["response"].map(
                {
                    "yes": "Responder",
                    "no": "Non-responder",
                }
            )
        )

        response_display["response"] = pd.Categorical(
            response_display["response"],
            categories=[
                "Responder",
                "Non-responder",
            ],
            ordered=True,
        )

        response_display = response_display.sort_values(
            "response"
        )

        response_chart = px.bar(
            response_display,
            x="response",
            y="subject_count",
            color="response",
            labels={
                "response": "Response",
                "subject_count": "Subjects",
            },
            text_auto=True,
            color_discrete_map={
                "Responder": "#5B8E7D",
                "Non-responder": "#7C8796",
            },
        )

        response_chart.update_layout(
            showlegend=False,
            xaxis_title=None,
        )

        st.plotly_chart(
            response_chart,
            width="stretch",
        )

    with sex_col:
        st.subheader("Subjects by Sex")

        sex_display = sexes.copy()

        sex_display["sex"] = sex_display["sex"].map(
            {
                "F": "Female",
                "M": "Male",
            }
        )

        sex_chart = px.bar(
            sex_display,
            x="sex",
            y="subject_count",
            labels={
                "sex": "Sex",
                "subject_count": "Subjects",
            },
            text_auto=True,
            color_discrete_sequence=["#6FA3D2"],
        )

        sex_chart.update_layout(
            showlegend=False,
            xaxis_title=None,
        )

        st.plotly_chart(
            sex_chart,
            width="stretch",
        )

    with st.expander("View All Baseline Samples"):
        cohort_display = cohort.copy()

        cohort_display["response"] = (
            cohort_display["response"].map(
                {
                    "yes": "Responder",
                    "no": "Non-responder",
                }
            )
        )

        cohort_display["sex"] = (
            cohort_display["sex"].map(
                {
                    "F": "Female",
                    "M": "Male",
                }
            )
        )

        cohort_display.columns = [
            "Sample",
            "Subject",
            "Project",
            "Response",
            "Sex",
        ]

        st.caption(
            f"{len(cohort_display):,} samples meet the "
            "baseline melanoma / miraclib / PBMC criteria."
        )

        st.dataframe(
            cohort_display,
            hide_index=True,
            width="stretch",
        )

    st.divider()

    st.subheader(
        "B-Cell Average for Melanoma Male Responders"
    )

    st.write(
        "Melanoma male responders at time 0 across all "
        "treatments and sample types."
    )

    subset_count_col, subset_average_col = st.columns(2)

    subset_count_col.metric(
        "Matching Samples",
        f"{b_cell_summary['sample_count']:,}",
    )

    average_b_cells = b_cell_summary[
        "average_b_cells"
    ]

    subset_average_col.metric(
        "Average B-Cell Count",
        (
            f"{average_b_cells:.2f}"
            if average_b_cells is not None
            else "N/A"
        ),
    )

    if average_b_cells is None:
        st.info(
            "No samples match this subset."
        )
