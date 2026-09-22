# Immune Cell Population Analysis

Reproducible SQLite pipeline, statistical analysis, and Streamlit dashboard for immune-cell population data from clinical-trial samples.

The project loads `cell-count.csv` into SQLite, computes relative cell-population frequencies, compares melanoma responder and non-responder cohorts, and presents the requested Part 2, Part 3, and Part 4 results.

## Dashboard

Public dashboard: https://teiko-immune-cell-dashboard.streamlit.app

Run locally or in GitHub Codespaces:

```bash
make setup
make dashboard
```

The dashboard automatically builds `cell-count.db` from `cell-count.csv` if the database is missing.

## Reproduce the Pipeline

```bash
make setup
make pipeline
```

`make pipeline` rebuilds the database and prints the analytical outputs.

## Data Model

The database has three normalized tables:

- `subjects`: subject-level metadata
- `samples`: sample type and treatment-time metadata
- `cell_counts`: one row per sample and measured population

Measured populations:

- `b_cell`
- `cd8_t_cell`
- `cd4_t_cell`
- `nk_cell`
- `monocyte`

`cell-count.db` is generated and intentionally ignored by git.

## Part 2

Relative frequency is calculated per sample as:

```text
population count / total measured count * 100
```

Output columns:

```text
sample,total_count,population,count,percentage
```

The supplied dataset produces 52,500 population-level relative-frequency rows.

## Part 3

The responder comparison uses melanoma patients receiving `miraclib` with PBMC samples and known response status.

Each subject contributes one value per population: the mean relative frequency across days 0, 7, and 14. This avoids treating repeated samples from the same subject as independent observations.

For each population, responders and non-responders are compared using a two-sided Mann-Whitney U test. The five p-values are adjusted with Benjamini-Hochberg FDR correction at 0.05. Rank-biserial correlation is reported as an effect-size measure.

No cell population remains statistically significant after Benjamini-Hochberg correction at FDR 0.05. CD4 T cells have the smallest adjusted p-value (raw p = 0.0124; adjusted p = 0.0621) and a higher median relative frequency among responders. Their unadjusted result meets the 0.05 threshold, but the adjusted result does not.

## Part 4

For melanoma patients receiving `miraclib` with PBMC samples at time 0:

- baseline samples: 656
- samples by project: `prj1=384`, `prj3=272`
- subjects by response: `no=325`, `yes=331`
- subjects by sex: `F=312`, `M=344`

For melanoma male responders at time 0 across all sample and treatment types, the average B-cell count is:

```text
10206.15
```

This final subset intentionally does not inherit the earlier Part 4 `miraclib` or PBMC filters.

## Tests

Run:

```bash
python -m unittest discover -s tests -v
```

Current suite: 27 tests.

## Files

- `cell-count.csv`: source data
- `schema.sql`: SQLite schema
- `data_validation.py`: CSV validation
- `load_data.py`: database build logic
- `analysis.py`: analytical queries and statistics
- `pipeline.py`: command-line output
- `dashboard.py`: Streamlit dashboard
- `tests/`: unit tests
- `Makefile`: setup, pipeline, dashboard commands
- `requirements.txt`: dependencies
