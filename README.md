# Teiko Technical Assessment

Analysis of immune-cell population data from clinical-trial samples.

The project will provide:

- a reproducible SQLite data pipeline,
- relative-frequency analysis of immune-cell populations,
- responder vs. non-responder statistical analysis,
- the requested baseline cohort queries,
- an interactive dashboard for the analytical results.

Setup and reproduction instructions will be added as the implementation is completed.

## Statistical approach

For the responder/non-responder comparison, each subject contributes one value per immune-cell population: the mean relative frequency across that subject's available longitudinal samples. This avoids treating repeated measurements from the same subject as independent observations.

Each population is compared between responders and non-responders using a two-sided Mann-Whitney U test. P-values across the five population comparisons are adjusted using the Benjamini-Hochberg procedure with an FDR threshold of 0.05. Rank-biserial correlation is reported as an effect-size measure.

In the supplied dataset, no population remains statistically significant after FDR correction. CD4 T cells show the strongest directional difference, with higher relative frequency among responders, but do not meet the adjusted significance threshold.