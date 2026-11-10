# Agronomic and genetic variability module (v1.2.0)

## Input
One observation per genotype × replication × environment, with one or more numeric traits. Download examples from the new tab. For a single environment, select `Single environment`. For multiple environments, use the multi-environment option.

## Outputs
Descriptive statistics, fixed-effect ANOVA, GCV, PCV, entry-mean broad-sense H², genetic advance, genotype means and genotype-mean correlations. Exports are available for each result.

## Assumptions
The first implementation requires a balanced complete experiment with at least two genotypes and replicates. The variance components are ANOVA method-of-moments estimators, not REML. For single-environment RCBD, the ANOVA model is `trait ~ replication + genotype`; for CRD it is `trait ~ genotype`. For multi-environment RCBD it is `trait ~ environment + environment:replication + genotype + environment:genotype`. For multi-environment CRD it is `trait ~ environment + genotype + environment:genotype`.

Negative estimated genetic components are retained as raw values and truncated at zero for derived coefficients. GCV and PCV use the absolute grand mean to avoid sign changes; near-zero means result in NA. The H² estimator is on an entry-mean basis, whereas PCV is plot-level; they must not be confused. Selection intensity uses the normal truncation selection differential `dnorm(qnorm(1-alpha))/alpha`. Genetic advance uses entry-mean phenotypic variance and entry-mean H². No narrow-sense heritability or BLUP is estimated.

**Not independently benchmarked against ASReml, lme4, sommer, agricolae, or other statistical packages.** Validate against trusted reference analyses before using in a publication.
