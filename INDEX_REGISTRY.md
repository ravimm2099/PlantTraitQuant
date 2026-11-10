# Functional index and analysis inventory (final v1.0.0)

This is a **feature inventory**, not a validated count of unique published functional diversity indices. One statistic computed in multiple modules is counted once conceptually. Some results are exploratory proxies.

## Distinct named functional diversity metrics (5 foundational)
- FDis (functional dispersion)
- RaoQ (Rao quadratic entropy)
- FRic_2D (two-dimensional convex hull functional richness proxy)
- FEve (functional evenness)
- FDiv_2D (two-dimensional functional divergence proxy)

## Additional trait-structure outputs (7, including exploratory proxies)
- Nearest-neighbor functional originality
- Centroid-distance functional specialization
- Normalized RaoQ redundancy proxy
- Threshold vulnerability proxy
- Threshold over-redundancy proxy
- Two-dimensional convex-hull area (not a hypervolume)
- CWM-based functional centroid dissimilarity

## Other related outputs, not counted as separate foundational indices
- SES-FDis, SES-RaoQ, functional beta and trait-tree partition outputs
- Species richness, Shannon, Simpson, evenness, taxonomic beta indices
- Trait PCA, Gower/PCoA, RLQ, fourth-corner, genetic marker summaries
- Core and mini-core representativeness, coverage and retention metrics
- Pooled group comparisons and site-accumulation curves

**Auditable count**: 5 foundational FD metrics, plus 7 additional operational trait-structure measures, of which several are proxies. This is **not** a claim of 12 independently validated, distinct published FD indices. A formal full-app index total requires source-level auditing and benchmark validation.

## Multidimensional FRic added in v4.4.0
- Multidimensional convex-hull FRic using 2 to 5 PCA axes, requiring optional `geometry`. This is an extension of FRic, not a new independent index.

## Missing
- Probabilistic / kernel functional hypervolumes
- Full multivariate FRic and FDiv benchmarks
- Formal redundancy/vulnerability index implementations
- Full null-model family and spatially explicit multiscale analysis

## v1.2.0 agronomic statistics

Added descriptive statistics, ANOVA, GCV, PCV, entry-mean broad-sense H², genetic advance, GA% of mean, genotype means and genotype-mean correlation. These are distinct statistics, not all independent diversity indices. See AGRONOMIC_METHODS.md.
