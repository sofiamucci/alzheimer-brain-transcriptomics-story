# 13_covariate_confounding_check.R
# Checks whether donor covariates (RIN, PMI) from 01b_add_donor_covariates.R
# are confounded with condition (Young/Old/AD) — i.e., whether they differ
# systematically between groups, which would bias DESeq2 results run with
# the published design (~ condition, no covariates).
#
# Result: RIN differs significantly by condition (AD lowest); PMI does not.
# Documented as a known limitation of the current design in the README —
# the published pipeline was NOT re-run with RIN as a covariate, to keep
# the main analysis reproducible as originally published.

library(ggplot2)

meta <- read.csv("data/sample_metadata_full.csv")
meta$condition <- factor(meta$condition, levels = c("Young", "Old", "AD"))

ggplot(meta, aes(x = condition, y = rin, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 2, alpha = 0.7) +
  theme_bw() + theme(legend.position = "none") +
  labs(title = "RIN by condition", y = "RIN", x = NULL)

ggplot(meta, aes(x = condition, y = pmi_hr, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 2, alpha = 0.7) +
  theme_bw() + theme(legend.position = "none") +
  labs(title = "PMI by condition", y = "PMI (hr)", x = NULL)

kruskal.test(rin ~ condition, data = meta)     # p = 0.003 — significant
kruskal.test(pmi_hr ~ condition, data = meta)  # p = 0.78 — not significant