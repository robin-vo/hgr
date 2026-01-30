# vog

An R package implementing methods for hierarchical claims reserving models in non-life insurance, currently including: 

- the Negative Binomial Chain-Ladder (NB-CL) model.

## Overview

The NB-CL model provides a full likelihood framework for claim count triangles, addressing limitations of classical Chain-Ladder methods:

- **Full likelihood**: Enables likelihood-based inference, AIC/BIC comparison, and LR tests
- **Overdispersion**: Models variance through dispersion parameter κ with structural interpretation
- **Bias correction**: REML-like correction for finite-sample bias in κ estimation
- **Prediction intervals**: Parametric bootstrap incorporating both process and parameter uncertainty

## Installation

```r
# Install from source
install.packages("vog", repos = NULL, type = "source")

# Or using devtools
devtools::install_github("username/vog")
```

## Quick Start

```r
library(vog)

# Fit model to triangle
fit <- fit_nbcl(triangle)
print(fit)

# Get reserve estimates
reserves <- reserve_nbcl(fit)

# Bootstrap prediction intervals (with REML correction)
boot <- bootstrap_nbcl(fit, B = 5000, correct_kappa = TRUE)
predict_interval(boot, level = 0.95)

# Diagnostics
plot_diagnostics(fit)
```

## Main Functions

| Function | Description |
|----------|-------------|
| `fit_nbcl()` | Fit NB-CL model to run-off triangle |
| `kappa_corrected()` | Extract bias-corrected dispersion parameter |
| `reserve_nbcl()` | Calculate reserve point estimates |
| `bootstrap_nbcl()` | Parametric bootstrap for prediction intervals |
| `predict_interval()` | Extract prediction intervals from bootstrap |
| `profile_kappa()` | Profile likelihood for κ with confidence interval |
| `plot_diagnostics()` | Diagnostic plots (residuals, profile likelihood) |
| `compare_models()` | Compare Poisson, ODP, and NB-CL models |
| `lr_test_overdispersion()` | LR test for overdispersion vs Poisson |

## References

Van Oirbeek, R. (2026). The Negative Binomial Chain-Ladder: A Full Likelihood Model for Claim Count Reserving. Work in progress.

## License

GPL (>= 3)
