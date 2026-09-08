# Tissue pO₂ and neuronal death in rTg4510 and wild-type mice

Analysis code accompanying **[manuscript title]**, [authors], [journal/preprint], [year].

This repository reproduces the statistical analysis and figure panel testing whether
neurons that disappear over a 4-week in vivo imaging window are located in relatively
hypoxic tissue, and whether that relationship differs between rTg4510 and wild-type mice.

## Summary of the analysis

Tissue pO₂ was measured at 2,282 neurons across 7 animals (4 rTg4510, 3 WT), of which
24 disappeared over the imaging period. Because baseline oxygenation varies several-fold
between animals (per-animal means 16.0–46.7 mmHg), pO₂ is expressed as each neuron's
deviation from its own animal's mean (ΔpO₂). The animal is the experimental unit
throughout, and all summary statistics are animal-weighted.

Four estimators are fitted, and only effects consistent across all four are interpreted:

| Model | Package | Role |
|---|---|---|
| Mixed-effects logistic regression, random intercept for animal | `lme4` | Primary |
| Firth penalized logistic regression, animal as fixed effect | `brglm2` | Rare-event / conditional alternative |
| GEE, exchangeable correlation, animal as cluster | `geepack` | Directional cross-check only (7 clusters) |
| Within-animal permutation test, 10,000 permutations | base R | Distribution-free |

## Repository contents

```
.
├── pO2_dying_neurons.Rmd    # full analysis; knit to HTML
├── companion_panel.R        # figure panel (within-animal comparison)
├── data/
│   ├── Tg_results.xlsx      # rTg4510 neurons
│   └── WT_results.xlsx      # wild-type neurons
└── figures/                 # output (created by companion_panel.R)
```

### Data dictionary

Both workbooks share one sheet with the following columns:

| Column | Description |
|---|---|
| `CaseID` | Animal identifier |
| `NeuronID` | Neuron identifier, **unique within an animal only** |
| `pO2` | Tissue pO₂ at that neuron (mmHg) |
| `z` | Within-animal z-score of `pO2` |
| `isDying` | `TRUE` if the neuron disappeared over the 4-week window |
| `pO2c` | Within-animal centred pO₂ (`pO2` minus that animal's mean) |

`z` and `pO2c` are recomputed from `pO2` at the top of the analysis and checked against
the supplied columns; the analysis does not depend on them being present.

## Reproducing the analysis

Requires **R ≥ 4.3**. Install dependencies:

```r
install.packages(c(
  "readxl", "dplyr", "tidyr", "ggplot2", "knitr", "rmarkdown",
  "lme4", "broom.mixed", "brglm2", "geepack", "broom", "Hmisc",
  "systemfonts", "ragg"
))
```

Then, from the repository root:

```r
# Full analysis -> pO2_dying_neurons.html
rmarkdown::render("pO2_dying_neurons.Rmd")

# Figure panel -> figures/panel_e.png and figures/panel_e.pdf
source("companion_panel.R")
```

Both scripts assume the repository root as the working directory. Exact package
versions used for the published analysis are recorded in the `sessionInfo()` output
at the end of the rendered HTML.

### Notes on reproducibility

- The permutation test is seeded (`set.seed`), so p-values reproduce exactly.
- `companion_panel.R` renders in Arial where available and falls back to a
  metric-compatible substitute (Liberation Sans, Helvetica, Nimbus Sans) otherwise;
  the layout is identical either way. The font family actually used is printed at
  the top of the run.
- The mixed model's random-intercept variance is estimated at the boundary (zero)
  once genotype is included. This is a boundary estimate, not a convergence failure;
  the Firth model with animal as a fixed effect is reported alongside it as the
  conditional alternative.

## Citation

If you use this code, please cite the manuscript above. [Add the Zenodo DOI badge here
once the first release is archived.]

## License

Code is released under the MIT License (see `LICENSE`). [Confirm the intended license
for the data files with your co-authors before publishing.]
