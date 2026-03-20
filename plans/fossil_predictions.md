# Plan: Fossil Climate Prediction via PIP

## Context

The PIP (Phylogenetically-Informed Predictions) model in `02_phy_regression.R` predicts climate variables (MAT, MAP) from leaf morphology while accounting for shared evolutionary history. Unlike a standard regression, PIP predictions for any new specimen require three things simultaneously: (1) the regression coefficients β and Pagel's lambda λ, (2) the full VCV matrix of all training species, and (3) the residuals of the training species. A fossil's prediction is:

```
ŷ_fossil = X_fossil · β  +  V[fossil, extant] · V[extant, extant]⁻¹ · residuals[extant]
```

This means fossils must be placed in the phylogeny to obtain `V[fossil, extant]`. Currently, `02_phy_regression.R` does not save the fitted PGLS model or its components, making fossil prediction impossible without modifications.

Target: predict both **MAT** and **log(MAP)** for fossil specimens.

---

## Changes Required

### 1. Modify `code/02_phy_regression.R` — save PIP model components

After fitting the PGLS models for MAT and MAP, save all components needed for fossil prediction:

```r
pip_components <- list(
  # MAT
  beta_mat       = beta,           # GLS coefficients (k × 1)
  lambda_mat     = lambda,         # Pagel's lambda (scalar)
  V_lam_mat      = V_lam,          # lambda-transformed VCV (n × n)
  resid_mat      = resid_incl_all, # training residuals (n × 1) — computed over full data
  X_mat          = X_all,          # design matrix (n × k)
  formula_mat    = mat_formula,

  # MAP (same structure — add analogous MAP block)
  beta_map       = ...,
  lambda_map     = ...,
  V_lam_map      = ...,
  resid_map      = ...,
  X_map          = ...,
  formula_map    = map_formula,

  # Shared
  dat_imputed    = dat_imputed,    # imputed extant data with row names = species
  impute_model   = impute_model,   # caret bagImpute model for trait imputation
  tree_pruned    = phy             # pruned extant tree (ape phylo object)
)
saveRDS(pip_components, file = "models/pip_components.rds")
```

Note: the full-data residuals (`resid_incl_all`) are `dat_imputed$mat - (X_all %*% beta)` computed on all extant species — distinct from the LOOCV loop which excludes one species at a time. These are the residuals used in fossil prediction.

---

### 2. Create `code/03_fossil_predictions.R` — new script

#### Input format

User provides `data/fossil_traits.csv` with columns:

| Column | Description |
|--------|-------------|
| `fossil_name` | Unique identifier (used as tip label) |
| `genus` | Genus for tree placement (mirrors `00_data_cleaning.R` logic) |
| `family` | Family fallback |
| `order` | Order fallback |
| `age_ma` | Fossil age in million years (used to set branch length) |
| `margin.score`, `ln.leaf.area.mm2`, … | Any subset of the trait columns; NAs are imputed |

#### Pipeline

**Step 1 — Load components**
```r
source("code/Phylogenetically-Informed_Predictions_Source.R")
pip <- readRDS("models/pip_components.rds")
fossils <- read.csv("data/fossil_traits.csv")
```

**Step 2 — Graft fossils onto the pruned extant tree**

Reuse the `bind.tip()` logic from `00_data_cleaning.R`, but set edge length so the tip terminates at the fossil's age rather than the present:

```r
h <- max(nodeHeights(tree))
for each fossil:
  find MRCA of genus tips → fallback to family → fallback to order
  edge_len <- nodeheight(tree, target_node) - (h - age_ma * scale_factor)
  tree <- bind.tip(tree, tip.label = fossil_name, where = target_node, edge.length = edge_len)
```

The age scaling assumes tree branch lengths are in the same time units as `age_ma`. This must be verified against `best_wcvp.tre_dated` units (confirm: Ma or relative).

**Step 3 — Recompute VCV for expanded tree**
```r
phylomat_expanded <- vcv(tree_with_fossils)
```

Extract V_lam for the expanded tree using the pre-fit lambda:
```r
V_lam_exp <- phylomat_expanded * lambda
diag(V_lam_exp) <- diag(V_lam_exp) + (1 - lambda)
```

**Step 4 — Impute missing fossil traits**
```r
fossil_imputed <- predict(pip$impute_model, newdata = fossils[, trait_cols])
rownames(fossil_imputed) <- fossils$fossil_name
```

**Step 5 — Build fossil design matrix**
```r
X_fossil <- model.matrix(formula_mat, fossil_imputed)
X_fossil <- X_fossil[, common_vars, drop = FALSE]  # align with beta column order
```

**Step 6 — Apply PIP formula**

For each target variable (MAT, then log MAP):
```r
idx_extant <- rownames(pip$dat_imputed)   # n extant species
idx_fossil  <- fossils$fossil_name         # m fossils

V_inv  <- solve(V_lam_exp[idx_extant, idx_extant])
V_cov  <- V_lam_exp[idx_extant, idx_fossil, drop = FALSE]   # n × m

phylo_adj <- t(V_cov) %*% V_inv %*% residuals   # m × 1
ŷ_fossil  <- X_fossil %*% beta + phylo_adj       # m × 1
```

Compute prediction SE (coefficient uncertainty only, from `pip_components$vcv_mat`):
```r
se_yh <- sqrt(diag(X_fossil %*% vcv_mat %*% t(X_fossil)))
```

**Step 7 — Output**
```r
results <- data.frame(
  fossil_name = fossils$fossil_name,
  mat_pred    = ŷ_mat,
  mat_se      = se_mat,
  mat_95lo    = ŷ_mat - 1.96 * se_mat,
  mat_95hi    = ŷ_mat + 1.96 * se_mat,
  map_pred    = exp(ŷ_map),      # back-transform log(MAP)
  map_se      = se_map,
  map_95lo    = exp(ŷ_map - 1.96 * se_map),
  map_95hi    = exp(ŷ_map + 1.96 * se_map)
)
write.csv(results, "tables/fossil_predictions.csv", row.names = FALSE)
```

---

## Files Modified / Created

| File | Action |
|------|--------|
| `code/02_phy_regression.R` | Add MAP PGLS block (mirrors MAT), save `models/pip_components.rds` |
| `code/03_fossil_predictions.R` | New script — full fossil prediction pipeline |
| `data/fossil_traits.csv` | New input file (user-provided) |
| `tables/fossil_predictions.csv` | New output |

---

## Known Issues to Fix During Implementation

- `addFossil()` in `Phylogenetically-Informed_Predictions_Source.R` (line 475) references `V` instead of `Vmat` — it is **not used** in this plan (we use `bind.tip()` + `vcv()` instead), but the bug should be noted.
- `02_phy_regression.R` currently only implements the MAT LOOCV, not MAP. The MAP PGLS block needs to be added before the save step.

---

## Verification

1. Run `02_phy_regression.R` — confirm `models/pip_components.rds` is written.
2. Create a toy `data/fossil_traits.csv` with 2–3 rows using known extant species (with `mat`/`map` values masked), run `03_fossil_predictions.R`, and confirm predictions are close to the PIP LOOCV estimates from script 02.
3. Check that CIs are wider for fossils with more missing traits (imputed) vs. fossils with complete trait data.
