# PART 1
The matrix algebra in the provided script calculates Phylogenetically-Informed Predictions (PIP). This method generates climate estimates for fossil species by combining a standard linear regression of physical leaf traits with an adjustment derived from the phylogenetic relationships between the fossil species and modern (extant) species.

The complete prediction is defined by the following equation:

$$\hat{y} = X\beta + V_{cross}^T V_{inv} e$$

This equation is computed in two distinct parts in the R script: the trait-based component and the phylogenetic adjustment.

### 1. The Trait-Based Component ($X\beta$)

This part calculates the expected climate value based strictly on the measured morphology of the fossil leaves, ignoring evolutionary history.

* **$X$ (`X_sp_mat` or `X_site_mat`)**: The design matrix of fossil traits. It has dimensions $m \times p$, where $p$ is the number of trait variables plus an intercept column and $m$ depends on the variant: $m$ = number of fossil species for PIP sp (and LM sp), $m$ = number of species-by-site rows for PIP+site, and $m$ = number of fossil sites for LM site. Each row contains the trait data for one fossil entity.
* **$\beta$ (`pip$beta_mat`)**: A column vector of dimensions $p \times 1$ containing the regression coefficients. These coefficients were estimated from the modern species data using Phylogenetic Generalized Least Squares (PGLS) with Pagel's $\lambda$ jointly estimated.

The operation $X\beta$ represents standard matrix multiplication. The trait values in each row of $X$ are multiplied by their corresponding coefficients in $\beta$ and summed, resulting in a base predicted value of dimensions $m \times 1$. 

### 2. The Phylogenetic Adjustment ($V_{cross}^T V_{inv} e$)

This part calculates the correction factor applied to the trait-based prediction. It uses the known errors of the modern species to adjust the predictions of the fossil species based on how closely related they are on the evolutionary tree.

* **$e$ (`resid_ord_mat`)**: A column vector of dimensions $n \times 1$, where $n$ is the number of modern species. This vector contains the raw residuals (actual climate minus $X\beta$) for the modern species.
* **$V_{inv}$ (`V_inv_mat`)**: A square matrix of dimensions $n \times n$. It is the inverse of the $\lambda$-transformed phylogenetic covariance matrix `V_lam` for the modern species: off-diagonal entries of the raw VCV are multiplied by Pagel's $\lambda$ while the diagonal is preserved (`V_lam <- V * lambda; diag(V_lam) <- diag(V)`). The raw VCV represents shared evolutionary branch lengths between modern species pairs; the $\lambda$-scaling matches the covariance structure under which $\beta$ was estimated, and the inversion accounts for non-independence among the modern species.
* **$V_{cross}^T$ (`t(V_cross_mat)`)**: $V_{cross}$ is a cross-covariance matrix of dimensions $n \times m$, representing the shared evolutionary branch lengths between the $n$ modern species and the $m$ fossil entries. Because no fossil species appears in the extant training dataset, the modern and fossil tip sets are disjoint with respect to the VCV the model was trained on — every entry of $V_{cross}$ is therefore an off-diagonal cell of the full VCV (a shared branch length between two different taxa, never a self-covariance). This means the entire block is uniformly scaled by $\lambda$ (`pip$lambda_mat`) with no diagonal exception, consistent with the $\lambda$-scaling already implicit in $V_{inv}$. The operation `t()` transposes it to dimensions $m \times n$.

**The Sequential Calculation:**

1.  **$V_{inv} e$**: The inverse covariance matrix is multiplied by the modern residuals. This is the GLS "whitening" step: it solves $V x = e$, decorrelating the residuals so that variance shared via common ancestry among the modern species is not double-counted when projected onto the fossils. The result is an $n \times 1$ vector.
2.  **$V_{cross}^T (V_{inv} e)$**: The transposed cross-covariance matrix is multiplied by the decorrelated modern errors. This projects them onto the fossil entries. If a fossil shares a long evolutionary branch with a specific modern species, a larger proportion of that modern species' (decorrelated) residual is applied as an adjustment to the fossil's prediction. The final result is an $m \times 1$ vector (`phylo_adj_mat`).

### 3. The Final Prediction

In Sections 9 and 10 of the script, the two components are added together:

$$\hat{y} = X\beta + \text{phylo\_adj}$$

* **PIP sp**: Uses one row per fossil species, with that species' grand-mean trait values (averaged across all sites it appears in) for matrix $X$. Per-species predictions are then averaged within each fossil site.
* **PIP+site**: Uses one row per species-by-site combination, with the species' site-specific trait values for matrix $X$. The phylogenetic adjustment is indexed by species, so all rows for the same species share the same adjustment regardless of site. Per-(species, site) predictions are averaged within each fossil site.

The final output of $X\beta + \text{phylo\_adj}$ is an $m \times 1$ vector containing the temperature (`yhat_sp_mat`) or precipitation (`yhat_sp_map`) prediction for each fossil input. Because the precipitation model operates on $\log(\text{MAP})$, the script applies the exponential function `exp()` to the final vector to return the values to a standard linear scale.

# PART 2
This example calculates a Phylogenetically-Informed Prediction (PIP) for the Mean Annual Temperature (MAT) of a single fossil species.

Assume the biological trait measured is the proportion of leaves with toothed margins, designated as Tooth Proportion (TP).

For simplicity, assume Pagel's $\lambda = 1$, so $V_{cross}$ and $V_{inv}$ are unscaled (scaling by $\lambda$ would multiply $V_{cross}$ entry-wise by $\lambda$ and replace $V$ with $V_{lam}$ before inversion).

### 1. Initial Parameters and Data

**Modern Species Data**
There are two extant (modern) species, A and B. Their exact environmental MAT and leaf trait (TP) values are known.
* **Species A:** $\text{TP} = 0.8$, $\text{MAT} = 10^\circ\text{C}$
* **Species B:** $\text{TP} = 0.2$, $\text{MAT} = 20^\circ\text{C}$

**Fossil Species Data**
There is one fossil species, F. Its trait is measured from the fossil record, but its MAT is unknown.
* **Species F:** $\text{TP} = 0.5$, $\text{MAT} = ?$

**The Regression Model ($\beta$)**
Take the PGLS coefficients as given for this illustration: intercept = 25 and slope = -15. (In the actual pipeline, $\beta$ comes from `pglmEstLambda()` fit on the full extant species dataset; here we just adopt fixed values so the rest of the algebra is concrete.)
$$\beta = \begin{bmatrix} 25 \\ -15 \end{bmatrix}$$

### 2. Computing the Trait-Based Prediction ($X\beta$)

The design matrix $X$ for the fossil species contains a column for the intercept (always 1) and a column for the measured trait (TP = 0.5).

$$X = \begin{bmatrix} 1 & 0.5 \end{bmatrix}$$

The expected MAT based strictly on the fossil's physical morphology is calculated via matrix multiplication:

$$X\beta = \begin{bmatrix} 1 & 0.5 \end{bmatrix} \begin{bmatrix} 25 \\ -15 \end{bmatrix} = (1 \times 25) + (0.5 \times -15) = 17.5^\circ\text{C}$$

If evolutionary history is ignored, the predicted MAT for Fossil F is 17.5°C.

### 3. Computing the Phylogenetic Adjustment ($V_{cross}^T V_{inv} e$)

To adjust the 17.5°C prediction, the model calculates the error in the modern species and applies it to the fossil species based on their evolutionary relatedness.

**Step 3a: Calculate Modern Residuals ($e$)**
First, calculate the base trait prediction for the modern species using the same $\beta$ coefficients.
* Predicted $\text{MAT}_A = 25 - 15(0.8) = 13^\circ\text{C}$
* Predicted $\text{MAT}_B = 25 - 15(0.2) = 22^\circ\text{C}$

The residual error ($e$) is the actual MAT minus the predicted MAT.
* $e_A = 10 - 13 = -3$
* $e_B = 20 - 22 = -2$

$$e = \begin{bmatrix} -3 \\ -2 \end{bmatrix}$$

**Step 3b: Define Evolutionary Relationships ($V$ and $V_{cross}$)**
Assume the total height of the evolutionary tree from root to tip is 1.0 unit of time.
* Species A and B share an evolutionary branch for 0.4 units of time before diverging.
* Fossil F is closely related to Species A, sharing 0.7 units of time with A, but only 0.4 units of time with B.

The variance-covariance matrix for the modern species ($V$) and its inverse ($V_{inv}$) are:
$$V = \begin{bmatrix} 1.0 & 0.4 \\ 0.4 & 1.0 \end{bmatrix}$$
$$V_{inv} \approx \begin{bmatrix} 1.19 & -0.48 \\ -0.48 & 1.19 \end{bmatrix}$$

The transposed cross-covariance matrix ($V_{cross}^T$) linking the fossil to the modern species is:
$$V_{cross}^T = \begin{bmatrix} 0.7 & 0.4 \end{bmatrix}$$

**Step 3c: Execute the Adjustment Algebra**
Multiply the inverse covariance matrix by the modern residuals to isolate the independent phylogenetic error.
$$V_{inv} e = \begin{bmatrix} 1.19 & -0.48 \\ -0.48 & 1.19 \end{bmatrix} \begin{bmatrix} -3 \\ -2 \end{bmatrix} = \begin{bmatrix} -2.61 \\ -0.94 \end{bmatrix}$$

Multiply the cross-covariance matrix by the isolated error to project the adjustment onto the fossil species.
$$V_{cross}^T (V_{inv} e) = \begin{bmatrix} 0.7 & 0.4 \end{bmatrix} \begin{bmatrix} -2.61 \\ -0.94 \end{bmatrix}$$
$$= (0.7 \times -2.61) + (0.4 \times -0.94) = -1.827 - 0.376 = -2.203^\circ\text{C}$$

The phylogenetic adjustment for Fossil F is -2.203°C. This negative adjustment occurs because its closest relative, Species A, resides in an environment 3°C colder than its leaf traits alone would predict.

### 4. The Final Prediction

Add the phylogenetic adjustment to the trait-based prediction.

$$\hat{y} = X\beta + V_{cross}^T V_{inv} e$$
$$\hat{y} = 17.5 + (-2.203) = 15.297^\circ\text{C}$$

The final Phylogenetically-Informed Prediction for the fossil species' Mean Annual Temperature is 15.3°C.