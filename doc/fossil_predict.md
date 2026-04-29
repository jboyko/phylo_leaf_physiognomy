# PART 1
The matrix algebra in the provided script calculates Phylogenetically-Informed Predictions (PIP). This method generates climate estimates for fossil species by combining a standard linear regression of physical leaf traits with an adjustment derived from the phylogenetic relationships between the fossil species and modern (extant) species.

The complete prediction is defined by the following equation:

$$\hat{y} = X\beta + V_{cross}^T V_{inv} e$$

This equation is computed in two distinct parts in the R script: the trait-based component and the phylogenetic adjustment.

### 1. The Trait-Based Component ($X\beta$)

This part calculates the expected climate value based strictly on the measured morphology of the fossil leaves, ignoring evolutionary history.

* **$X$ (`X_sp_mat` or `X_site_mat`)**: The design matrix of fossil traits. It has dimensions $m \times p$, where $m$ is the number of fossil observations and $p$ is the number of trait variables plus an intercept column. Each row contains the trait data for one fossil entity.
* **$\beta$ (`pip$beta_mat`)**: A column vector of dimensions $p \times 1$ containing the regression coefficients. These coefficients were previously estimated from the modern species data using Phylogenetic Generalized Least Squares (PGLS).

The operation $X\beta$ represents standard matrix multiplication. The trait values in each row of $X$ are multiplied by their corresponding coefficients in $\beta$ and summed, resulting in a base predicted value of dimensions $m \times 1$. 

### 2. The Phylogenetic Adjustment ($V_{cross}^T V_{inv} e$)

This part calculates the correction factor applied to the trait-based prediction. It uses the known errors of the modern species to adjust the predictions of the fossil species based on how closely related they are on the evolutionary tree.

* **$e$ (`resid_ord_mat`)**: A column vector of dimensions $n \times 1$, where $n$ is the number of modern species. This vector contains the residuals (actual climate minus trait-predicted climate) for the modern species.
* **$V_{inv}$ (`V_inv_mat`)**: A square matrix of dimensions $n \times n$. It is the inverse of the variance-covariance matrix of the modern species. The original covariance matrix represents the shared evolutionary branch lengths between modern species pairs. The inverse matrix is used to account for the non-independence of the modern species data.
* **$V_{cross}^T$ (`t(V_cross_mat)`)**: $V_{cross}$ is a cross-covariance matrix of dimensions $n \times m$, representing the shared evolutionary branch lengths between the $n$ modern species and the $m$ fossil species. The operation `t()` transposes it to dimensions $m \times n$. The script scales this matrix by Pagel's $\lambda$ (`pip$lambda_mat`), a parameter that quantifies the strength of the phylogenetic signal.

**The Sequential Calculation:**

1.  **$V_{inv} e$**: The inverse covariance matrix is multiplied by the modern residuals. This isolates the independent phylogenetic error for each modern species, removing the redundant error caused by shared ancestry among the modern species themselves. The result is an $n \times 1$ vector.
2.  **$V_{cross}^T (V_{inv} e)$**: The transposed cross-covariance matrix is multiplied by the isolated modern errors. This projects the modern errors onto the fossil species. If a fossil species shares a long evolutionary branch with a specific modern species, a larger proportion of that modern species' residual error is applied as an adjustment to the fossil's prediction. The final result is an $m \times 1$ vector (`phylo_adj_mat`).

### 3. The Final Prediction

In Sections 9 and 10 of the script, the two components are added together:

$$\hat{y} = X\beta + \text{phylo\_adj}$$

* **PIP sp**: Uses the aggregated species-level trait means for matrix $X$.
* **PIP+site**: Uses the specific site-level trait values for matrix $X$, but utilizes the exact same species-level phylogenetic adjustment vector.

The final output is an $m \times 1$ vector containing the temperature (`yhat_sp_mat`) or precipitation (`yhat_sp_map`) prediction for each fossil input. Because the precipitation model operates on $\log(\text{MAP})$, the script applies the exponential function `exp()` to the final vector to return the values to a standard linear scale.

# PART 2
This example calculates a Phylogenetically-Informed Prediction (PIP) for the Mean Annual Temperature (MAT) of a single fossil species.

Assume the biological trait measured is the proportion of leaves with toothed margins, designated as Tooth Proportion (TP). 

### 1. Initial Parameters and Data

**Modern Species Data**
There are two extant (modern) species, A and B. Their exact environmental MAT and leaf trait (TP) values are known.
* **Species A:** $\text{TP} = 0.8$, $\text{MAT} = 10^\circ\text{C}$
* **Species B:** $\text{TP} = 0.2$, $\text{MAT} = 20^\circ\text{C}$

**Fossil Species Data**
There is one fossil species, F. Its trait is measured from the fossil record, but its MAT is unknown.
* **Species F:** $\text{TP} = 0.5$, $\text{MAT} = ?$

**The Regression Model ($\beta$)**
Prior analysis of all modern plant species yielded a standard linear regression model predicting MAT from TP. The coefficients are an intercept of 25 and a slope of -15.
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