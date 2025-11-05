# Shrinkage Mean Estimators

This repository provides **R**, **Python**, and **MATLAB** implementations for four shrinkage mean estimators. These methods are designed for high-dimensional data, such as predicting expected returns or constructing mean-variance portfolios in finance.

These estimators share several desirable properties: they are **distribution-free** (not requiring a normal distribution) and use simple shrinkage intensities. They are shown to be more accurate than benchmarks, especially for data with tail-asymmetries, making them highly effective for real-world financial data.

These implementations are based on the research by **[Vali Asimit](https://sites.google.com/view/valiasimit), [Ziwei Chen](https://sites.google.com/view/ziwei-chen), and [Nathan Lassance](https://sites.google.com/view/nathanlassance)**. For full details on the methodology, please see the paper:
**[Distribution-free shrinkage of high-dimensional mean vector with an application to stock returns](https://www.researchgate.net/publication/391900068_Nonparametric_shrinkage_of_high-dimensional_mean_vector_for_return_prediction_and_portfolio_selection)**.

---

## Implemented Estimators

The following four shrinkage estimators are available as functions in each language module:

* `St_MSh_estimator`
* `D_MSh_estimator`
* `O_LSh_estimator`
* `T_LSh_estimator`

---

## Repository Structure

The code is organized by language into three main folders:

* **Matlab/**: Contains `.m` function files for each individual estimator.
* **Python/**: Contains a single script, `mean_shrinkage_estimators.py`, which includes all four functions.
* **R/**: Contains a single script, `mean_shrinkage_estimators.R`, which includes all four functions.

---

## How to Use

These scripts provide the *functions* for the estimators. You must provide your own data matrix (e.g., `X`) to use them.

### Python

1.  Make sure `mean_shrinkage_estimators.py` is in your project directory.
2.  Load your data (e.g., using `pandas` or `numpy`).
3.  Import the functions and call them.

```python
from mean_shrinkage_estimators import St_MSh_estimator, D_MSh_estimator # etc.
X_data = pd.read_csv("my_data.csv").values
theta_St = St_MSh_estimator(X_data)
print(theta_St)
```
### R
1. Load the functions into your R session using `source()`.
2. Load your data (e.g., using `read.csv()`).
3. Call the functions.

```R
source("R/mean_shrinkage_estimators.R")
X_data <- as.matrix(read.csv("my_data.csv", header = FALSE))
theta_St <- St_MSh_estimator(X_data)
print(theta_St)
```
### Matlab
1. Add the **Matlab/** folder to your MATLAB path.
2. Load your data into the workspace (e.g., using readmatrix).
3. Call the functions directly.

```Matlab
addpath('Matlab');
X_data = readmatrix('my_data.csv');
theta_St = St_MSh_estimator(X_data);
disp(theta_St);
```
   



