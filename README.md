# XMR

The XMR package implements the XMR (Cross-Population Mendelian Randomization) approach to estimate the causal effect between an exposure and an outcome using genome-wide summary statistics from multiple populations.

XMR is a probabilistic approach for MR analysis in underrepresented populations which improves the power and robustness of causal inference in the target population by leveraging a large-sample auxiliary population.
Specifically, XMR decomposes the observed SNP-trait effects into true causla effects and confounding factors (e.g. pleiotropy, sample structure) hidden in summary statistics.
By explicitly modelling the genetic correlation between two populations, XMR effectively borrows information from the large-sample group.
XMR further corrects bias introduced by IV selection and LD clumping to reduce false positive rates.


<img width="2314" height="1484" alt="649dc9cb-1" src="https://github.com/user-attachments/assets/ee19cdd5-d6ca-4441-b420-641d7bb052df" />


Overview of the XMR method. XMR estimates the causal effect \(\beta\) between exposure \(X_2\) and outcome \(Y_2\) in a small-sample population by leveraging data on the same exposure \(X_1\) from a large-sample population. The method involves several key elements:
\textbf{(A)} IVs are selected from the large-sample population (\(X_1\)) to improve power compared to the limited IVs available from the small-sample population (\(X_2\)). The distributions of observed \(-\log_{10}(p)\) values for SNP-exposure associations across chromosomes are shown.
\textbf{(B)} The XMR model diagram. The arrowed lines represent effects with directions. The blue dashed line indicates the correlation between \(X_1\) and \(X_2\). 
\textbf{(C)} Selection bias and confounding factors contribute to observed SNP-trait associations.
\textbf{(D)} An illustrative example of causal inference between SHBG (sex hormone-binding globulin) and T2D (type 2 diabetes) in an African population, using conventional two-sample MR methods (left) and XMR (right), respectively. The estimated causal effect is shown as a red line, with the 95\% confidence interval shaded in transparent red. Triangles represent observed SNP effect sizes (\(\hat{\gamma}_{2,j}\) and \(\hat{\Gamma}_{2,j}\)), colored by their posterior probability of IV validity (\(Z_j = 1\) in dark blue; \(Z_j = 0\) in light blue).
Figure 1 created with BioRender.com, with permission.

    
## Installation

```r
#install.packages("devtools")
devtools::install_github("yudan121/XMR")
```

## Usage

We illustrate how to perform cross-population MR analysis using XMR by a real example, i.e. LDL cholesterol (LDLC, exposure) and myocardial infarction (MI, outcome), with EUR as the auxiliary population and EAS as the target population.

The XMR analysis comprises the following steps:

**Step 1:** Prepare data and estimate background parameters (C matrix and Omega matrix via cross-population LD score regression).

**Step 2:** Fit XMR for causal inference.

To have a quick look at XMR, you can skip Step 1 and directly jump to Step 2 using the outputs we have prepared.

```r
library(XMR)

exposure <- "LDLC"
outcome  <- "MI"

# Sample sizes
N1 <- 343621  # EUR (auxiliary population)
N2 <- 72866   # EAS (target population)

# Modified IV selection threshold for correction of selection bias
threshold <- 5e-05
t0 <- abs(qnorm(threshold / 2))
dt <- 0.13 / (sqrt(N2 / N1))
modified_threshold <- 2 * (1 - pnorm(abs(t0 + dt)))

# Load example data
data(C)
data(Omega)
data(clumped_data)

# Fit XMR
XMR_res <- fit_XMR(
  data = clumped_data,
  C = C,
  Omega0 = Omega,
  Threshold = modified_threshold,
  tol1 = 1e-07,
  tol2 = 1e-07,
  min_thres = 1e-2
)
```

### Input data format

The input `data` should be a data.frame containing the following columns:

| Column | Description |
|--------|-------------|
| `b.exp.pop1` | SNP-exposure effect in the auxiliary population (pop1) |
| `b.exp.pop2` | SNP-exposure effect in the target population (pop2) |
| `b.out.pop2` | SNP-outcome effect in the target population (pop2) |
| `se.exp.pop1` | Standard error of `b.exp.pop1` |
| `se.exp.pop2` | Standard error of `b.exp.pop2` |
| `se.out.pop2` | Standard error of `b.out.pop2` |
| `L2.pop1` | LD score in the auxiliary population; if NULL, XMR treats it as all 1 |
| `L2.pop2` | LD score in the target population; if NULL, XMR treats it as all 1 |
| `L12` | Cross-population LD score between pop1 and pop2; if NULL, XMR treats it as all 1 |

### Parameters related to confounding factors

- **C matrix**: 3x3 matrix; Captures the effects of sample structure (population stratification, cryptic relatedness, sample overlap, etc.).
- **Omega matrix**: 3x3 matrix; Variance-covariance matrix of polygenic effects.

Both can be estimated using bi-variate LD score regression.

## Reproducibility



## Reference

[Xinrui Huang, Zitong Chao, Zhiwei Wang, Xianghong Hu, and Can Yang]. XMR: A cross-population Mendelian randomization method for causal inference using genome-wide summary statistics. [2026].


## Contact information

Please feel free to contact [Xinrui Huang](xhuangcn@connect.ust.hk), [Prof. Xianghong Hu](huxh@szu.edu.cn), or [Prof. Can Yang](macyang@ust.hk) if you have any questions.
