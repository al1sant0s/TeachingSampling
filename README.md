# TeachingSampling <img src="Teachingsamplinglogo.png" align="right" height="139" alt="TeachingSampling logo" />

### An R package that draws complex samples and estimates complex parameters


`TeachingSampling` allows you to select samples from the most common probabilistic sampling designs and estimate complex parameters such as totals, means, ratios, regression coefficients, and quantiles.

The package is based on:

> Gutierrez, H. A. (2009). *Estrategias de muestreo: diseño de encuestas y estimación de parámetros*. Editorial Universidad Santo Tomás.

---

## Installation

### Stable version from CRAN

```r
install.packages("TeachingSampling")
```

### Development version from GitHub

```r
install.packages("devtools")
devtools::install_github("psirusteam/TeachingSampling")
```

---

## Functions

### Sampling designs

| Function | Description |
|---|---|
| `S.SI()` | Simple random sampling without replacement |
| `S.SY()` | Systematic sampling |
| `S.BE()` | Bernoulli sampling |
| `S.PO()` | Poisson sampling |
| `S.WR()` | Simple random sampling with replacement |
| `S.PPS()` | PPS sampling with replacement |
| `S.piPS()` | PPS sampling without replacement |
| `S.STSI()` | Stratified simple random sampling |
| `S.STPPS()` | Stratified PPS sampling with replacement |
| `S.STpiPS()` | Stratified PPS sampling without replacement |

### Inclusion probabilities

| Function | Description |
|---|---|
| `PikPPS()` | Inclusion probabilities proportional to size |
| `PikSTPPS()` | Inclusion probabilities for stratified PPS |
| `PikHol()` | Optimal inclusion probabilities (Holmberg) |
| `Pik()` | First-order inclusion probabilities from design |
| `Pikl()` | Second-order inclusion probabilities |

### Estimation

| Function | Description |
|---|---|
| `E.SI()` | Estimation under simple random sampling |
| `E.SY()` | Estimation under systematic sampling |
| `E.BE()` | Estimation under Bernoulli sampling |
| `E.PO()` | Estimation under Poisson sampling |
| `E.WR()` | Estimation under with-replacement sampling |
| `E.PPS()` | Hansen-Hurwitz estimator under PPS-WR |
| `E.piPS()` | HT estimator under piPS sampling |
| `E.STSI()` | Estimation under stratified SI |
| `E.STPPS()` | Estimation under stratified PPS-WR |
| `E.STpiPS()` | Estimation under stratified piPS |
| `E.1SI()` | Estimation under single-stage cluster sampling |
| `E.2SI()` | Estimation under two-stage SI sampling |
| `E.UC()` | Estimation using the Ultimate Cluster method |
| `E.Quantile()` | Weighted quantile estimation |
| `E.Trim()` | Weight trimming and redistribution |

### Regression and calibration

| Function | Description |
|---|---|
| `E.Beta()` | Regression coefficient estimation |
| `GREG.SI()` | Generalised regression estimator |
| `Wk()` | GREG calibration weights |
| `IPFP()` | Iterative proportional fitting (raking) |

### Variance estimation

| Function | Description |
|---|---|
| `VarHT()` | Exact Horvitz-Thompson variance |
| `VarSYGHT()` | HT and Sen-Yates-Grundy variance estimators |
| `HT()` | Horvitz-Thompson estimator |
| `Deltakl()` | Matrix of joint inclusion probability differences |

### Sampling support (small populations)

| Function | Description |
|---|---|
| `Support()` | Sampling support for SI designs |
| `SupportWR()` | Sampling support for WR designs |
| `SupportRS()` | Complete support for all sample sizes |
| `Ik()` | Sample membership indicator matrix |
| `IkWR()` | Frequency indicator matrix for WR sampling |
| `IkRS()` | Indicator matrix for all sample sizes |
| `OrderWR()` | Ordered WR sampling support |
| `nk()` | Frequency matrix for WR sampling |
| `p.WR()` | Sample probabilities under WR sampling |

### Allocation

| Function | Description |
|---|---|
| `kish_allocation()` | Kish compromise allocation for stratified sampling |

### Utilities

| Function | Description |
|---|---|
| `Domains()` | Domain indicator matrix |
| `T.SIC()` | Cluster totals for single-stage sampling |

---

## Usage example

```r
library(TeachingSampling)

data("Lucy")
N <- nrow(Lucy)
n <- 400

# Draw a simple random sample without replacement
sam <- S.SI(N, n)
sam <- sam[sam != 0]

# Estimate population totals
y <- data.frame(Income  = Lucy$Income[sam],
                Taxes    = Lucy$Taxes[sam])

E.SI(N, n, y)
```

---

## Authors

**Hugo Andrés Gutiérrez Rojas** — Package author and maintainer  
Email: hagutierrezro@gmail.com  
GitHub: [@psirusteam](https://github.com/psirusteam)

**Yury Vanessa Ochoa Montes** 
Email: yury.ochoa@urosario.edu.co

---

## Support

- 📖 [Reference manual (CRAN)](https://cran.r-project.org/web/packages/TeachingSampling/TeachingSampling.pdf)
- [CRAN page](https://cran.r-project.org/web/packages/TeachingSampling)
- [Report an issue](https://github.com/psirusteam/TeachingSampling/issues)

Comments, corrections, and suggestions are always welcome.
