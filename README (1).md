
# FDA Breast Cancer Analysis using Functional Data Analysis (R)

## 📌 Project Overview
This project applies advanced **Functional Data Analysis (FDA)** techniques in R to analyze breast cancer incidence trends across different age groups. Developed as part of the **MSc in AI for Business Intelligence** at the **University of Leicester**, this research demonstrates smoothing, GCV optimization, covariance visualization, functional principal component analysis (FPCA), and VARIMAX rotation.

## 🧠 Objective
To uncover meaningful patterns and trends in **breast cancer rates by age group** using FDA techniques. The analysis provides deeper insight into variability, trends, and correlations across age-based data curves.

## 📊 Dataset
- **Dataset**: Built-in `Cancerrate` dataset from the `fda` R package.
- **Focus**: Breast cancer incidence rates across **age groups from 45 to 85+**.

## 🛠️ Tools & Libraries Used
- Language: **R**
- Core Package: `fda`
- Plotting: Base R graphics
- Functional Data Techniques: B-spline basis, GCV, PCA, VARIMAX

## ⚙️ Methodology
1. **Smoothing** using B-spline basis with automatic lambda selection via **Generalized Cross-Validation (GCV)**
2. **Descriptive statistics** including mean and standard deviation functions
3. **Covariance structure** visualized through 3D surface and contour plots
4. **Functional Principal Component Analysis (FPCA)** to reduce dimensionality
5. **VARIMAX Rotation** to improve interpretability of principal components

## 🔍 Key Functional Data Analysis Techniques
- Smoothing and basis function expansion
- Cross-validation for optimal smoothing parameter
- Covariance surface estimation
- Eigen decomposition via FPCA
- Rotated harmonics using VARIMAX

## 📈 Results and Visualizations
- Clear smoothed trends in breast cancer rates by age
- Optimal lambda identified through GCV curve
- Covariance surface indicating inter-age relationships
- Scree plot revealing dominant principal components
- Phase-plane plots showcasing dynamic behavior

## 📁 Repository Structure
```
FDA_Breast_Cancer_Analysis.R       # R Script for full analysis
FDA_Breast_Cancer_Report.pdf       # Final project report with visuals
README.md                          # Project overview and documentation
```

## 🚀 How to Run the Code
1. Install R and the `fda` package.
2. Open `FDA_Breast_Cancer_Analysis.R` in RStudio.
3. Run all code blocks to reproduce the analysis and plots.

```r
install.packages("fda")  # if not already installed
```

## 📚 Skills Demonstrated
- Functional data smoothing and basis construction
- Advanced R programming and FDA pipelines
- Statistical interpretation of multivariate data
- Scientific visualization and reporting

## 🏁 Final Remarks
This project showcases the power of functional data analysis in health-related time-series problems. It's ideal for domains like epidemiology, biostatistics, and healthcare analytics.

## 📎 Resources
- Ramsay & Silverman, *Functional Data Analysis* textbook
- FDA R Package documentation
- MSc AI for Business Intelligence curriculum (University of Leicester)

## 👨‍💻 Author
**Satyawan Singh**  
MSc AI for Business Intelligence  
University of Leicester  
🔗 GitHub: [ss1738](https://github.com/ss1738)

---
