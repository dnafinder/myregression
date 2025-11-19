 [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/myregression)

📐 myregression — Linear Regression Toolkit for MATLAB

myregression is a unified MATLAB toolbox for linear regression, inverse calibration, and regression comparison.
It includes three harmonized functions originally released on MATLAB FileExchange and now updated, cleaned, and consolidated into a single repository.

✨ What's Included
1. myregr — Linear Regression Analysis

A comprehensive least-squares routine providing:

slope & intercept (value, SE, 95% CI)

Pearson & Spearman correlation

regression standard error

variability decomposition

tests on slope/intercept

homoscedasticity test

regression power

Deming regression

diagnostic plots

2. myregrinv — Inverse Regression & Calibration

Designed for analytical and biostatistical applications.
Features:

inverse prediction (x̂ and 95% CI)

limit of detection (LOD)

limit of quantification (LOQ)

calibration quality index

formatted summary table

3. myregrcomp — Compare Two Regressions

Implements Stanton A. Glantz’s method for comparing two independent linear models:

global F-test

slope comparison (parallelism)

intercept comparison

intersection point

visual comparison plot

🚀 Usage

Place the three .m files on your MATLAB path:

myregr.m

myregrinv.m

myregrcomp.m

Then call them directly in MATLAB using your data vectors.

📦 Requirements

MATLAB

No toolboxes required

All functions fully self-contained

📚 Citation

If you use this toolbox in academic work, please cite:

Cardillo G. (2007–2025).
myregression — A MATLAB toolbox for linear, inverse, and comparative regression.
GitHub: https://github.com/dnafinder/myregression

🔑 License

Released under the GNU GPL-3.0 license.
See the LICENSE file for details.

👤 Author

Giuseppe Cardillo
📧 giuseppe.cardillo.75@gmail.com
