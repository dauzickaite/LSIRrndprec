# LSIRrndprec

Code for reproducing results in E. Carson and I. Dauzickaite, `Mixed precision sketching for least-squares problems and its application in GMRES-based iterative refinement', 2024, https://doi.org/10.48550/arXiv.2410.06319. 

MATLAB R2023b is used.

Run plot1_norms.m for the plots in the first figure, plot2_cond.m for the plots in the second figure. To obtain results with the working precision set to single and double run, respectively, main_u_single.m and main_u_double.com.

Note: requires chop (https://github.com/higham/chop) and advanpix toolbox (https://www.advanpix.com/).
