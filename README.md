# Tensor-mortality-prediction
Multi-population Mortality Forecasting using  Tensor Decomposition (Matlab Code)

Software Dependency
===================
The plots/tables were produced under the following software
1. MATLAB R2018b
2. Operating system: Windows 10

[Important Notes]:
2. The data used here are central mortality rates (Mx_1x1), which should be downloaded and processed from The Human Mortality Database (https://www.mortality.org/)

Guidlines for Reproduction
=========================
1. Before reproduction
   - Setup matlab
     Before using the matlab scripts here, one needs to first install the trmf-exp-0.1 by running the "install.m" script in it.
  
2. For specific tables and figures
    - Table 2
      For the results of M1-M3, run "./Table3.R"; For the result of RMF, run "./trmf-exp-0.1/table_3.m" (need to first install trmf-exp).
    
    - Table 4
      First run "./trmf-exp-0.1/table_5.m" to generate "./trmf-exp-0.1/Rolling Forecast Evaluation.csv" then run "Table4.R”.

    - Table 5 and 6
      First run "Table_5_6.m” to generate "TRMF_future_Female95.csv". Then run "Table_5_6.R” to generate the numbers in the table. Please note that you need to change the value of "age_ax" in "Table_5_6.R”
      to 35,45,55,65 and 75 to get the different values corresponding to the five ages, respectively.
    
    - Figure 6 and 7
      First, need to run "./trmf-exp-0.1/table_3.m" to generate "TRMF_US1933_norm_testmx.csv", then run "Figure_6_7.R".

      
    - Figure 8, 9, and 10
      FIrst, run "normalization_overall_mean_0726.m" to generate "TRMF_future.csv", "TRMFPI_female_future_20.csv" and "TRMFPI_female_future_80.csv" then run "Figure_8_9_10.R”.
