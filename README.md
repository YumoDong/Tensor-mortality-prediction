# Tensor-mortality-prediction
Multi-population Mortality Forecasting using  Tensor Decomposition (Matlab Code)

Software Dependency
===================
The plots/tables were produced under the following software
1. MATLAB R2018b
2. Operating system: Windows 10

[Important Notes]:
1. "Tensor Toolbox" is used to do the tensor decomposition. The homepage is http://www.sandia.gov/~tgkolda/TensorToolbox/.
   "DSP System Toolbox" is used to calculate the RMSE/RMSFE. It can be downloaded in MATLAB.
   "Econometrics Toolbox" is used to the fit time series model. It can be downloaded in MATLAB.
2. The MATLAB code for the Lee-Carter model is refered in the Appendix (figure 9) of this paper  https://pdfs.semanticscholar.org/427c/09210feca24b4b2cabc0d351231a0ffaa2a5.pdf.
3. The data used here are central mortality rates (Mx_1x1), which should be downloaded and processed from The Human Mortality Database (https://www.mortality.org/)

Guidlines for Reproduction
=========================
We use the 10 European countries as an example to show how codes work.
The codes for Frame A and Frame B are separate. e.g. "CPD_A_Eur" is the CPD method for 10 European countries in Frame A. "CPD_B_Eur" is the CPD method for 10 European countries in Frame B.

#################################################################

For CPD method, there is one rank to be determined. 
You only need to change "qq=5" at the begining, which is the length of validation set/testing set.
Short-term forecsat: qq=5
Mid-term forecsat:   qq=10
Long-term forecsat:  qq=20
Extra-term forecsat: qq=30
The rank determined in the validation set is recorded in the "rank_validation".
The decomposed vectors(year/age/countries) are recorded in "CP_new". e.g. "CP_new{1}.U{1}" is the year vector.
The country-specific RMSFEs from CPD model are recorded in the "RMSE_countries".

##################################################################

For the Tucker method, there are three ranks to be determined, so we split the whole task into 2 steps.

* First Step:
Open the "Tucker_abc_A_Eur", which is the code to determine the three ranks. In this file, you only need to change the "qq=5".
For Tucker(aaa) model, "a" is recorded in the "Tucker_aaa_a".
For Tucker(aab) model, "a" is recorded in the "Tucker_aab_a", "b" is recorded in the "Tucker_aab_b". 
For Tucker(abc) model, "a" is recorded in the "Tucker_abc_a", "b" is recorded in the "Tucker_abc_b", "c" is recorded in the "Tucker_abc_c".
 
* Second Step:
Open the "Tucker_testing_A_Eur", which  is the code to forecast mortality. 
In this file, you need to change both "rank=[8 14 9]" and "qq=5".
Remember the "rank" is determined in the * First Step.
The decomposed vectors(year/age/countries) are recorded in "CP_new". e.g. "CP_new{1}.U{1}" is the year vector.
The country-specific RMSFEs from Tucker model are recorded in the "RMSE_countries".
