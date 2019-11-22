
ag=91-20;%starting age is 20
% short-term 5 years
cpd_5=load('rmse_A5_UK_CPD.mat');
svd_5=load('rmse_A5_UK_SVD.mat');
tucker_5aaa=load('rmse_A5_UK_tuckeraaa.mat');
tucker_5aab=load('rmse_A5_UK_tuckeraab.mat');
tucker_5abc=load('rmse_A5_UK_tuckerabc.mat');
cpd5=cpd_5.RMSEage_TS;
svd5=svd_5.RMSEage_LC;
tucker5aaa=tucker_5aaa.RMSEage_TS;
tucker5aab=tucker_5aab.RMSEage_TS;
tucker5abc=tucker_5abc.RMSEage_TS;

figure(1);
plot(20:(ag-1+20),svd5(:,1),'r--*',20:(ag-1+20),cpd5(:,1),'k',20:(ag-1+20),tucker5aaa(:,1),'b--',20:(ag-1+20),tucker5aab(:,1),'m--',20:(ag-1+20),tucker5abc(:,1),'g--')
xlabel('Ages') 
ylabel('RMSFE')
title('The RMSFEs for Total population, UK (Short-term 5 yesrs, Frame A)','FontSize',12)
legend('SVD','CPD','Tucker(aaa)','Tucker(aab)','Tucker(abc)')

% Mid-term 10 years
cpd_10=load('rmse_A10_UK_CPD.mat');
svd_10=load('rmse_A10_UK_SVD.mat');
tucker_10aaa=load('rmse_A10_UK_tuckeraaa.mat');
tucker_10aab=load('rmse_A10_UK_tuckeraab.mat');
tucker_10abc=load('rmse_A10_UK_tuckerabc.mat');
cpd10=cpd_10.RMSEage_TS;
svd10=svd_10.RMSEage_LC;
tucker10aaa=tucker_10aaa.RMSEage_TS;
tucker10aab=tucker_10aab.RMSEage_TS;
tucker10abc=tucker_10abc.RMSEage_TS;

figure(2);
plot(20:(ag-1+20),svd10(:,1),'r--*',20:(ag-1+20),cpd10(:,1),'k',20:(ag-1+20),tucker10aaa(:,1),'b--',20:(ag-1+20),tucker10aab(:,1),'m--',20:(ag-1+20),tucker10abc(:,1),'g--')
xlabel('Ages') 
ylabel('RMSFE')
title('The RMSFEs for Total population, UK (Mid-term 10 yesrs, Frame A)','FontSize',12)
legend('SVD','CPD','Tucker(aaa)','Tucker(aab)','Tucker(abc)')

% Long-term 20 years
cpd_20=load('rmse_A20_UK_CPD.mat');
svd_20=load('rmse_A20_UK_SVD.mat');
tucker_20aaa=load('rmse_A20_UK_tuckeraaa.mat');
tucker_20aab=load('rmse_A20_UK_tuckeraab.mat');
tucker_20abc=load('rmse_A20_UK_tuckerabc.mat');
cpd20=cpd_20.RMSEage_TS;
svd20=svd_20.RMSEage_LC;
tucker20aaa=tucker_20aaa.RMSEage_TS;
tucker20aab=tucker_20aab.RMSEage_TS;
tucker20abc=tucker_20abc.RMSEage_TS;

figure(3);
plot(20:(ag-1+20),svd20(:,1),'r--*',20:(ag-1+20),cpd20(:,1),'k',20:(ag-1+20),tucker20aaa(:,1),'b--',20:(ag-1+20),tucker20aab(:,1),'m--',20:(ag-1+20),tucker20abc(:,1),'g--')
xlabel('Ages') 
ylabel('RMSFE')
title('The RMSFEs for Total population, UK (Long-term 20 yesrs, Frame A)','FontSize',12)
legend('SVD','CPD','Tucker(aaa)','Tucker(aab)','Tucker(abc)')
