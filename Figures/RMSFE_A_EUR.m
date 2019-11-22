ag=91;
% short-term 5years
cpd_5=load('rmse_A5_EUR_CPD.mat');
svd_5=load('rmse_A5_EUR_SVD.mat');
tucker_5222=load('rmse_A5_EUR_tucker222.mat');
tucker_5aaa=load('rmse_A5_EUR_tuckeraaa.mat');
tucker_5aab=load('rmse_A5_EUR_tuckeraab.mat');
tucker_5abc=load('rmse_A5_EUR_tuckerabc.mat');
cpd5=cpd_5.RMSEage_TS;
svd5=svd_5.RMSEage_LC;
tucker5222=tucker_5222.RMSEage_TS;
tucker5aaa=tucker_5aaa.RMSEage_TS;
tucker5aab=tucker_5aab.RMSEage_TS;
tucker5abc=tucker_5abc.RMSEage_TS;

figure(1);
plot(0:(ag-1),svd5(:,1),'r--*',0:(ag-1),cpd5(:,1),'k',0:(ag-1),tucker5222(:,1),'c--*',0:(ag-1),tucker5aaa(:,1),'b--',0:(ag-1),tucker5aab(:,1),'m--',0:(ag-1),tucker5abc(:,1),'g--')
xlabel('Ages') 
ylabel('RMSFE')
title('The RMSFEs for Total population, Europe (Short-term 5 yesrs, Frame A)','FontSize',12)
legend('SVD','CPD','Tucker(222)','Tucker(aaa)','Tucker(aab)','Tucker(abc)')

% mid-term 10 years
cpd_10=load('rmse_A10_EUR_CPD.mat');
svd_10=load('rmse_A10_EUR_SVD.mat');
tucker_10222=load('rmse_A10_EUR_tucker222.mat');
tucker_10aaa=load('rmse_A10_EUR_tuckeraaa.mat');
tucker_10aab=load('rmse_A10_EUR_tuckeraab.mat');
tucker_10abc=load('rmse_A10_EUR_tuckerabc.mat');

cpd10=cpd_10.RMSEage_TS;
svd10=svd_10.RMSEage_LC;
tucker10222=tucker_10222.RMSEage_TS;
tucker10aaa=tucker_10aaa.RMSEage_TS;
tucker10aab=tucker_10aab.RMSEage_TS;
tucker10abc=tucker_10abc.RMSEage_TS;

figure(2);
plot(0:(ag-1),svd10(:,1),'r--*',0:(ag-1),cpd10(:,1),'k',0:(ag-1),tucker10222(:,1),'c--*',0:(ag-1),tucker10aaa(:,1),'b--',0:(ag-1),tucker10aab(:,1),'m--',0:(ag-1),tucker10abc(:,1),'g--')
xlabel('Ages') 
ylabel('RMSFE')
title('The RMSFEs for Total population, Europe (Mid-term 10 yesrs, Frame A)','FontSize',12)
legend('SVD','CPD','Tucker(222)','Tucker(aaa)','Tucker(aab)','Tucker(abc)')

% long-term 20 years
cpd_20=load('rmse_A20_EUR_CPD.mat');
svd_20=load('rmse_A20_EUR_SVD.mat');
tucker_20222=load('rmse_A20_EUR_tucker222.mat');
tucker_20aaa=load('rmse_A20_EUR_tuckeraaa.mat');
tucker_20aab=load('rmse_A20_EUR_tuckeraab.mat');
tucker_20abc=load('rmse_A20_EUR_tuckerabc.mat');

cpd20=cpd_20.RMSEage_TS;
svd20=svd_20.RMSEage_LC;
tucker20222=tucker_20222.RMSEage_TS;
tucker20aaa=tucker_20aaa.RMSEage_TS;
tucker20aab=tucker_20aab.RMSEage_TS;
tucker20abc=tucker_20abc.RMSEage_TS;

figure(3);
plot(0:(ag-1),svd20(:,1),'r--*',0:(ag-1),cpd20(:,1),'k',0:(ag-1),tucker20222(:,1),'c--*',0:(ag-1),tucker20aaa(:,1),'b--',0:(ag-1),tucker20aab(:,1),'m--',0:(ag-1),tucker20abc(:,1),'g--')
xlabel('Ages') 
ylabel('RMSFE')
title('The RMSFEs for Total population, Europe (Long-term 20 yesrs, Frame A)','FontSize',12)
legend('SVD','CPD','Tucker(222)','Tucker(aaa)','Tucker(aab)','Tucker(abc)')

% Extra-term 30 years
cpd_30=load('rmse_A20_EUR_CPD.mat');
svd_30=load('rmse_A20_EUR_SVD.mat');
tucker_30222=load('rmse_A30_EUR_tucker222.mat');
tucker_30aaa=load('rmse_A30_EUR_tuckeraaa.mat');
tucker_30aab=load('rmse_A30_EUR_tuckeraab.mat');
tucker_30abc=load('rmse_A30_EUR_tuckerabc.mat');

cpd30=cpd_30.RMSEage_TS;
svd30=svd_30.RMSEage_LC;
tucker30222=tucker_30222.RMSEage_TS;
tucker30aaa=tucker_30aaa.RMSEage_TS;
tucker30aab=tucker_30aab.RMSEage_TS;
tucker30abc=tucker_30abc.RMSEage_TS;

figure(4);
plot(0:(ag-1),svd30(:,1),'r--*',0:(ag-1),cpd30(:,1),'k',0:(ag-1),tucker30222(:,1),'c--*',0:(ag-1),tucker30aaa(:,1),'b--',0:(ag-1),tucker30aab(:,1),'m--',0:(ag-1),tucker30abc(:,1),'g--')
xlabel('Ages') 
ylabel('RMSFE')
title('The RMSFEs for Total population, Europe (Extra-term 30 yesrs, Frame A)','FontSize',12)
legend('SVD','CPD','Tucker(222)','Tucker(aaa)','Tucker(aab)','Tucker(abc)')
