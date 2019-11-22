function [fore_arma, mse]=armafunction(x,qq,arimapdq); %qq is forecasting period
 % arimapdq =arima(),make a ARIMA(p,d,q) model for later use
c=size(x,2); %column
for i=1:c
[fore_arma(:,i),mse(:,i)]=forecast(estimate(arimapdq,x(:,i)),qq,'Y0',x(:,i));
end
