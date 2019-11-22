%STEP 1: Import Data
%Entering some parameters
rank=[3 7 1];
qq=20; %length of valiadation set

T_begin=1922; %the begining year
T_end=2016; % the ending year
tt=T_end-T_begin+1-2*qq; %length of test set
ag=91-20; %length of ages
cc=2; %number of genders
rms3d= dsp.RMS; %RMS fuction
rms3d.Dimension='ALL'; %calculate the whole tensor's RMS
arimapdq=arima(0,1,0); %make a ARIMA(0,1,0) model for later use
trys=300; %This this is used to determine best fit of final model in test set
male_d=4;
female_d=3;
total_d=5;  %3,4,5 represent different death rate 
xx=15;

% UK data
%In this file,CHE represents UK female, DNK represents UK male

CHE1=load('GBR.mat'); %UK female
yyy=CHE1.GBRNP;
xxx=yyy(:,female_d);  
che1=table2array(xxx); %table to array
zzz=zeros(ag,tt); 
for n=1:tt   
    zzz(1:ag,n)=che1(1+0+20+(n-1)*111:0+111*n-20,1);
end
lnCHE_t=log(zzz);  %ln(m)
MeanCHE1=mean(lnCHE_t,2); % the mean of ln(mx) at a given age across test set.
lnCHE_ct=zeros(ag,tt);
for m=1:tt
    lnCHE_ct(:,m)=lnCHE_t(:,m)-MeanCHE1;
end
lnCHE_ct; %Centering ln(mx) 
% Validation Set 
qqq =zeros(ag,qq); %Validation Matrix
for n=1:qq
    qqq(1:ag,n)=che1(1+20+111*tt+(n-1)*111:111*tt+111*n-20,1);
end
lnCHE_v=log(qqq);  %the ln(m) 
MeanCHE2=mean(lnCHE_v,2); % the mean of ln(mx) at a given age across test set.
lnCHE_cv=zeros(ag,qq);
for m=1:qq
    lnCHE_cv(:,m)=lnCHE_v(:,m)-MeanCHE2;
end
lnCHE_cv; %Centering ln(mx) 
%Remove variables
clear CHE1
clear xxx 
clear zzz
clear qqq
clear n
clear i
clear yyy

DNK1=load('GBR.mat'); %UK male
yyy=DNK1.GBRNP;
xxx=yyy(:,male_d);   
dnk1=table2array(xxx); %table to array
zzz=zeros(ag,tt);
for n=1:tt   
    zzz(1:ag,n)=dnk1(1+0+20+(n-1)*111:0+111*n-20,1);
end
lnDNK_t=log(zzz);  %ln(m)
MeanDNK1=mean(lnDNK_t,2); % the mean of ln(mx) at a given age across test set.
lnDNK_ct=zeros(ag,tt);
for m=1:tt
    lnDNK_ct(:,m)=lnDNK_t(:,m)-MeanDNK1;
end
lnDNK_ct; %Centering ln(mx) 
% Validation Set 
qqq =zeros(ag,qq); % Validation Matrix
for n=1:qq
    qqq(1:ag,n)=dnk1(1+20+111*tt+(n-1)*111:111*tt+111*n-20,1);
end
lnDNK_v=log(qqq);  %the ln(m) 
MeanDNK2=mean(lnDNK_v,2); % the mean of ln(mx) at a given age across test set.
lnDNK_cv=zeros(ag,qq);
for m=1:qq
    lnDNK_cv(:,m)=lnDNK_v(:,m)-MeanDNK2;
end
lnDNK_cv; %Centering ln(mx) 
%Remove variables
clear DNK1
clear xxx 
clear zzz
clear qqq
clear n
clear i
clear yyy

% new training set_ Target tensor.
lnCHE_new=[lnCHE_t,lnCHE_v]; 
lnDNK_new=[lnDNK_t,lnDNK_v]; 

% get means of the training+val
MeanCHE=mean(lnCHE_new,2);% the Mean of the year vectors
MeanDNK=mean(lnDNK_new,2);% the Mean of the year vectors

%
lnCHE_newc=zeros(ag,tt+qq); %Centering log(mx) 
lnDNK_newc=zeros(ag,tt+qq); %Centering log(mx) 

for m=1:(tt+qq);
    
    lnCHE_newc(:,m)=lnCHE_new(:,m)-MeanCHE;
    lnDNK_newc(:,m)=lnDNK_new(:,m)-MeanDNK;
     
end
clear m
%
%We use after centering ln(mx) to consruct new Tensor.
TT_new=zeros(ag,tt+qq,cc); %ag*(tt+qq)*cc tensor
TT_new(:,:,1)=lnCHE_newc;
TT_new(:,:,2)=lnDNK_newc;
TT_new=tensor(TT_new); %This is after centering tensor

% The real log(mx) before centering
TT_real=zeros(ag,tt+qq,cc);

TT_real(:,:,1)=lnCHE_new; %This is ag*tt+qq*cc before centering data
TT_real(:,:,2)=lnDNK_new;

%Mean of tt+qq
Mean_tq=[MeanCHE,MeanDNK]

%choose best fit [a b c] decomposition(training+validation)                   
 CP_newf=cell(trys,1);
parfor rr=1:trys;
    rng(rr) % random seeds
    CP_newf(rr)= {tucker_als(TT_new,rank)};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CP_ranknew=zeros(ag,tt+qq,cc);
RMSE_ranknew=zeros(trys,1);
for rr=1:trys;
    CP_ranknewf=double(CP_newf{rr,1});% imput decomposition after centering
   for ci=1:cc
    for m=1:tt+qq
    CP_ranknew(:,m,ci)=CP_ranknewf(:,m,ci)+Mean_tq(:,ci);
    end
    end
    %we reverse-centering the tesnor,which is log(mx) now.
    CP_ranknewd=CP_ranknew-TT_real;
    RMSE_ranknew(rr,1)=rms3d(CP_ranknewd); %300 times of RMSEs
end
clear m;
[Minnew,Seqnew]=min(RMSE_ranknew);% the best one 
CP_new=CP_newf(Seqnew); % Now we get best rank decomposition.

cpnew_age=CP_new{1}.U{1};  % age vector of tensor
cpnew_country=CP_new{1}.U{3}; %gender vector of tensor
cpnew_lambda=CP_new{1}.core; %the lambda
cpnew_year=arimafunction(CP_new{1}.U{2},qq,arimapdq); %the forecast year vectors,random walk

FC_cpnew=ttm(cpnew_lambda,{cpnew_age,cpnew_year,cpnew_country});%the forecast tensor in ktensor format
FCT_cpnew=full(FC_cpnew); %Transfer ktensor to tensor
FCD_cpnew=zeros(ag,qq,cc); %transfer tensor to double format
for i=1:cc
FCD_cpnew(:,:,i)=FCT_cpnew(:,:,i);
end         %remember this is centering log(mx).
%
FCD_cpnewb=FCD_cpnew;
parfor ci=1:cc;
    FCD_cpnewb(:,:,ci)=FCD_cpnew(:,:,ci)+Mean_tq(:,ci);
end%we reverse-centering the tesnor,which is log(mx) now.
clear m
FCD_cpnewb; 
%this is our forecast tensor which will be compared with the Lee-Carter model.

%Compare with the Lee-Carter model by RMSE
%Import real observations 
 % Forecasting Matrix
CHE_new =zeros(ag,qq);
DNK_new =zeros(ag,qq);
parfor n=1:qq
    
   CHE_new(:,n)=che1(1+20+111*(tt+qq)+(n-1)*111:111*(tt+qq)+111*n-20,1);
    DNK_new(:,n)=dnk1(1+20+111*(tt+qq)+(n-1)*111:111*(tt+qq)+111*n-20,1);
    
end
%
%impute real log(observations) 
Realob_new=zeros(ag,qq,cc);

Realob_new(:,:,1)=log(CHE_new); %they are log(mx) before the centering process
Realob_new(:,:,2)=log(DNK_new);

%RMSE of forecast and gender specific RMSE.
RMSE_cpnew= Realob_new-FCD_cpnewb; %get the difference between the real and the forecst
RMSE_new=rms3d(RMSE_cpnew); %the RMSE of tensor add the mean back.
% RMSE of each genders
RMSE_countries=zeros(cc,1);
parfor ci=1:cc;
    gap=Realob_new(:,:,ci)-FCD_cpnewb(:,:,ci);
    RMSE_countries(ci)=rms3d(gap);
end
                
