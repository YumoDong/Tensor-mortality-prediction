
 %STEP 1: Import Data
% Initial settings

qq=5; %length of valiadation/testing set (year), you can change it to change the forecasting horizon.

T_begin=1950; %the begining year
T_end=2014; % the ending year
tt=T_end-T_begin+1-2*qq; %length of test set
sa=0; % starting age
ag=91-sa; %length of ages
cc=10; %number of countries/genders
rms3d= dsp.RMS; %RMS fuction
rms3d.Dimension='ALL'; %calculate the whole tensor's RMS
arimapdq=arima(0,1,0); %make a random model for later use
seeds=300; %This this is used to determine best fit of final model in validation set
trys=300; %This this is used to determine best fit of final model in test set
male_d=4;
female_d=3;
total_d=5;  %3,4,5 represent different death rate total_d=5;  %3 is female, 4 is male, 5 is total 

% UK data
GBR1=load('GBR.mat');
UK=GBR1.GBRNP;
uk=UK(:,total_d);   %UK'total death rate in column 5
uk1=table2array(uk); %table to array
UK_t=zeros(ag,tt); %(91 ages * tt years)
for n=1:tt   
    UK_t(1:ag,n)=uk1(1+3108+(n-1)*111:3108+111*n-20,1);
end
lnUK_t=log(UK_t);  %ln(m) 1922-1922+tt-1 year with 0-90 age
MeanUK1=mean(lnUK_t,2); % the mean of log(mx) at a given age across test set.
lnUK_ct=zeros(ag,tt);
for m=1:tt
    lnUK_ct(:,m)=lnUK_t(:,m)-MeanUK1;
end
lnUK_ct; %Centering ln(mx) of UK from 1922-1922+tt-1  
% Validation Set (1922+tt to 1922+tt+qq-1)
UK_v =zeros(ag,qq); %(91 ages * qq years) Validation Matrix
for n=1:qq
    UK_v(1:ag,n)=uk1(1+3108+111*tt+(n-1)*111:111*tt+3108+111*n-20,1);
end
lnUK_v=log(UK_v);  %the ln(m) of 1922+tt to 1922+tt+qq-1 *age 0-90
MeanUK2=mean(lnUK_v,2); % the mean of log(mx) at a given age across test set.
lnUK_cv=zeros(ag,qq);
for m=1:qq
    lnUK_cv(:,m)=lnUK_v(:,m)-MeanUK2;
end
lnUK_cv; %Centering ln(mx) of UK from 1990-1999 
%Remove variables
clear GBR1
clear uk
clear num
clear n
clear i
clear UK_t
clear UK_v

CHE1=load('CHE.mat'); %Switherland
yyy=CHE1.CHE;
xxx=yyy(:,total_d);   %total death rate in column 5
che1=table2array(xxx); %table to array
zzz=zeros(ag,tt); %(91 ages * tt years)
for n=1:tt   
    zzz(1:ag,n)=che1(1+8214+(n-1)*111:8214+111*n-20,1);
end
lnCHE_t=log(zzz);  %ln(m) 1922-1922+tt-1 year with 0-90 age
MeanCHE1=mean(lnCHE_t,2); % the mean of log(mx) at a given age across test set.
lnCHE_ct=zeros(ag,tt);
for m=1:tt
    lnCHE_ct(:,m)=lnCHE_t(:,m)-MeanCHE1;
end
lnCHE_ct; %Centering ln(mx) from 1922-1922+tt-1  
% Validation Set (1922+tt to 1922+tt+qq-1)
qqq =zeros(ag,qq); %(91 ages * qq years) Validation Matrix
for n=1:qq
    qqq(1:ag,n)=che1(1+8214+111*tt+(n-1)*111:8214+111*tt+111*n-20,1);
end
lnCHE_v=log(qqq);  %the ln(m) of 1922+tt to 1922+tt+qq-1 *age 0-90
MeanCHE2=mean(lnCHE_v,2); % the mean of log(mx) at a given age across test set.
lnCHE_cv=zeros(ag,qq);
for m=1:qq
    lnCHE_cv(:,m)=lnCHE_v(:,m)-MeanCHE2;
end
lnCHE_cv; %Centering ln(mx) of 1922+tt to 1922+tt+qq-1
%Remove variables
clear CHE1
clear xxx 
clear zzz
clear qqq
clear n
clear i
clear yyy

DNK1=load('DNK.mat'); %Denmark
yyy=DNK1.DNK1;
xxx=yyy(:,total_d);   %total death rate in column 5
dnk1=table2array(xxx); %table to array
zzz=zeros(ag,tt); %(91 ages * tt years)
for n=1:tt   
    zzz(1:ag,n)=dnk1(1+12765+(n-1)*111:12765+111*n-20,1);
end
lnDNK_t=log(zzz);  %ln(m) 1922-1922+tt-1 year with 0-90 age
MeanDNK1=mean(lnDNK_t,2); % the mean of log(mx) at a given age across test set.
lnDNK_ct=zeros(ag,tt);
for m=1:tt
    lnDNK_ct(:,m)=lnDNK_t(:,m)-MeanDNK1;
end
lnDNK_ct; %Centering ln(mx) from 1922-1922+tt-1  
% Validation Set (1922+tt to 1922+tt+qq-1)
qqq =zeros(ag,qq); %(91 ages * qq years) Validation Matrix
for n=1:qq
    qqq(1:ag,n)=dnk1(1+12765+111*tt+(n-1)*111:12765+111*tt+111*n-20,1);
end
lnDNK_v=log(qqq);  %the ln(m) of 1922+tt to 1922+tt+qq-1 *age 0-90
MeanDNK2=mean(lnDNK_v,2); % the mean of log(mx) at a given age across test set.
lnDNK_cv=zeros(ag,qq);
for m=1:qq
    lnDNK_cv(:,m)=lnDNK_v(:,m)-MeanDNK2;
end
lnDNK_cv; %Centering ln(mx) of 1922+tt to 1922+tt+qq-1
%Remove variables
clear DNK1
clear xxx 
clear zzz
clear qqq
clear n
clear i
clear yyy

ESP1=load('ESP.mat'); %Spain
yyy=ESP1.ESP;
xxx=yyy(:,total_d);   %total death rate in column 5
esp1=table2array(xxx); %table to array
zzz=zeros(ag,tt); %(91 ages * tt years)
for n=1:tt   
    zzz(1:ag,n)=esp1(1+4662+(n-1)*111:4662+111*n-20,1);
end
lnESP_t=log(zzz);  %ln(m) 1922-1922+tt-1 year with 0-90 age
MeanESP1=mean(lnESP_t,2); % the mean of log(mx) at a given age across test set.
lnESP_ct=zeros(ag,tt);
for m=1:tt
    lnESP_ct(:,m)=lnESP_t(:,m)-MeanESP1;
end
lnESP_ct; %Centering ln(mx) from 1922-1922+tt-1  
% Validation Set (1922+tt to 1922+tt+qq-1)
qqq =zeros(ag,qq); %(91 ages * qq years) Validation Matrix
for n=1:qq
    qqq(1:ag,n)=esp1(1+4662+111*tt+(n-1)*111:4662+111*tt+111*n-20,1);
end
lnESP_v=log(qqq);  %the ln(m) of 1922+tt to 1922+tt+qq-1 *age 0-90
MeanESP2=mean(lnESP_v,2); % the mean of log(mx) at a given age across test set.
lnESP_cv=zeros(ag,qq);
for m=1:qq
    lnESP_cv(:,m)=lnESP_v(:,m)-MeanESP2;
end
lnESP_cv; %Centering ln(mx) of 1922+tt to 1922+tt+qq-1
%Remove variables
clear ESP1
clear xxx 
clear zzz
clear qqq
clear n
clear i
clear yyy

FIN1=load('FIN.mat'); %Finland
yyy=FIN1.FIN;
xxx=yyy(:,total_d);   %total death rate in column 5
fin1=table2array(xxx); %table to array
zzz=zeros(ag,tt); %(91 ages * tt years)
for n=1:tt   
    zzz(1:ag,n)=fin1(1+7992+(n-1)*111:7992+111*n-20,1);
end
lnFIN_t=log(zzz);  %ln(m) 1922-1922+tt-1 year with 0-90 age
MeanFIN1=mean(lnFIN_t,2); % the mean of log(mx) at a given age across test set.
lnFIN_ct=zeros(ag,tt);
for m=1:tt
    lnFIN_ct(:,m)=lnFIN_t(:,m)-MeanFIN1;
end
lnFIN_ct; %Centering ln(mx) from 1922-1922+tt-1  
% Validation Set (1922+tt to 1922+tt+qq-1)
qqq =zeros(ag,qq); %(91 ages * qq years) Validation Matrix
for n=1:qq
    qqq(1:ag,n)=fin1(1+7992+111*tt+(n-1)*111:7992+111*tt+111*n-20,1);
end
lnFIN_v=log(qqq);  %the ln(m) of 1922+tt to 1922+tt+qq-1 *age 0-90
MeanFIN2=mean(lnFIN_v,2); % the mean of log(mx) at a given age across test set.
lnFIN_cv=zeros(ag,qq);
for m=1:qq
    lnFIN_cv(:,m)=lnFIN_v(:,m)-MeanFIN2;
end
lnFIN_cv; %Centering ln(mx) of 1922+tt to 1922+tt+qq-1
%Remove variables
clear FIN1
clear xxx 
clear zzz
clear qqq
clear n
clear i
clear yyy

ITA1=load('ITA.mat');%Italy
yyy=ITA1.ITA;
xxx=yyy(:,total_d);   %total death rate in column 5
ita1=table2array(xxx); %table to array
zzz=zeros(ag,tt); %(91 ages * tt years)
for n=1:tt   
    zzz(1:ag,n)=ita1(1+8658+(n-1)*111:8658+111*n-20,1);
end
lnITA_t=log(zzz);  %ln(m) 1922-1922+tt-1 year with 0-90 age
MeanITA1=mean(lnITA_t,2); % the mean of log(mx) at a given age across test set.
lnITA_ct=zeros(ag,tt);
for m=1:tt
    lnITA_ct(:,m)=lnITA_t(:,m)-MeanITA1;
end
lnITA_ct; %Centering ln(mx) from 1922-1922+tt-1  
% Validation Set (1922+tt to 1922+tt+qq-1)
qqq =zeros(ag,qq); %(91 ages * qq years) Validation Matrix
for n=1:qq
    qqq(1:ag,n)=ita1(1+8658+111*tt+(n-1)*111:8658+111*tt+111*n-20,1);
end
lnITA_v=log(qqq);  %the ln(m) of 1922+tt to 1922+tt+qq-1 *age 0-90
MeanITA2=mean(lnITA_v,2); % the mean of log(mx) at a given age across test set.
lnITA_cv=zeros(ag,qq);
for m=1:qq
    lnITA_cv(:,m)=lnITA_v(:,m)-MeanITA2;
end
lnITA_cv; %Centering ln(mx) of 1922+tt to 1922+tt+qq-1
%Remove variables
clear ITA1
clear xxx 
clear zzz
clear qqq
clear n
clear i
clear yyy

NLD1=load('NLD.mat');%Netherlands
yyy=NLD1.NLD;
xxx=yyy(:,total_d);   %total death rate in column 5
nld1=table2array(xxx); %table to array
zzz=zeros(ag,tt); %(91 ages * tt years)
for n=1:tt   
    zzz(1:ag,n)=nld1(1+11100+(n-1)*111:11100+111*n-20,1);
end
lnNLD_t=log(zzz);  %ln(m) 1922-1922+tt-1 year with 0-90 age
MeanNLD1=mean(lnNLD_t,2); % the mean of log(mx) at a given age across test set.
lnNLD_ct=zeros(ag,tt);
for m=1:tt
    lnNLD_ct(:,m)=lnNLD_t(:,m)-MeanNLD1;
end
lnNLD_ct; %Centering ln(mx) from 1922-1922+tt-1  
% Validation Set (1922+tt to 1922+tt+qq-1)
qqq =zeros(ag,qq); %(91 ages * qq years) Validation Matrix
for n=1:qq
    qqq(1:ag,n)=nld1(1+11100+111*tt+(n-1)*111:11100+111*tt+111*n-20,1);
end
lnNLD_v=log(qqq);  %the ln(m) of 1922+tt to 1922+tt+qq-1 *age 0-90
MeanNLD2=mean(lnNLD_v,2); % the mean of log(mx) at a given age across test set.
lnNLD_cv=zeros(ag,qq);
for m=1:qq
    lnNLD_cv(:,m)=lnNLD_v(:,m)-MeanNLD2;
end
lnNLD_cv; %Centering ln(mx) of 1922+tt to 1922+tt+qq-1
%Remove variables
clear NLD1
clear xxx 
clear zzz
clear qqq
clear n
clear i
clear yyy


NOR1=load('NOR.mat'); %Norway
yyy=NOR1.NOR;
xxx=yyy(:,total_d);   %total death rate in column 5
nor1=table2array(xxx); %table to array
zzz=zeros(ag,tt); %(91 ages * tt years)
for n=1:tt   
    zzz(1:ag,n)=nor1(1+11544+(n-1)*111:11544+111*n-20,1);
end
lnNOR_t=log(zzz);  %ln(m) 1922-1922+tt-1 year with 0-90 age
MeanNOR1=mean(lnNOR_t,2); % the mean of log(mx) at a given age across test set.
lnNOR_ct=zeros(ag,tt);
for m=1:tt
    lnNOR_ct(:,m)=lnNOR_t(:,m)-MeanNOR1;
end
lnNOR_ct; %Centering ln(mx) from 1922-1922+tt-1  
% Validation Set (1922+tt to 1922+tt+qq-1)
qqq =zeros(ag,qq); %(91 ages * qq years) Validation Matrix
for n=1:qq
    qqq(1:ag,n)=nor1(1+11544+111*tt+(n-1)*111:11544+111*tt+111*n-20,1);
end
lnNOR_v=log(qqq);  %the ln(m) of 1922+tt to 1922+tt+qq-1 *age 0-90
MeanNOR2=mean(lnNOR_v,2); % the mean of log(mx) at a given age across test set.
lnNOR_cv=zeros(ag,qq);
for m=1:qq
    lnNOR_cv(:,m)=lnNOR_v(:,m)-MeanNOR2;
end
lnNOR_cv; %Centering ln(mx) of 1922+tt to 1922+tt+qq-1
%Remove variables
clear NOR1
clear xxx 
clear zzz
clear qqq
clear n
clear i
clear yyy

SWE1=load('SWE.mat');%Sweden
yyy=SWE1.SWE;
xxx=yyy(:,total_d);   %total death rate in column 5
swe1=table2array(xxx); %table to array
zzz=zeros(ag,tt); %(91 ages * tt years)
for n=1:tt   
    zzz(1:ag,n)=swe1(1+22089+(n-1)*111:22089+111*n-20,1);
end
lnSWE_t=log(zzz);  %ln(m) 1922-1922+tt-1 year with 0-90 age
MeanSWE1=mean(lnSWE_t,2); % the mean of log(mx) at a given age across test set.
lnSWE_ct=zeros(ag,tt);
for m=1:tt
    lnSWE_ct(:,m)=lnSWE_t(:,m)-MeanSWE1;
end
lnSWE_ct; %Centering ln(mx) from 1922-1922+tt-1  
% Validation Set (1922+tt to 1922+tt+qq-1)
qqq =zeros(ag,qq); %(91 ages * qq years) Validation Matrix
for n=1:qq
    qqq(1:ag,n)=swe1(1+22089+111*tt+(n-1)*111:22089+111*tt+111*n-20,1);
end
lnSWE_v=log(qqq);  %the ln(m) of 1922+tt to 1922+tt+qq-1 *age 0-90
MeanSWE2=mean(lnSWE_v,2); % the mean of log(mx) at a given age across test set.
lnSWE_cv=zeros(ag,qq);
for m=1:qq
    lnSWE_cv(:,m)=lnSWE_v(:,m)-MeanSWE2;
end
lnSWE_cv; %Centering ln(mx) of 1922+tt to 1922+tt+qq-1
%Remove variables
clear SWE1
clear xxx 
clear zzz
clear qqq
clear n
clear i
clear yyy

FRATNP1=load('FRATNP.mat');%France
yyy=FRATNP1.FRATNP;
xxx=yyy(:,total_d);   %total death rate in column 5
fra1=table2array(xxx); %table to array
zzz=zeros(ag,tt); %(91 ages * tt years)
for n=1:tt   
    zzz(1:ag,n)=fra1(1+14874+(n-1)*111:14874+111*n-20,1);
end
lnFRA_t=log(zzz);  %ln(m) 1922-1922+tt-1 year with 0-90 age
MeanFRA1=mean(lnFRA_t,2); % the mean of log(mx) at a given age across test set.
lnFRA_ct=zeros(ag,tt);
for m=1:tt
    lnFRA_ct(:,m)=lnFRA_t(:,m)-MeanFRA1;
end
lnFRA_ct; %Centering ln(mx) from 1922-1922+tt-1  
% Validation Set (1922+tt to 1922+tt+qq-1)
qqq =zeros(ag,qq); %(91 ages * qq years) Validation Matrix
for n=1:qq
    qqq(1:ag,n)=fra1(1+14874+111*tt+(n-1)*111:14874+111*tt+111*n-20,1);
end
lnFRA_v=log(qqq);  %the ln(m) of 1922+tt to 1922+tt+qq-1 *age 0-90
MeanFRA2=mean(lnFRA_v,2); % the mean of log(mx) at a given age across test set.
lnFRA_cv=zeros(ag,qq);
for m=1:qq
    lnFRA_cv(:,m)=lnFRA_v(:,m)-MeanFRA2;
end
lnFRA_cv; %Centering ln(mx) of 1922+tt to 1922+tt+qq-1
%Remove variables
clear FRA1
clear xxx 
clear zzz
clear qqq
clear n
clear i
clear yyy

%                 STEP 2: Construct Target Tensor
%We use centering log(mx) to consruct Tensor.
TT=zeros(ag,tt,cc);
TT(:,:,1)=lnSWE_ct;
TT(:,:,2)=lnUK_ct;
TT(:,:,3)=lnCHE_ct;
TT(:,:,4)=lnDNK_ct;
TT(:,:,5)=lnESP_ct;
TT(:,:,6)=lnFIN_ct;
TT(:,:,7)=lnITA_ct;
TT(:,:,8)=lnNLD_ct;
TT(:,:,9)=lnNOR_ct;
TT(:,:,10)=lnFRA_ct;
TT=tensor(TT);

%                 STEP 3: Tensor Decomposition
%We will test all these 300 seeds's RMSE to determine the best fit.
TT_t=zeros(ag,tt,cc);
TT_t(:,:,1)=lnSWE_t;
TT_t(:,:,2)=lnUK_t;
TT_t(:,:,3)=lnCHE_t;%This is the ln(mx) before centering data.
TT_t(:,:,4)=lnDNK_t;
TT_t(:,:,5)=lnESP_t;
TT_t(:,:,6)=lnFIN_t;
TT_t(:,:,7)=lnITA_t;
TT_t(:,:,8)=lnNLD_t;
TT_t(:,:,9)=lnNOR_t;
TT_t(:,:,10)=lnFRA_t;

Mean1=[MeanSWE1,MeanUK1,MeanCHE1,MeanDNK1,MeanESP1,MeanFIN1,MeanITA1,MeanNLD1,MeanNOR1,MeanFRA1];
%'rank' is the hyperparameter in this part
% We try 'ALS optimization for CP and Tucker tensor decompositions'
CP=cell(10,1); % CP method

%Step 3.1: Choose the best rank 1 fit

CP1=cell(seeds,1);
parfor rr=1:seeds;
    rng(rr);
    CP1(rr)= {parafac_als(TT,1)};
end

RMSE_rank1=zeros(seeds,1);
CP_rank1c=zeros(ag,tt,cc);
for rr=1:seeds;
    CP_rank1=double(CP1{rr,1});% imput every cp-rank1 decomposition after centering
    for ci=1:cc
    for m=1:tt
    CP_rank1c(:,m,ci)=CP_rank1(:,m,ci)+Mean1(:,ci);
    end
    end
    %we reverse-centering the rank 1 tesnor,which is log(mx) now.
    
    CP_rank1d=CP_rank1c-TT_t;
    RMSE_rank1(rr,1)=rms3d(CP_rank1d); %300 times of RMSEs
end
clear m;
[Min1,Seq1]=min(RMSE_rank1);% the best 

CP(1)=CP1(Seq1); % Now we get best rank1 decomposition.

%Step 3.2: Choose the best rank 2 fit

CP2=cell(seeds,1);
parfor rr=1:seeds;
    rng(rr);
    CP2(rr)= {parafac_als(TT,2)};
end

CP_rank2c=zeros(ag,tt,cc);
RMSE_rank2=zeros(seeds,1);
for rr=1:seeds;
    CP_rank2=double(CP2{rr,1});% imput every cp-rank2 decomposition after centering
    for ci=1:cc
    for m=1:tt
    CP_rank2c(:,m,ci)=CP_rank2(:,m,ci)+Mean1(:,ci);
    end
    end
    %we reverse-centering the rank 2 tesnor,which is log(mx) now.
    
    CP_rank2d=CP_rank2c-TT_t;
    RMSE_rank2(rr,1)=rms3d(CP_rank2d); %300 times of RMSEs
end
clear m;
[Min2,Seq2]=min(RMSE_rank2);% the best one 

CP(2)=CP2(Seq2); % Now we get best rank 2 decomposition.

%Step 3.3: Choose the best rank 3 fit
CP3=cell(seeds,1);
parfor rr=1:seeds;
    rng(rr);
    CP3(rr)= {parafac_als(TT,3)};
end
CP_rank3c=zeros(ag,tt,cc);
RMSE_rank3=zeros(seeds,1);
for rr=1:seeds;
    CP_rank3=double(CP3{rr,1});% imput every cp-rank3 decomposition after centering
    for ci=1:cc
    for m=1:tt
    CP_rank3c(:,m,ci)=CP_rank3(:,m,ci)+Mean1(:,ci);
    end
    end
    %we reverse-centering the rank 3 tesnor,which is log(mx) now.
    
    CP_rank3d=CP_rank3c-TT_t;
    RMSE_rank3(rr,1)=rms3d(CP_rank3d); %300 times of RMSEs
end
clear m;
[Min3,Seq3]=min(RMSE_rank3);% the best one 

CP(3)=CP3(Seq3); % Now we get best rank 3 decomposition.

%Step 3.4: Choose the best rank 4 fit
CP4=cell(seeds,1);
parfor rr=1:seeds;
    rng(rr);
    CP4(rr)= {parafac_als(TT,4)};
end
CP_rank4c=zeros(ag,tt,cc);
RMSE_rank4=zeros(seeds,1);
for rr=1:seeds;
    CP_rank4=double(CP4{rr,1});% imput every cp-rank4 decomposition after centering
   
    for ci=1:cc
    for m=1:tt
    CP_rank4c(:,m,ci)=CP_rank4(:,m,ci)+Mean1(:,ci);
    end
    end
    %we reverse-centering the rank 4 tesnor,which is log(mx) now.
    
    CP_rank4d=CP_rank4c-TT_t;
    RMSE_rank4(rr,1)=rms3d(CP_rank4d); % seeds times of RMSEs
end
clear m;
[Min4,Seq4]=min(RMSE_rank4);% the best one 

CP(4)=CP4(Seq4); % Now we get best rank 4 decomposition.
corl=CP4{1,1}.U{2}

%step 3.5 choose the best rank 5 fit
CP5=cell(seeds,1);
parfor rr=1:seeds;
    rng(rr);
    CP5(rr)= {parafac_als(TT,5)};
end
CP_rank5c=zeros(ag,tt,cc);
RMSE_rank5=zeros(seeds,1);
for rr=1:seeds;
    CP_rank5=double(CP5{rr,1});% imput every cp-rank5 decomposition after centering
   
    for ci=1:cc
    for m=1:tt
    CP_rank5c(:,m,ci)=CP_rank5(:,m,ci)+Mean1(:,ci);
    end
    end
    %we reverse-centering the rank 5 tesnor,which is log(mx) now.
    CP_rank5d=CP_rank5c-TT_t;
    RMSE_rank5(rr,1)=rms3d(CP_rank5d); % seeds times of RMSEs
end
clear m;
[Min5,Seq5]=min(RMSE_rank5);% the best one

CP(5)=CP5(Seq5); % Now we get best rank 5 decomposition.

%step 3.6 choose the best rank 6 fit
CP6=cell(seeds,1);
parfor rr=1:seeds;
    rng(rr);
    CP6(rr)= {parafac_als(TT,6)};
end
CP_rank6c=zeros(ag,tt,cc);
RMSE_rank6=zeros(seeds,1);
for rr=1:seeds;
    CP_rank6=double(CP6{rr,1});% imput every cp-rank6 decomposition after centering
   
    for ci=1:cc
    for m=1:tt
    CP_rank6c(:,m,ci)=CP_rank6(:,m,ci)+Mean1(:,ci);
    end
    end
    %we reverse-centering the rank 6 tesnor,which is log(mx) now.
    
    CP_rank6d=CP_rank6c-TT_t;
    RMSE_rank6(rr,1)=rms3d(CP_rank6d); % seeds times of RMSEs
end
clear m;
[Min6,Seq6]=min(RMSE_rank6);% the best one 

CP(6)=CP6(Seq6); % Now we get best rank 6 decomposition.

%step 3.7 choose the best rank 7 fit
CP7=cell(seeds,1);
parfor rr=1:seeds;
    rng(rr);
    CP7(rr)= {parafac_als(TT,7)};
end
CP_rank7c=zeros(ag,tt,cc);

RMSE_rank7=zeros(seeds,1);
for rr=1:seeds;
    CP_rank7=double(CP7{rr,1});% imput every cp-rank7 decomposition after centering
   
    for ci=1:cc
    for m=1:tt
    CP_rank7c(:,m,ci)=CP_rank7(:,m,ci)+Mean1(:,ci);
    end
    end
    %we reverse-centering the rank 7 tesnor,which is log(mx) now.
    
    CP_rank7d=CP_rank7c-TT_t;
    RMSE_rank7(rr,1)=rms3d(CP_rank7d); % seeds times of RMSEs
end
clear m;
[Min7,Seq7]=min(RMSE_rank7);% the best 

CP(7)=CP7(Seq7); % Now we get best rank 7 decomposition.

%step 3.8 choose the best rank 8 fit
CP8=cell(seeds,1);
parfor rr=1:seeds;
    rng(rr);
    CP8(rr)= {parafac_als(TT,8)};
end
CP_rank8c=zeros(ag,tt,cc);
RMSE_rank8=zeros(seeds,1);
for rr=1:seeds;
    CP_rank8=double(CP8{rr,1});% imput every cp-rank8 decomposition after centering
   
    for ci=1:cc
    for m=1:tt
    CP_rank8c(:,m,ci)=CP_rank8(:,m,ci)+Mean1(:,ci);
    end
    end
    %we reverse-centering the rank 8 tesnor,which is log(mx) now.
    
    CP_rank8d=CP_rank8c-TT_t;
    RMSE_rank8(rr,1)=rms3d(CP_rank8d); % seeds times of RMSEs
end
clear m;
[Min8,Seq8]=min(RMSE_rank8);% the best
CP(8)=CP8(Seq8); % Now we get best rank 8 decomposition.

%step 3.9 choose the best rank 9 fit
CP9=cell(seeds,1);
parfor rr=1:seeds;
    rng(rr);
    CP9(rr)= {parafac_als(TT,9)};
end
CP_rank9c=zeros(ag,tt,cc);
RMSE_rank9=zeros(seeds,1);
for rr=1:seeds;
    CP_rank9=double(CP9{rr,1});% imput every cp-rank9 decomposition after centering
   
    for ci=1:cc
    for m=1:tt
    CP_rank9c(:,m,ci)=CP_rank9(:,m,ci)+Mean1(:,ci);
    end
    end
    %we reverse-centering the rank 9 tesnor,which is log(mx) now. 
    CP_rank9d=CP_rank9c-TT_t;
    RMSE_rank9(rr,1)=rms3d(CP_rank9d); % seeds times of RMSEs
end
clear m;
[Min9,Seq9]=min(RMSE_rank9);% the best one 
CP(9)=CP9(Seq9); % Now we get best rank 9 decomposition.

%step 3.10 choose the best rank 10 fit
CP10=cell(seeds,1);
parfor rr=1:seeds;
    rng(rr);
    CP10(rr)= {parafac_als(TT,10)};
end
CP_rank10c=zeros(ag,tt,cc);
RMSE_rank10=zeros(seeds,1);
for rr=1:seeds;
    CP_rank10=double(CP10{rr,1});% imput every cp-rank10 decomposition after centering
   
    for ci=1:cc
    for m=1:tt
    CP_rank10c(:,m,ci)=CP_rank10(:,m,ci)+Mean1(:,ci);
    end
    end
    %we reverse-centering the rank 10 tesnor,which is log(mx) now. 
    CP_rank10d=CP_rank10c-TT_t;
    RMSE_rank10(rr,1)=rms3d(CP_rank10d); % seeds times of RMSEs
end
clear m;
[Min10,Seq10]=min(RMSE_rank10);% the best one 
CP(10)=CP10(Seq10); % Now we get best rank 10 decomposition.

%                STEP 4: Random walk to the year vectors
cp1_year=arimafunction(CP{1,1}.U{2},qq,arimapdq); %the forecast year vectors,random walk
cp2_year=arimafunction(CP{2,1}.U{2},qq,arimapdq); %the forecast year vectors,random walk
cp3_year=arimafunction(CP{3,1}.U{2},qq,arimapdq); %the forecast year vectors,random walk
cp4_year=arimafunction(CP{4,1}.U{2},qq,arimapdq); %the forecast year vectors,random walk
cp5_year=arimafunction(CP{5,1}.U{2},qq,arimapdq); %the forecast year vectors,random walk
cp6_year=arimafunction(CP{6,1}.U{2},qq,arimapdq); %the forecast year vectors,random walk
cp7_year=arimafunction(CP{7,1}.U{2},qq,arimapdq); %the forecast year vectors,random walk
cp8_year=arimafunction(CP{8,1}.U{2},qq,arimapdq); %the forecast year vectors,random walk
cp9_year=arimafunction(CP{9,1}.U{2},qq,arimapdq); %the forecast year vectors,random walk
cp10_year=arimafunction(CP{10,1}.U{2},qq,arimapdq); %the forecast year vectors,random walk

%    STEP 5:Forecast qq years mortality by tensor and arma model
%rank 1 forecast
cp1_age=CP{1,1}.U{1};  % age vector of rank 1 tensor
cp1_country=CP{1,1}.U{3}; %country  vector of rank 1 tensor
cp1_lambda=CP{1,1}.lambda; %the lambda
FC_cp1=ktensor(cp1_lambda,{cp1_age,cp1_year,cp1_country});%the forecast tensor in ktensor format
FCT_cp1=full(FC_cp1); %Transfer ktensor to tensor
FCD_cp1=zeros(ag,qq,cc); %transfer tensor to double format
for i=1:cc
FCD_cp1(:,:,i)=FCT_cp1(:,:,i);
end %remember this is centering log(mx).
clear i;
FCD_cp1b=FCD_cp1;
for ci=1:cc
for m=1:qq
    FCD_cp1b(:,m,ci)=FCD_cp1(:,m,ci)+Mean1(:,ci);
end
end
    FCD_cp1b; %we reverse-centering the rank 1 tesnor,which is log(mx) now.

%rank 2 forecast;
cp2_age=CP{2,1}.U{1};% age vectors of rank 2 tensor
cp2_country=CP{2,1}.U{3}; %country vectors of rank2 tensor
cp2_lambda=CP{2,1}.lambda; %the lambda
FC_cp2=ktensor(cp2_lambda,{cp2_age,cp2_year,cp2_country}); %forecast tensor in ktensor format
FCT_cp2=full(FC_cp2); %Transfer to tensor format
FCD_cp2=zeros(ag,qq,cc); %transfer tensor to double format
for i=1:cc
FCD_cp2(:,:,i)=FCT_cp2(:,:,i);
end %remember this is centering log(mx).

FCD_cp2b=FCD_cp2;
for ci=1:cc
for m=1:qq
    FCD_cp2b(:,m,ci)=FCD_cp2(:,m,ci)+Mean1(:,ci);
end
end
    FCD_cp2b; %we reverse-centering the rank 2 tesnor,which is log(mx) now.
clear i
clear m

%Rank 3 forecast
cp3_age=CP{3,1}.U{1};% age vectors of rank 3 tensor
cp3_country=CP{3,1}.U{3}; %country vectors of rank 3 tensor
cp3_lambda=CP{3,1}.lambda; %the lambda
FC_cp3=ktensor(cp3_lambda,{cp3_age,cp3_year,cp3_country}); %forecast tensor in ktensor format
FCT_cp3=full(FC_cp3); %Transfer to tensor format
FCD_cp3=zeros(ag,qq,cc); %transfer tensor to double format
for i=1:cc
FCD_cp3(:,:,i)=FCT_cp3(:,:,i);
end %remember this is centering log(mx).

FCD_cp3b=FCD_cp3;
for ci=1:cc
for m=1:qq
    FCD_cp3b(:,m,ci)=FCD_cp3(:,m,ci)+Mean1(:,ci);
end
end
    FCD_cp3b; %we reverse-centering the rank 3 tesnor,which is log(mx) now.
clear i
clear m

%rank 4 forecast
cp4_age=CP{4,1}.U{1};% age vectors of rank 4 tensor
cp4_country=CP{4,1}.U{3}; %country vectors of rank 4 tensor
cp4_lambda=CP{4,1}.lambda; %the lambda
FC_cp4=ktensor(cp4_lambda,{cp4_age,cp4_year,cp4_country}); %forecast tensor in ktensor format
FCT_cp4=full(FC_cp4); %Transfer to tensor format
FCD_cp4=zeros(ag,qq,cc); %transfer tensor to double format
for i=1:cc
FCD_cp4(:,:,i)=FCT_cp4(:,:,i);
end %remember this is centering log(mx).

FCD_cp4b=FCD_cp4;
for ci=1:cc
for m=1:qq
    FCD_cp4b(:,m,ci)=FCD_cp4(:,m,ci)+Mean1(:,ci);
end
end
    FCD_cp4b; %we reverse-centering the rank 4 tesnor,which is log(mx) now.
clear i
clear m

%rank 5 forecast
cp5_age=CP{5,1}.U{1};% age vectors of rank 5 tensor
cp5_country=CP{5,1}.U{3}; %country vectors of rank 5 tensor
cp5_lambda=CP{5,1}.lambda; %the lambda
FC_cp5=ktensor(cp5_lambda,{cp5_age,cp5_year,cp5_country}); %forecast tensor in ktensor format
FCT_cp5=full(FC_cp5); %Transfer ktensor to tensor format
FCD_cp5=zeros(ag,qq,cc); %transfer tensor to double format
for i=1:cc
FCD_cp5(:,:,i)=FCT_cp5(:,:,i);
end %remember this is after centering log(mx).
% reverse centering
FCD_cp5b=FCD_cp5;
for ci=1:cc
for m=1:qq
    FCD_cp5b(:,m,ci)=FCD_cp5(:,m,ci)+Mean1(:,ci);
end
end
    FCD_cp5b; %we reverse-centering the rank 5 tesnor,which is log(mx) now.
clear i
clear m

%rank 6 forecast
cp6_age=CP{6,1}.U{1};% age vectors of rank 6 tensor
cp6_country=CP{6,1}.U{3}; %country vectors of rank 6 tensor
cp6_lambda=CP{6,1}.lambda; %the lambda
FC_cp6=ktensor(cp6_lambda,{cp6_age,cp6_year,cp6_country}); %forecast tensor in ktensor format
FCT_cp6=full(FC_cp6); %Transfer ktensor to tensor format
FCD_cp6=zeros(ag,qq,cc); %transfer tensor to double format
for i=1:cc
FCD_cp6(:,:,i)=FCT_cp6(:,:,i);
end %remember this is after centering log(mx).
% reverse centering
FCD_cp6b=FCD_cp6;
for ci=1:cc
for m=1:qq
    FCD_cp6b(:,m,ci)=FCD_cp6(:,m,ci)+Mean1(:,ci);
end
end
    FCD_cp6b; %we reverse-centering the rank 6 tesnor,which is log(mx) now.
clear i
clear m

%rank 7 forecast
cp7_age=CP{7,1}.U{1};% age vectors of rank 7 tensor
cp7_country=CP{7,1}.U{3}; %country vectors of rank 7 tensor
cp7_lambda=CP{7,1}.lambda; %the lambda
FC_cp7=ktensor(cp7_lambda,{cp7_age,cp7_year,cp7_country}); %forecast tensor in ktensor format
FCT_cp7=full(FC_cp7); %Transfer ktensor to tensor format
FCD_cp7=zeros(ag,qq,cc); %transfer tensor to double format
for i=1:cc
FCD_cp7(:,:,i)=FCT_cp7(:,:,i);
end %remember this is after centering log(mx).
% reverse centering
FCD_cp7b=FCD_cp7;
for ci=1:cc
for m=1:qq
    FCD_cp7b(:,m,ci)=FCD_cp7(:,m,ci)+Mean1(:,ci);
end
end
    FCD_cp7b; %we reverse-centering the rank 7 tesnor,which is log(mx) now.
clear i
clear m

%rank 8 forecast
cp8_age=CP{8,1}.U{1};% age vectors of rank 8 tensor
cp8_country=CP{8,1}.U{3}; %country vectors of rank 8 tensor
cp8_lambda=CP{8,1}.lambda; %the lambda
FC_cp8=ktensor(cp8_lambda,{cp8_age,cp8_year,cp8_country}); %forecast tensor in ktensor format
FCT_cp8=full(FC_cp8); %Transfer ktensor to tensor format
FCD_cp8=zeros(ag,qq,cc); %transfer tensor to double format
for i=1:cc
FCD_cp8(:,:,i)=FCT_cp8(:,:,i);
end %remember this is after centering log(mx).
% reverse centering
FCD_cp8b=FCD_cp8;
for ci=1:cc
for m=1:qq
    FCD_cp8b(:,m,ci)=FCD_cp8(:,m,ci)+Mean1(:,ci);
end
end
    FCD_cp8b; %we reverse-centering the rank 8 tesnor,which is log(mx) now.
clear i
clear m

%rank 9 forecast
cp9_age=CP{9,1}.U{1};% age vectors of rank 9 tensor
cp9_country=CP{9,1}.U{3}; %country vectors of rank 9 tensor
cp9_lambda=CP{9,1}.lambda; %the lambda
FC_cp9=ktensor(cp9_lambda,{cp9_age,cp9_year,cp9_country}); %forecast tensor in ktensor format
FCT_cp9=full(FC_cp9); %Transfer ktensor to tensor format
FCD_cp9=zeros(ag,qq,cc); %transfer tensor to double format
for i=1:cc
FCD_cp9(:,:,i)=FCT_cp9(:,:,i);
end %remember this is after centering log(mx).
% reverse centering
FCD_cp9b=FCD_cp9;
for ci=1:cc
for m=1:qq
    FCD_cp9b(:,m,ci)=FCD_cp9(:,m,ci)+Mean1(:,ci);
end
end
    FCD_cp9b; %we reverse-centering the rank 9 tesnor,which is log(mx) now.
clear i
clear m

%rank 10 forecast
cp10_age=CP{10,1}.U{1};% age vectors of rank 10 tensor
cp10_country=CP{10,1}.U{3}; %country vectors of rank 10 tensor
cp10_lambda=CP{10,1}.lambda; %the lambda
FC_cp10=ktensor(cp10_lambda,{cp10_age,cp10_year,cp10_country}); %forecast tensor in ktensor format
FCT_cp10=full(FC_cp10); %Transfer ktensor to tensor format
FCD_cp10=zeros(ag,qq,cc); %transfer tensor to double format
for i=1:cc
FCD_cp10(:,:,i)=FCT_cp10(:,:,i);
end %remember this is after centering log(mx).
% reverse centering
FCD_cp10b=FCD_cp10;
for ci=1:cc
for m=1:qq
    FCD_cp10b(:,m,ci)=FCD_cp10(:,m,ci)+Mean1(:,ci);
end
end
    FCD_cp10b; %we reverse-centering the rank 10 tesnor,which is log(mx) now.
clear i
clear m


%  STEP 6: Compare forecasts and determine the hyperparameters
%impute real observations 
Realob=zeros(ag,qq,cc);
Realob(:,:,1)=lnSWE_v;
Realob(:,:,2)=lnUK_v;
Realob(:,:,3)=lnCHE_v;%they are log(mx) before the centering process
Realob(:,:,4)=lnDNK_v;
Realob(:,:,5)=lnESP_v;
Realob(:,:,6)=lnFIN_v;
Realob(:,:,7)=lnITA_v;
Realob(:,:,8)=lnNLD_v;
Realob(:,:,9)=lnNOR_v;
Realob(:,:,10)=lnFRA_v;


%RMSE of rank 1 forecast
RMSE_cpr1d= Realob-FCD_cp1b; %get the difference between the real and the forecst 1.
RMSE_cpr1=rms3d(RMSE_cpr1d); %the RMSE of rank 1 tensor add the mean back

%RMSE of rank 2 forecast
RMSE_cpr2d= Realob-FCD_cp2b; %get the difference between the real and the forecst 2.
RMSE_cpr2=rms3d(RMSE_cpr2d); %the RMSE of rank 2 tensor add the mean back

%RMSE of rank 3 forecast
RMSE_cpr3d= Realob-FCD_cp3b; %get the difference between the real and the forecst 3.
RMSE_cpr3=rms3d(RMSE_cpr3d); %the RMSE of rank 3 tensor add the mean back

%RMSE of rank 4 forecast
RMSE_cpr4d= Realob-FCD_cp4b; %get the difference between the real and the forecst 4.
RMSE_cpr4=rms3d(RMSE_cpr4d); %the RMSE of rank 4 tensor add the mean back

%RMSE of rank 5 foreacst
RMSE_cpr5d= Realob-FCD_cp5b; %get the difference between the real and the forecst 5.
RMSE_cpr5=rms3d(RMSE_cpr5d); %the RMSE of rank 5 tensor add the mean back

%RMSE of rank 6 foreacst
RMSE_cpr6d= Realob-FCD_cp6b; %get the difference between the real and the forecst 6.
RMSE_cpr6=rms3d(RMSE_cpr6d); %the RMSE of rank 6 tensor add the mean back

%RMSE of rank 7 foreacst
RMSE_cpr7d= Realob-FCD_cp7b; %get the difference between the real and the forecst 7.
RMSE_cpr7=rms3d(RMSE_cpr7d); %the RMSE of rank 7 tensor add the mean back

%RMSE of rank 8 foreacst
RMSE_cpr8d= Realob-FCD_cp8b; %get the difference between the real and the forecst 8.
RMSE_cpr8=rms3d(RMSE_cpr8d); %the RMSE of rank 8 tensor add the mean back

%RMSE of rank 9 foreacst
RMSE_cpr9d= Realob-FCD_cp9b; %get the difference between the real and the forecst 9.
RMSE_cpr9=rms3d(RMSE_cpr9d); %the RMSE of rank 9 tensor add the mean back

%RMSE of rank 10 foreacst
RMSE_cpr10d= Realob-FCD_cp10b; %get the difference between the real and the forecst 10.
RMSE_cpr10=rms3d(RMSE_cpr10d); %the RMSE of rank 10 tensor add the mean back
RMSE_ranks= zeros(10);
RMSE_ranks = [RMSE_cpr1,RMSE_cpr2,RMSE_cpr3,RMSE_cpr4,RMSE_cpr5,RMSE_cpr6,RMSE_cpr7,RMSE_cpr8,RMSE_cpr9,RMSE_cpr10]
[M, rank_validation] = min(RMSE_ranks);
%******************************************************************************************************************

%Step 7: We use rank=rank_validation as our final model and do the forecast the total:

%We need to use 1922 to 1922+tt+qq-1 as our test set and fit a new Target tensor.
lnSWE_new=[lnSWE_t,lnSWE_v]; %log(mx) of SWE from 
lnUK_new=[lnUK_t,lnUK_v]; %log(mx) of UK 
lnCHE_new=[lnCHE_t,lnCHE_v]; 
lnDNK_new=[lnDNK_t,lnDNK_v]; 
lnESP_new=[lnESP_t,lnESP_v]; 
lnFIN_new=[lnFIN_t,lnFIN_v]; 
lnITA_new=[lnITA_t,lnITA_v]; 
lnNLD_new=[lnNLD_t,lnNLD_v]; 
lnNOR_new=[lnNOR_t,lnNOR_v]; 
lnFRA_new=[lnFRA_t,lnFRA_v]; 
% get means of the training+val
MeanSWE=mean(lnSWE_new,2);
MeanUK=mean(lnUK_new,2);
MeanCHE=mean(lnCHE_new,2);% the Mean of the year vectors
MeanDNK=mean(lnDNK_new,2);% the Mean of the year vectors
MeanESP=mean(lnESP_new,2);% the Mean of the year vectors
MeanFIN=mean(lnFIN_new,2);% the Mean of the year vectors
MeanITA=mean(lnITA_new,2);% the Mean of the year vectors
MeanNLD=mean(lnNLD_new,2);% the Mean of the year vectors
MeanNOR=mean(lnNOR_new,2);% the Mean of the year vectors
MeanFRA=mean(lnFRA_new,2);% the Mean of the year vectors

%
lnSWE_newc=zeros(ag,tt+qq); %Centering log(mx) 
lnUK_newc=zeros(ag,tt+qq); %Centering log(mx) 
lnCHE_newc=zeros(ag,tt+qq); %Centering log(mx) 
lnDNK_newc=zeros(ag,tt+qq); %Centering log(mx) 
lnESP_newc=zeros(ag,tt+qq); %Centering log(mx) 
lnFIN_newc=zeros(ag,tt+qq); %Centering log(mx) 
lnITA_newc=zeros(ag,tt+qq); %Centering log(mx) 
lnNLD_newc=zeros(ag,tt+qq); %Centering log(mx) 
lnNOR_newc=zeros(ag,tt+qq); %Centering log(mx) 
lnFRA_newc=zeros(ag,tt+qq); %Centering log(mx) 
for m=1:(tt+qq);
    lnSWE_newc(:,m)=lnSWE_new(:,m)-MeanSWE; 
    lnUK_newc(:,m)=lnUK_new(:,m)-MeanUK;
    lnCHE_newc(:,m)=lnCHE_new(:,m)-MeanCHE;
    lnDNK_newc(:,m)=lnDNK_new(:,m)-MeanDNK;
    lnESP_newc(:,m)=lnESP_new(:,m)-MeanESP;
    lnFIN_newc(:,m)=lnFIN_new(:,m)-MeanFIN;
    lnITA_newc(:,m)=lnITA_new(:,m)-MeanITA;
    lnNLD_newc(:,m)=lnNLD_new(:,m)-MeanNLD;
    lnNOR_newc(:,m)=lnNOR_new(:,m)-MeanNOR;
    lnFRA_newc(:,m)=lnFRA_new(:,m)-MeanFRA;       
end
clear m
%
%We use after centering log(mx) to consruct new Tensor.
TT_new=zeros(ag,tt+qq,cc); %ag*(tt+qq)*cc tensor
TT_new(:,:,1)=lnSWE_newc;
TT_new(:,:,2)=lnUK_newc;
TT_new(:,:,3)=lnCHE_newc;
TT_new(:,:,4)=lnDNK_newc;
TT_new(:,:,5)=lnESP_newc;
TT_new(:,:,6)=lnFIN_newc;
TT_new(:,:,7)=lnITA_newc;
TT_new(:,:,8)=lnNLD_newc;
TT_new(:,:,9)=lnNOR_newc;
TT_new(:,:,10)=lnFRA_newc;
TT_new=tensor(TT_new); %This is ag*(tt+qq)*cc after centering tensor

% The real log(mx) before centering
TT_real=zeros(ag,tt+qq,cc);
TT_real(:,:,1)=lnSWE_new;
TT_real(:,:,2)=lnUK_new;
TT_real(:,:,3)=lnCHE_new; %This is ag*tt+qq*cc before centering data
TT_real(:,:,4)=lnDNK_new;
TT_real(:,:,5)=lnESP_new;
TT_real(:,:,6)=lnFIN_new;
TT_real(:,:,7)=lnITA_new;
TT_real(:,:,8)=lnNLD_new;
TT_real(:,:,9)=lnNOR_new;
TT_real(:,:,10)=lnFRA_new;

%Mean of tt+qq
Mean_tq=[MeanSWE,MeanUK,MeanCHE,MeanDNK,MeanESP,MeanFIN,MeanITA,MeanNLD,MeanNOR,MeanFRA]

%choose best fit decomposition(training+validation)                   
 CP_newf=cell(trys,1);
parfor rr=1:trys;
    rng(rr); %random seeds
    CP_newf(rr)= {parafac_als(TT_new,rank_validation)};% rank is determined as rank_validation
end
CP_ranknew=zeros(ag,tt+qq,cc);
RMSE_ranknew=zeros(trys,1);
for rr=1:trys;
    CP_ranknewf=double(CP_newf{rr,1});% imput every decomposition after centering
   for ci=1:cc
    for m=1:tt+qq
    CP_ranknew(:,m,ci)=CP_ranknewf(:,m,ci)+Mean_tq(:,ci);
    end
    end
    %we reverse-centering the rank_validation tesnor,which is log(mx) now.
    CP_ranknewd=CP_ranknew-TT_real;
    RMSE_ranknew(rr,1)=rms3d(CP_ranknewd); %300 times of RMSEs
end
clear m;
[Minnew,Seqnew]=min(RMSE_ranknew);% the best one i
CP_new=CP_newf(Seqnew); % Now we get best decomposition.

% random walk to the year vector
% rank_validation ARIMA
cpnew_age=CP_new{1}.U{1};  % age vector of rank_validation tensor
cpnew_country=CP_new{1}.U{3}; %country  vector of rank_validation tensor
cpnew_lambda=CP_new{1}.lambda; %the lambda
cpnew_year=arimafunction(CP_new{1}.U{2},qq,arimapdq); %the forecast year vectors,random walk

FC_cpnew=ktensor(cpnew_lambda,{cpnew_age,cpnew_year,cpnew_country});%the forecast tensor in ktensor format
FCT_cpnew=full(FC_cpnew); %Transfer ktensor to tensor
FCD_cpnew=zeros(ag,qq,cc); %transfer tensor to double format
for i=1:cc
FCD_cpnew(:,:,i)=FCT_cpnew(:,:,i);
end         %remember this is centering log(mx).
%
FCD_cpnewb=FCD_cpnew;
for ci=1:cc;
for m=1:qq;
    FCD_cpnewb(:,m,ci)=FCD_cpnew(:,m,ci)+Mean_tq(:,ci);
end
end%we reverse-centering the tesnor,which is log(mx) now.
clear m
FCD_cpnewb; 
%this is our forecast tensor which will be compared with the Lee-Carter model.

%Step 8: Compare with the Lee-Carter model by RMSE(1922+tt+qq)
%Import real observations from   of 10 countries.
SWE_new =zeros(ag,qq); %(91 ages * qq years) Forecasting Matrix
UK_new =zeros(ag,qq);
CHE_new =zeros(ag,qq);
DNK_new =zeros(ag,qq);
ESP_new =zeros(ag,qq);
FIN_new =zeros(ag,qq);
ITA_new =zeros(ag,qq);
NLD_new =zeros(ag,qq);
NOR_new =zeros(ag,qq);
FRA_new =zeros(ag,qq);
for n=1:qq
    SWE_new(:,n)=swe1(1+22089+111*(tt+qq)+(n-1)*111:22089+111*(tt+qq)+111*n-20,1);
    UK_new(:,n)=uk1(1+3108+111*(tt+qq)+(n-1)*111:3108+111*(tt+qq)+111*n-20,1);
    CHE_new(:,n)=che1(1+8214+111*(tt+qq)+(n-1)*111:8214+111*(tt+qq)+111*n-20,1);
    DNK_new(:,n)=dnk1(1+12765+111*(tt+qq)+(n-1)*111:12765+111*(tt+qq)+111*n-20,1);
    ESP_new(:,n)=esp1(1+4662+111*(tt+qq)+(n-1)*111:4662+111*(tt+qq)+111*n-20,1);
    FIN_new(:,n)=fin1(1+7992+111*(tt+qq)+(n-1)*111:7992+111*(tt+qq)+111*n-20,1);
    ITA_new(:,n)=ita1(1+8658+111*(tt+qq)+(n-1)*111:8658+111*(tt+qq)+111*n-20,1);
    NLD_new(:,n)=nld1(1+11100+111*(tt+qq)+(n-1)*111:11100+111*(tt+qq)+111*n-20,1);
    NOR_new(:,n)=nor1(1+11544+111*(tt+qq)+(n-1)*111:11544+111*(tt+qq)+111*n-20,1);
    FRA_new(:,n)=fra1(1+14874+111*(tt+qq)+(n-1)*111:14874+111*(tt+qq)+111*n-20,1);
end

%
%impute real log(observations) 
Realob_new=zeros(ag,qq,cc);
Realob_new(:,:,1)=log(SWE_new);
Realob_new(:,:,2)=log(UK_new);
Realob_new(:,:,3)=log(CHE_new); %they are log(mx) before the centering process
Realob_new(:,:,4)=log(DNK_new);
Realob_new(:,:,5)=log(ESP_new);
Realob_new(:,:,6)=log(FIN_new);
Realob_new(:,:,7)=log(ITA_new);
Realob_new(:,:,8)=log(NLD_new);
Realob_new(:,:,9)=log(NOR_new);
Realob_new(:,:,10)=log(FRA_new);

%RMSE of forecast and country specific RMSE.
RMSE_cpnew= Realob_new-FCD_cpnewb; %get the difference between the real and the forecst 
RMSE_new=rms3d(RMSE_cpnew); %the RMSE of tensor add the mean back.

% RMSE of each country
RMSE_countries=zeros(cc,1);
for ci=1:cc;
    gap=Realob_new(:,:,ci)-FCD_cpnewb(:,:,ci);
    RMSE_countries(ci)=rms3d(gap);
end

