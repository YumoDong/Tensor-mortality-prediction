   %STEP 1: Import Data
%Entering some parameters
rank=[2 2 2]; % rank determined in the validation set
qq=20; %length of valiadation set

T_begin=1950; %the begining year
T_end=2014; % the ending year
tt=T_end-T_begin+1-2*qq; %length of test set
ag=91; %length of ages
cc=10; %number of countries
rms3d= dsp.RMS; %RMS fuction
rms3d.Dimension='ALL'; %calculate the whole tensor's RMS
arimapdq=arima(0,1,0); %make a ARIMA(p,d,q) model for later use
trys=300; %This this is used to determine best fit of final model in test set(rank 3)
male_d=4;
female_d=3;
total_d=5;  %3,4,5 represent different death rate 
xx=15;
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
lnUK_cv; %Centering ln(mx) of UK
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

%******************************************************************************************************************

%We use rank[a a b] as our final model and do the forecast the total:

%We need to use 1950 to 1950+tt+qq-1 as our new trainging set and fit a new Target tensor.
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

% get means of the new trainging
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
parfor m=1:(tt+qq);
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
TT_real(:,:,3)=lnCHE_new; %This is ag*(tt+qq)*cc before centering data
TT_real(:,:,4)=lnDNK_new;
TT_real(:,:,5)=lnESP_new;
TT_real(:,:,6)=lnFIN_new;
TT_real(:,:,7)=lnITA_new;
TT_real(:,:,8)=lnNLD_new;
TT_real(:,:,9)=lnNOR_new;
TT_real(:,:,10)=lnFRA_new;

%Mean of tt+qq
Mean_tq=[MeanSWE,MeanUK,MeanCHE,MeanDNK,MeanESP,MeanFIN,MeanITA,MeanNLD,MeanNOR,MeanFRA]
  
%choose best fit rank [a b c] decomposition(training+validation)                              
 CP_newf=cell(trys,1);
parfor rr=1:trys;
    rng(rr) % random seeds 1:trys
    CP_newf(rr)= {tucker_als(TT_new,rank)};
end
CP_ranknew=zeros(ag,tt+qq,cc);
RMSE_ranknew=zeros(trys,1);
for rr=1:trys;
    CP_ranknewf=double(CP_newf{rr,1});% imput every cp-rank6 decomposition after centering
   for ci=1:cc
    for m=1:tt+qq
    CP_ranknew(:,m,ci)=CP_ranknewf(:,m,ci)+Mean_tq(:,ci);
    end
    end
    %we reverse-centering the rank 10 tesnor,which is log(mx) now.
    CP_ranknewd=CP_ranknew-TT_real;
    RMSE_ranknew(rr,1)=rms3d(CP_ranknewd); %300 times of RMSEs
end
clear m;
[Minnew,Seqnew]=min(RMSE_ranknew);% the best one is Seq=9, Min=0.0513
CP_new=CP_newf(Seqnew); % Now we get best rank 3 decomposition.

%rank[a b c] forecast
cpnew_age=CP_new{1}.U{1};  % age vector of tensor
cpnew_country=CP_new{1}.U{3}; %country  vector of tensor
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
%Compare with the Lee-Carter model by RMSE(1950+tt+qq)
%Import real observations of 10 countries.
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
parfor n=1:qq
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
%impute real log(observations) ag*(1950+tt+qq to 2014)*cc
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
RMSE_new=rms3d(RMSE_cpnew); %the RMSE of tensor 
% RMSE of each country
RMSE_countries=zeros(cc,1);
parfor ci=1:cc;
    gap=Realob_new(:,:,ci)-FCD_cpnewb(:,:,ci);
    RMSE_countries(ci)=rms3d(gap);
end
