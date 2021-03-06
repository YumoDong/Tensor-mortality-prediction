  %STEP 1: Import Data
% Initial settings
qq=5; %length of valiadation/testing set (year), you can change it to change the forecasting horizon.

T_begin=1922; %the begining year
T_end=2014; % the ending year
tt=T_end-T_begin+1-qq*2; %length of test set
sa=0; % starting age(age range: sa to 90)
ag=91-sa; %length of ages
cc=10; %number of countries
rms3d= dsp.RMS; %RMS fuction
rms3d.Dimension='ALL'; %calculate the whole tensor's RMS
arimapdq=arima(0,1,0); %make a random walk for later use
seeds=50; %This this is used to determine best fit of final model in test set
male_d=4;
female_d=3;
total_d=5;  %3,4,5 represent different death rate 
xx=15; %the upperbound of the vector year and age.

% UK data
GBR1=load('GBR.mat');
UK=GBR1.GBRNP;
uk=UK(:,total_d);   %UK'total death rate in column 5
uk1=table2array(uk); %table to array
UK_t=zeros(ag,tt); %(ag ages * tt years)
for n=1:tt   
    UK_t(1:ag,n)=uk1(sa+1+0+(n-1)*111:0+111*n-20,1);

end
lnUK_t=log(UK_t);  %ln(m) 
MeanUK1=mean(lnUK_t,2); % the mean of log(mx) at a given age across test set.
lnUK_ct=zeros(ag,tt);
for m=1:tt
    lnUK_ct(:,m)=lnUK_t(:,m)-MeanUK1;
end
lnUK_ct; %Centering ln(mx) 
% Validation Set 
UK_v =zeros(ag,qq); % Validation Matrix
for n=1:qq
    UK_v(1:ag,n)=uk1(sa+1+111*tt+(n-1)*111:111*tt+111*n-20,1);
end
lnUK_v=log(UK_v);  %the ln(m) 
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
zzz=zeros(ag,tt);
for n=1:tt   
    zzz(1:ag,n)=che1(sa+1+5106+(n-1)*111:5106+111*n-20,1);
end
lnCHE_t=log(zzz);  %ln(m) 
MeanCHE1=mean(lnCHE_t,2); % the mean of log(mx) at a given age across test set.
lnCHE_ct=zeros(ag,tt);
for m=1:tt
    lnCHE_ct(:,m)=lnCHE_t(:,m)-MeanCHE1;
end
lnCHE_ct; %Centering ln(mx) 
% Validation Set 
qqq =zeros(ag,qq); % Validation Matrix
for n=1:qq
    qqq(1:ag,n)=che1(sa+1+5106+111*tt+(n-1)*111:5106+111*tt+111*n-20,1);
end
lnCHE_v=log(qqq);  %the ln(m) 
MeanCHE2=mean(lnCHE_v,2); % the mean of log(mx) at a given age across test set.
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

DNK1=load('DNK.mat'); %Denmark
yyy=DNK1.DNK1;
xxx=yyy(:,total_d);   %total death rate in column 5
dnk1=table2array(xxx); %table to array
zzz=zeros(ag,tt);
for n=1:tt   
    zzz(1:ag,n)=dnk1(sa+1+9657+(n-1)*111:9657+111*n-20,1);
end
lnDNK_t=log(zzz);  %ln(m) 
MeanDNK1=mean(lnDNK_t,2); % the mean of log(mx) at a given age across test set.
lnDNK_ct=zeros(ag,tt);
for m=1:tt
    lnDNK_ct(:,m)=lnDNK_t(:,m)-MeanDNK1;
end
lnDNK_ct; %Centering ln(mx) 
% Validation Set 
qqq =zeros(ag,qq); % Validation Matrix
for n=1:qq
    qqq(1:ag,n)=dnk1(sa+1+9657+111*tt+(n-1)*111:9657+111*tt+111*n-20,1);
end
lnDNK_v=log(qqq);  %the ln(m) 
MeanDNK2=mean(lnDNK_v,2); % the mean of log(mx) at a given age across test set.
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

ESP1=load('ESP.mat'); %Spain
yyy=ESP1.ESP;
xxx=yyy(:,total_d);   %total death rate in column 5
esp1=table2array(xxx); %table to array
zzz=zeros(ag,tt); 
for n=1:tt   
    zzz(1:ag,n)=esp1(sa+1+1554+(n-1)*111:1554+111*n-20,1);
end
lnESP_t=log(zzz);  %ln(m)
MeanESP1=mean(lnESP_t,2); % the mean of log(mx) at a given age across test set.
lnESP_ct=zeros(ag,tt);
for m=1:tt
    lnESP_ct(:,m)=lnESP_t(:,m)-MeanESP1;
end
lnESP_ct; %Centering ln(mx) 
% Validation Set 
qqq =zeros(ag,qq); %Validation Matrix
for n=1:qq
    qqq(1:ag,n)=esp1(sa+1+1554+111*tt+(n-1)*111:1554+111*tt+111*n-20,1);
end
lnESP_v=log(qqq);  %the ln(m) 
MeanESP2=mean(lnESP_v,2); % the mean of log(mx) at a given age across test set.
lnESP_cv=zeros(ag,qq);
for m=1:qq
    lnESP_cv(:,m)=lnESP_v(:,m)-MeanESP2;
end
lnESP_cv; %Centering ln(mx)
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
zzz=zeros(ag,tt); 
for n=1:tt   
    zzz(1:ag,n)=fin1(sa+1+4884+(n-1)*111:4884+111*n-20,1);
end
lnFIN_t=log(zzz);  %ln(m) 
MeanFIN1=mean(lnFIN_t,2); % the mean of log(mx) at a given age across test set.
lnFIN_ct=zeros(ag,tt);
for m=1:tt
    lnFIN_ct(:,m)=lnFIN_t(:,m)-MeanFIN1;
end
lnFIN_ct; %Centering ln(mx)
% Validation Set 
qqq =zeros(ag,qq); %(Validation Matrix
for n=1:qq
    qqq(1:ag,n)=fin1(sa+1+4884+111*tt+(n-1)*111:4884+111*tt+111*n-20,1);
end
lnFIN_v=log(qqq);  %the ln(m) 
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
zzz=zeros(ag,tt); %
for n=1:tt   
    zzz(1:ag,n)=ita1(sa+1+5550+(n-1)*111:5550+111*n-20,1);
end
lnITA_t=log(zzz);  %ln(m) 
MeanITA1=mean(lnITA_t,2); % the mean of log(mx) at a given age across test set.
lnITA_ct=zeros(ag,tt);
for m=1:tt
    lnITA_ct(:,m)=lnITA_t(:,m)-MeanITA1;
end
lnITA_ct; %Centering ln(mx) 
% Validation Set 
qqq =zeros(ag,qq); %Validation Matrix
for n=1:qq
    qqq(1:ag,n)=ita1(sa+1+5550+111*tt+(n-1)*111:5550+111*tt+111*n-20,1);
end
lnITA_v=log(qqq);  %the ln(m) 
MeanITA2=mean(lnITA_v,2); % the mean of log(mx) at a given age across test set.
lnITA_cv=zeros(ag,qq);
for m=1:qq
    lnITA_cv(:,m)=lnITA_v(:,m)-MeanITA2;
end
lnITA_cv; %Centering ln(mx) 
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
zzz=zeros(ag,tt); 
for n=1:tt   
    zzz(1:ag,n)=nld1(sa+1+7992+(n-1)*111:7992+111*n-20,1);
end
lnNLD_t=log(zzz);  %ln(m) 
MeanNLD1=mean(lnNLD_t,2); % the mean of log(mx) at a given age across test set.
lnNLD_ct=zeros(ag,tt);
for m=1:tt
    lnNLD_ct(:,m)=lnNLD_t(:,m)-MeanNLD1;
end
lnNLD_ct; %Centering ln(mx) 
% Validation Set 
qqq =zeros(ag,qq); %(91 ages * qq years) Validation Matrix
for n=1:qq
    qqq(1:ag,n)=nld1(sa+1+7992+111*tt+(n-1)*111:7992+111*tt+111*n-20,1);
end
lnNLD_v=log(qqq);  %the ln(m) 
MeanNLD2=mean(lnNLD_v,2); % the mean of log(mx) at a given age across test set.
lnNLD_cv=zeros(ag,qq);
for m=1:qq
    lnNLD_cv(:,m)=lnNLD_v(:,m)-MeanNLD2;
end
lnNLD_cv; %Centering ln(mx)
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
zzz=zeros(ag,tt); 
for n=1:tt   
    zzz(1:ag,n)=nor1(sa+1+8436+(n-1)*111:8436+111*n-20,1);
end
lnNOR_t=log(zzz);  %ln(m)
MeanNOR1=mean(lnNOR_t,2); % the mean of log(mx) at a given age across test set.
lnNOR_ct=zeros(ag,tt);
for m=1:tt
    lnNOR_ct(:,m)=lnNOR_t(:,m)-MeanNOR1;
end
lnNOR_ct; %Centering ln(mx) 
% Validation Set 
qqq =zeros(ag,qq); %Validation Matrix
for n=1:qq
    qqq(1:ag,n)=nor1(sa+1+8436+111*tt+(n-1)*111:8436+111*tt+111*n-20,1);
end
lnNOR_v=log(qqq);  %the ln(m) 
MeanNOR2=mean(lnNOR_v,2); % the mean of log(mx) at a given age across test set.
lnNOR_cv=zeros(ag,qq);
for m=1:qq
    lnNOR_cv(:,m)=lnNOR_v(:,m)-MeanNOR2;
end
lnNOR_cv; %Centering ln(mx) 
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
zzz=zeros(ag,tt); %
for n=1:tt   
    zzz(1:ag,n)=swe1(sa+1+18981+(n-1)*111:18981+111*n-20,1);
end
lnSWE_t=log(zzz);  %ln(m) 
MeanSWE1=mean(lnSWE_t,2); % the mean of log(mx) at a given age across test set.
lnSWE_ct=zeros(ag,tt);
for m=1:tt
    lnSWE_ct(:,m)=lnSWE_t(:,m)-MeanSWE1;
end
lnSWE_ct; %Centering ln(mx) 
% Validation Set
qqq =zeros(ag,qq); % Validation Matrix
for n=1:qq
    qqq(1:ag,n)=swe1(sa+1+18981+111*tt+(n-1)*111:18981+111*tt+111*n-20,1);
end
lnSWE_v=log(qqq);  %the ln(m) 
MeanSWE2=mean(lnSWE_v,2); % the mean of log(mx) at a given age across test set.
lnSWE_cv=zeros(ag,qq);
for m=1:qq
    lnSWE_cv(:,m)=lnSWE_v(:,m)-MeanSWE2;
end
lnSWE_cv; %Centering ln(mx)
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
zzz=zeros(ag,tt); 
for n=1:tt   
    zzz(1:ag,n)=fra1(sa+1+11766+(n-1)*111:11766+111*n-20,1);
end
lnFRA_t=log(zzz);  %ln(m) 
MeanFRA1=mean(lnFRA_t,2); % the mean of log(mx) at a given age across test set.
lnFRA_ct=zeros(ag,tt);
for m=1:tt
    lnFRA_ct(:,m)=lnFRA_t(:,m)-MeanFRA1;
end
lnFRA_ct; %Centering ln(mx) 
% Validation Set 
qqq =zeros(ag,qq); % Validation Matrix
for n=1:qq
    qqq(1:ag,n)=fra1(sa+1+11766+111*tt+(n-1)*111:11766+111*tt+111*n-20,1);
end
lnFRA_v=log(qqq);  %the ln(m) 
MeanFRA2=mean(lnFRA_v,2); % the mean of log(mx) at a given age across test set.
lnFRA_cv=zeros(ag,qq);
for m=1:qq
    lnFRA_cv(:,m)=lnFRA_v(:,m)-MeanFRA2;
end
lnFRA_cv; %Centering ln(mx) 
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
%'rank' is the hyperparameter in this part
% We try 'ALS optimization for Tucker tensor decompositions'
CP=cell(xx,xx,10); %  cell of decomposition: 15*15*10
%It contains the Tucker decompositions result of all ranks.

%Step 3.1: Choose the best rank  fit
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

% Tucker decomposition start
%c=1
CP1=cell(xx,xx,seeds);
CP_rank1c=zeros(ag,tt,cc);
RMSE_rank1=zeros(xx,xx,seeds);

parfor i=1:xx;
    for j=1:xx,
    for rr=1:seeds;
    rng(rr)
    CP1(i,j,rr)= {tucker_als(TT,[i j 1])};
    end
    end
end  %get decomposition result

for i=1:xx;
    for j=1:xx;
    for rr=1:seeds;
    CP_rank1=double(CP1{i,j,rr});% imput every decomposition after centering
    
    for ci=1:cc
     CP_rank1c(:,:,ci)=CP_rank1(:,:,ci)+Mean1(:,ci);
    end
    %we reverse-centering the tesnor,which is log(mx) now.
    CP_rank1d=CP_rank1c-TT_t;
    RMSE_rank1(i,j,rr)=rms3d(CP_rank1d); %300 times of RMSEs
    end
    end
end
clear m;
[Min1,Seq1]=min(RMSE_rank1,[],3);
for i=1:xx;
for  j=1:xx;
CP(i,j,1)=CP1(i,j,Seq1(i,j)); % Now we get best rank [a b 1] decomposition.
end
end

% c = 2
CP2=cell(xx,xx,seeds);
CP_rank2c=zeros(ag,tt,cc);
RMSE_rank2=zeros(xx,xx,seeds);
parfor i=1:xx;
    for j=1:xx,
    for rr=1:seeds;
     rng(rr)
    CP2(i,j,rr)= {tucker_als(TT,[i j 2])};
    end
    end
end  %get decomposition result
for i=1:xx;
    for j=1:xx;
    for rr=1:seeds;
    CP_rank2=double(CP2{i,j,rr});% imput every decomposition after centering
    
    for ci=1:cc
     CP_rank2c(:,:,ci)=CP_rank2(:,:,ci)+Mean1(:,ci);
    end
    %we reverse-centering the tesnor,which is log(mx) now.
    CP_rank2d=CP_rank2c-TT_t;
    RMSE_rank2(i,j,rr)=rms3d(CP_rank2d); %300 times of RMSEs
    end
    end
end
clear m;
[Min2,Seq2]=min(RMSE_rank2,[],3);
for i=1:xx;
for  j=1:xx;
CP(i,j,2)=CP2(i,j,Seq2(i,j)); % Now we get best rank [a b 2] decomposition.
end
end

% c = 3
CP3=cell(xx,xx,seeds);
CP_rank3c=zeros(ag,tt,cc);
RMSE_rank3=zeros(xx,xx,seeds);
parfor i=1:xx;
    for j=1:xx,
    for rr=1:seeds;
        rng(rr)
    CP3(i,j,rr)= {tucker_als(TT,[i j 3])};
    end
    end
end  %get decomposition result
for i=1:xx;
    for j=1:xx;
    for rr=1:seeds;
    CP_rank3=double(CP3{i,j,rr});% imput every decomposition after centering
    
    for ci=1:cc
     CP_rank3c(:,:,ci)=CP_rank3(:,:,ci)+Mean1(:,ci);
    end
    %we reverse-centering the tesnor,which is log(mx) now.
    CP_rank3d=CP_rank3c-TT_t;
    RMSE_rank3(i,j,rr)=rms3d(CP_rank3d); %300 times of RMSEs
    end
    end
end
clear m;
[Min3,Seq3]=min(RMSE_rank3,[],3);
for i=1:xx;
for  j=1:xx;
CP(i,j,3)=CP3(i,j,Seq3(i,j)); % Now we get best rank [a b 3] decomposition.
end
end

% c= 4
CP4=cell(xx,xx,seeds);
CP_rank4c=zeros(ag,tt,cc);
RMSE_rank4=zeros(xx,xx,seeds);
parfor i=1:xx;
    for j=1:xx,
    for rr=1:seeds;
        rng(rr)
    CP4(i,j,rr)= {tucker_als(TT,[i j 4])};
    end
    end
end  %get decomposition result
for i=1:xx;
    for j=1:xx;
    for rr=1:seeds;
    CP_rank4=double(CP4{i,j,rr});% imput every  decomposition after centering
    
    for ci=1:cc
     CP_rank4c(:,:,ci)=CP_rank4(:,:,ci)+Mean1(:,ci);
    end
    %we reverse-centering the tesnor,which is log(mx) now.
    CP_rank4d=CP_rank4c-TT_t;
    RMSE_rank4(i,j,rr)=rms3d(CP_rank4d); %300 times of RMSEs
    end
    end
end
clear m;
[Min4,Seq4]=min(RMSE_rank4,[],3);
for i=1:xx;
for  j=1:xx;
CP(i,j,4)=CP4(i,j,Seq4(i,j)); % Now we get best rank [a b 4] decomposition.
end
end

% c = 5

CP5=cell(xx,xx,seeds);
CP_rank5c=zeros(ag,tt,cc);
RMSE_rank5=zeros(xx,xx,seeds);
parfor i=1:xx;
    for j=1:xx,
    for rr=1:seeds;
        rng(rr)
    CP5(i,j,rr)= {tucker_als(TT,[i j 5])};
    end
    end
end  %get decomposition result
for i=1:xx;
    for j=1:xx;
    for rr=1:seeds;
    CP_rank5=double(CP5{i,j,rr});% imput every  decomposition after centering
    
    for ci=1:cc
     CP_rank5c(:,:,ci)=CP_rank5(:,:,ci)+Mean1(:,ci);
    end
    %we reverse-centering the  tesnor,which is log(mx) now.
    CP_rank5d=CP_rank5c-TT_t;
    RMSE_rank5(i,j,rr)=rms3d(CP_rank5d); %300 times of RMSEs
    end
    end
end
clear m;
[Min5,Seq5]=min(RMSE_rank5,[],3);
for i=1:xx;
for  j=1:xx;
CP(i,j,5)=CP5(i,j,Seq5(i,j)); % Now we get best rank [a b 5] decomposition.
end
end

%c=6
CP6=cell(xx,xx,seeds);
CP_rank6c=zeros(ag,tt,cc);
RMSE_rank6=zeros(xx,xx,seeds);
parfor i=1:xx;
    for j=1:xx,
    for rr=1:seeds;
        rng(rr)
    CP6(i,j,rr)= {tucker_als(TT,[i j 6])};
    end
    end
end  %get decomposition result
for i=1:xx;
    for j=1:xx;
    for rr=1:seeds;
    CP_rank6=double(CP6{i,j,rr});% imput every  decomposition after centering
    
    for ci=1:cc
     CP_rank6c(:,:,ci)=CP_rank6(:,:,ci)+Mean1(:,ci);
    end
    %we reverse-centering the tesnor,which is log(mx) now.
    CP_rank6d=CP_rank6c-TT_t;
    RMSE_rank6(i,j,rr)=rms3d(CP_rank6d); %300 times of RMSEs
    end
    end
end
clear m;
[Min6,Seq6]=min(RMSE_rank6,[],3);
for i=1:xx;
for  j=1:xx;
CP(i,j,6)=CP6(i,j,Seq6(i,j)); % Now we get best rank [a b 6] decomposition.
end
end

%c=7
CP7=cell(xx,xx,seeds);
CP_rank7c=zeros(ag,tt,cc);
RMSE_rank7=zeros(xx,xx,seeds);
parfor i=1:xx;
    for j=1:xx,
    for rr=1:seeds;
        rng(rr)
    CP7(i,j,rr)= {tucker_als(TT,[i j 7])};
    end
    end
end  %get decomposition result
for i=1:xx;
    for j=1:xx;
    for rr=1:seeds;
    CP_rank7=double(CP7{i,j,rr});% imput every decomposition after centering
    
    for ci=1:cc
     CP_rank7c(:,:,ci)=CP_rank7(:,:,ci)+Mean1(:,ci);
    end
    %we reverse-centering tesnor,which is log(mx) now.
    CP_rank7d=CP_rank7c-TT_t;
    RMSE_rank7(i,j,rr)=rms3d(CP_rank7d); %300 times of RMSEs
    end
    end
end
clear m;
[Min7,Seq7]=min(RMSE_rank7,[],3);
for i=1:xx;
for  j=1:xx;
CP(i,j,7)=CP7(i,j,Seq7(i,j)); % Now we get best rank [a b 7] decomposition.
end
end

% c=8
CP8=cell(xx,xx,seeds);
CP_rank8c=zeros(ag,tt,cc);
RMSE_rank8=zeros(xx,xx,seeds);
parfor i=1:xx;
    for j=1:xx,
    for rr=1:seeds;
        rng(rr)
    CP8(i,j,rr)= {tucker_als(TT,[i j 8])};
    end
    end
end  %get decomposition result
for i=1:xx;
    for j=1:xx;
    for rr=1:seeds;
    CP_rank8=double(CP8{i,j,rr});% imput every decomposition after centering
    
    for ci=1:cc
     CP_rank8c(:,:,ci)=CP_rank8(:,:,ci)+Mean1(:,ci);
    end
    %we reverse-centering the tesnor,which is log(mx) now.
    CP_rank8d=CP_rank8c-TT_t;
    RMSE_rank8(i,j,rr)=rms3d(CP_rank8d); %300 times of RMSEs
    end
    end
end
clear m;
[Min8,Seq8]=min(RMSE_rank8,[],3);
for i=1:xx;
for  j=1:xx;
CP(i,j,8)=CP8(i,j,Seq8(i,j)); % Now we get best rank [a b 8] decomposition.
end
end

%c= 9
CP9=cell(xx,xx,seeds);
CP_rank9c=zeros(ag,tt,cc);
RMSE_rank9=zeros(xx,xx,seeds);
parfor i=1:xx;
    for j=1:xx,
    for rr=1:seeds;
        rng(rr)
    CP9(i,j,rr)= {tucker_als(TT,[i j 9])};
    end
    end
end  %get decomposition result
for i=1:xx;
    for j=1:xx;
    for rr=1:seeds;
    CP_rank9=double(CP9{i,j,rr});% imput every decomposition after centering
    
    for ci=1:cc
     CP_rank9c(:,:,ci)=CP_rank9(:,:,ci)+Mean1(:,ci);
    end
    %we reverse-centering the tesnor,which is log(mx) now.
    CP_rank9d=CP_rank9c-TT_t;
    RMSE_rank9(i,j,rr)=rms3d(CP_rank9d); %300 times of RMSEs
    end
    end
end
clear m;
[Min9,Seq9]=min(RMSE_rank9,[],3);
for i=1:xx;
for  j=1:xx;
CP(i,j,9)=CP9(i,j,Seq9(i,j)); % Now we get best rank [a b 9] decomposition.
end
end

%c= 10
CP10=cell(xx,xx,seeds);
CP_rank10c=zeros(ag,tt,cc);
RMSE_rank10=zeros(xx,xx,seeds);
parfor i=1:xx;
    for j=1:xx,
    for rr=1:seeds;
        rng(rr)
    CP10(i,j,rr)= {tucker_als(TT,[i j 10])};
    end
    end
end  %get decomposition result
for i=1:xx;
    for j=1:xx;
    for rr=1:seeds;
    CP_rank10=double(CP10{i,j,rr});% imput every decomposition after centering
    
    for ci=1:cc
     CP_rank10c(:,:,ci)=CP_rank10(:,:,ci)+Mean1(:,ci);
    end
    %we reverse-centering the  tesnor,which is log(mx) now.
    CP_rank10d=CP_rank10c-TT_t;
    RMSE_rank10(i,j,rr)=rms3d(CP_rank10d); %300 times of RMSEs
    end
    end
end
clear m;
[Min10,Seq10]=min(RMSE_rank10,[],3);
for i=1:xx;
for  j=1:xx;
CP(i,j,10)=CP10(i,j,Seq10(i,j)); % Now we get best rank [a b 10] decomposition.
end
end

%    STEP 4: Forecast qq years mortality by tensor and random walk model

FCD=cell(xx,xx,cc); %15*15*10 forecasts
arma_year=cell(xx,xx,cc); %armia for year
parfor pp=1:cc;
for j=1:xx;
for i=1:xx;
arma_year(i,j,pp)={arimafunction(CP{i,j,pp}.U{2},qq,arimapdq)}; %the forecast year vectors
end
end
end

parfor pp=1:cc;
for j=1:xx;
for i=1:xx;
cp1_age=CP{i,j,pp}.U{1};% age vectors of rank (i,j,pp) tensor
cp1_country=CP{i,j,pp}.U{3}; %country vectors of rank tensor
cp1_year=double(arma_year{i,j,pp}); %the forecast year vectors
cp1_lambda=CP{i,j,pp}.core; %the lambda
FC_cp1=ttm(cp1_lambda,{cp1_age,cp1_year,cp1_country}); %forecast tensor in ktensor format
FCT_cp1=full(FC_cp1); %Transfer to tensor format
FCD_cp1=zeros(ag,qq,cc); %transfer tensor to double format

for ii=1:cc
FCD_cp1(:,:,ii)=FCT_cp1(:,:,ii);
end %remember this is centering log(mx).

FCD_cp1b=FCD_cp1;
for ci=1:cc
    FCD_cp1b(:,:,ci)=FCD_cp1(:,:,ci)+Mean1(:,ci);
end

FCD(i,j,pp)={FCD_cp1b}; %we reverse-centering the tesnor,which is log(mx) now.
end
end
end

%  STEP 5: Compare forecasts and determine the hyperparameters
%impute real observations 
Realob=zeros(ag,qq,cc);
Realob(:,:,1)=lnSWE_v;
Realob(:,:,2)=lnUK_v;
Realob(:,:,3)=lnCHE_v;%they are log(mx) before the centering process
Realob(:,:,4)=lnDNK_v;%they are log(mx) before the centering process
Realob(:,:,5)=lnESP_v;%they are log(mx) before the centering process
Realob(:,:,6)=lnFIN_v;%they are log(mx) before the centering process
Realob(:,:,7)=lnITA_v;%they are log(mx) before the centering process
Realob(:,:,8)=lnNLD_v;%they are log(mx) before the centering process
Realob(:,:,9)=lnNOR_v;%they are log(mx) before the centering process
Realob(:,:,10)=lnFRA_v;%they are log(mx) before the centering process

%RMSE of all rank (a b c)  forecast
RMSE_tk=zeros(xx,xx,cc);
parfor pp=1:cc;
for j=1:xx;
for i=1:xx;
RMSE_cpr1d= Realob-FCD{i,j,pp}; %get the difference between the real and the forecst.
RMSE_tk(i,j,pp)=rms3d(RMSE_cpr1d); %the RMSE of Tucker(a b c) add the mean back 
end
end
end
minimum = min(min(min(RMSE_tk)));
[Tucker_abc_a, mixed]=find(RMSE_tk==minimum);
if rem(mixed,xx)==0
    Tucker_abc_b = xx;
else
    Tucker_abc_b = rem(mixed,xx);
end

Tucker_abc_c= fix((mixed-1)/xx)+1;

%RMSE of all rank (a a b)  forecast
RMSE_tk_aab=zeros(xx,cc);
parfor pp=1:cc;
for i=1:xx;
RMSE_cpr1d_aab= Realob-FCD{i,i,pp}; %get the difference between the real and the forecst.
RMSE_tk_aab(i,pp)=rms3d(RMSE_cpr1d_aab); %the RMSE of Tucker(a a b) add the mean back 
end
end
minimum = min(min(RMSE_tk_aab));
[Tucker_aab_a,Tucker_aab_b]=find(RMSE_tk_aab==minimum)

%RMSE of all rank (a a a)  forecast
RMSE_tk_aaa=zeros(cc,1);
parfor pp=1:cc;
RMSE_cpr1d_aaa= Realob-FCD{pp,pp,pp}; %get the difference between the real and the forecst.
RMSE_tk_aaa(pp,1)=rms3d(RMSE_cpr1d_aaa); %the RMSE of Tucker(a a a) add the mean back 
end
[M,Tucker_aaa_a] = min(RMSE_tk_aaa);
