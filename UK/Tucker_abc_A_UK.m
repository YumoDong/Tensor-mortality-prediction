    %STEP 1: Import Data
% Initial settings 
 qq=5; %length of valiadation/testing set (year), you can change it to change the forecasting horizon. 

T_begin=1922; %the begining year
T_end=2016; % the ending year
tt=T_end-T_begin+1-2*qq; %length of test set
ag=91-20; %length of ages
cc=2; %number of genders
rms3d= dsp.RMS; %RMS fuction
rms3d.Dimension='ALL'; %calculate the whole tensor's RMS
arimapdq=arima(0,1,0); %make a ARIMA(0,1,0) model for later use
seeds=50; %This this is used to determine best fit of final model in validation
male_d=4;
female_d=3;
total_d=5;  %3,4,5 represent different death rate 
xx=15;

%In this file,CHE represents UK female, DNK represents UK male
% UK data
CHE1=load('GBR.mat'); 
yyy=CHE1.GBRNP;
xxx=yyy(:,female_d);   %UK female
che1=table2array(xxx); %table to array
zzz=zeros(ag,tt); 
for n=1:tt   
    zzz(1:ag,n)=che1(1+0+20+(n-1)*111:0+111*n-20,1);
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
    qqq(1:ag,n)=che1(1+20+111*tt+(n-1)*111:111*tt+111*n-20,1);
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

DNK1=load('GBR.mat'); %UK male
yyy=DNK1.GBRNP;
xxx=yyy(:,male_d);   %male
dnk1=table2array(xxx); %table to array
zzz=zeros(ag,tt);
for n=1:tt   
    zzz(1:ag,n)=dnk1(1+0+20+(n-1)*111:0+111*n-20,1);
end
lnDNK_t=log(zzz);  %ln(m)
MeanDNK1=mean(lnDNK_t,2); % the mean of log(mx) at a given age across test set.
lnDNK_ct=zeros(ag,tt);
for m=1:tt
    lnDNK_ct(:,m)=lnDNK_t(:,m)-MeanDNK1;
end
lnDNK_ct; %Centering ln(mx) 
% Validation Set 
qqq =zeros(ag,qq); %Validation Matrix
for n=1:qq
    qqq(1:ag,n)=dnk1(1+20+111*tt+(n-1)*111:111*tt+111*n-20,1);
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


%                 STEP 2: Construct Target Tensor
%We use centering log(mx) to consruct Tensor.
TT=zeros(ag,tt,cc);
TT(:,:,1)=lnCHE_ct;
TT(:,:,2)=lnDNK_ct;
TT=tensor(TT);

%                 STEP 3: Tensor Decomposition

%We will test all these seeds's RMSE to determine the best fit.
TT_t=zeros(ag,tt,cc);
TT_t(:,:,1)=lnCHE_t;%This is the ln(mx) before centering data.
TT_t(:,:,2)=lnDNK_t;
Mean1=[MeanCHE1,MeanDNK1];

CP=cell(xx,xx,cc); %  cell of decomposition: 15*15*2 
% Tucker decomposition start
CP1=cell(xx,xx,seeds);
CP_rank1c=zeros(ag,tt,cc);
RMSE_rank1=zeros(xx,xx,seeds);

%c=1
parfor i=1:xx;
    for j=1:xx,
    for rr=1:seeds;
        rng(rr);
    CP1(i,j,rr)= {tucker_als(TT,[i j 1])};
    end
    end
end  %get decomposition result

for i=1:xx;
    for j=1:xx;
    for rr=1:seeds;
    CP_rank1=double(CP1{i,j,rr});% imput every cp-rank2 decomposition after centering
    
    for ci=1:cc
     CP_rank1c(:,:,ci)=CP_rank1(:,:,ci)+Mean1(:,ci);
    end
    %we reverse-centering the rank 2 tesnor,which is log(mx) now.
    CP_rank1d=CP_rank1c-TT_t;
    RMSE_rank1(i,j,rr)=rms3d(CP_rank1d); %50 times of RMSEs
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

% c= 2
CP2=cell(xx,xx,seeds);
CP_rank2c=zeros(ag,tt,cc);
RMSE_rank2=zeros(xx,xx,seeds);
parfor i=1:xx;
    for j=1:xx,
    for rr=1:seeds;
        rng(rr);
    CP2(i,j,rr)= {tucker_als(TT,[i j 2])};
    end
    end
end  %get decomposition result
for i=1:xx;
    for j=1:xx;
    for rr=1:seeds;
    CP_rank2=double(CP2{i,j,rr});% imput every cp-rank2 decomposition after centering
    
    for ci=1:cc
     CP_rank2c(:,:,ci)=CP_rank2(:,:,ci)+Mean1(:,ci);
    end
    %we reverse-centering the rank 2 tesnor,which is log(mx) now.
    CP_rank2d=CP_rank2c-TT_t;
    RMSE_rank2(i,j,rr)=rms3d(CP_rank2d); %50 times of RMSEs
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


%    STEP 4: Forecast qq years mortality by tensor and arma model
FCD=cell(xx,xx,cc); %15*15*2 forecasts
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
cp1_country=CP{i,j,pp}.U{3}; %gender vectors of rank tensor
cp1_year=double(arma_year{i,j,pp}); %the forecast rank year vectors
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
Realob(:,:,1)=lnCHE_v;%they are log(mx) before the centering process
Realob(:,:,2)=lnDNK_v;%they are log(mx) before the centering process

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
