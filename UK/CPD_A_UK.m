%STEP 1: Import Data
% Initial settings 
 qq=5; %length of valiadation/testing set (year), you can change it to change the forecasting horizon. 

T_begin=1922; %the begining year
T_end=2016; % the ending year
tt=T_end-T_begin+1-2*qq; %length of test set
sa=20;% starting age
ag=91-sa; %length of ages
cc=2; %number of countries/genders
rng(1);%random seed 1
rms3d= dsp.RMS; %RMS fuction
rms3d.Dimension='ALL'; %calculate the whole tensor's RMS
arimapdq=arima(0,1,0); %make a random model for later use
seeds=300; %this is used to determine best high-rank fit in validation set
trys=300; %This this is used to determine best fit of final model in test set
male_d=4;
female_d=3;
total_d=5;  %3,4,5 represent different death rate 

CHE1=load('GBR.mat'); %CHE represents UK female
yyy=CHE1.GBRNP;
xxx=yyy(:,female_d);   
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

DNK1=load('GBR.mat'); %DNK represents UK male
yyy=DNK1.GBRNP;
xxx=yyy(:,male_d);   
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
qqq =zeros(ag,qq); % Validation Matrix
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
%'rank' is the hyperparameter in this part
% We try 'ALS optimization for CP tensor decompositions
TT_t=zeros(ag,tt,cc);
TT_t(:,:,1)=lnCHE_t;%This is the ln(mx) before centering data.
TT_t(:,:,2)=lnDNK_t;
Mean1=[MeanCHE1,MeanDNK1];

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
Realob(:,:,1)=lnCHE_v;%they are log(mx) before the centering process
Realob(:,:,2)=lnDNK_v;%they are log(mx) before the centering process

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

%Step 7: We use rank=10 as our final model and do the forecast the total:

%We need to use 1922 to 1922+tt+qq-1 as our test set and fit a new Target tensor.
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
%We use after centering log(mx) to consruct new Tensor.
TT_new=zeros(ag,tt+qq,cc); %ag*(tt+qq)*cc tensor
TT_new(:,:,1)=lnCHE_newc;
TT_new(:,:,2)=lnDNK_newc;
TT_new=tensor(TT_new); %This is ag*tt+qq*cc after centering tensor

% The real log(mx) before centering
TT_real=zeros(ag,tt+qq,cc);

TT_real(:,:,1)=lnCHE_new; %This is ag*tt+qq*cc before centering data
TT_real(:,:,2)=lnDNK_new;

%Mean of tt+qq
Mean_tq=[MeanCHE,MeanDNK]
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
[Minnew,Seqnew]=min(RMSE_ranknew);% the best one 
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

%Step 8: Compare with the Lee-Carter model by RMSE
%Import real observations 
 %Forecasting Matrix
CHE_new =zeros(ag,qq);
DNK_new =zeros(ag,qq);
for n=1:qq
    
    CHE_new(:,n)=che1(1+20+111*(tt+qq)+(n-1)*111:111*(tt+qq)+111*n-20,1);
    DNK_new(:,n)=dnk1(1+20+111*(tt+qq)+(n-1)*111:111*(tt+qq)+111*n-20,1);
    
end
%
%impute real log(observations) 
Realob_new=zeros(ag,qq,cc);

Realob_new(:,:,1)=log(CHE_new); %they are log(mx) before the centering process
Realob_new(:,:,2)=log(DNK_new);


%RMSE of gender specific RMSE.
RMSE_cpnew= Realob_new-FCD_cpnewb; %get the difference between the real and the forecst
RMSE_new=rms3d(RMSE_cpnew); %the RMSE of tensor add the mean back.
% RMSE of each country
RMSE_countries=zeros(cc,1);
for ci=1:cc;
    gap=Realob_new(:,:,ci)-FCD_cpnewb(:,:,ci);
    RMSE_countries(ci)=rms3d(gap);
end

               
