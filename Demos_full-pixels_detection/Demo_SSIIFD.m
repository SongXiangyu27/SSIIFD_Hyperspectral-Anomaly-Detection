%--------------------------------------------------------------------------
% Run an Anomaly Detection based on Improved-iForest on the given data set.
% (C) Xiangyu Song,  University of Chinese Academy of Sciences
%
% This is a demo of using improved isolation forest-based spectral-spatial 
% anomaly detection framework for anomaly detection in both spectral and spatial domains
%
% If you use this code, please kindly cite the paper: <br>
%  @ARTICLE{9521674,  
%  author={Song, Xiangyu and Aryal, Sunil and Ting, Kai Ming and Liu, Zhen and He, Bin},  
%  journal={IEEE Transactions on Geoscience and Remote Sensing},   
%  title={Spectral¨CSpatial Anomaly Detection of Hyperspectral Data Based on Improved Isolation Forest},   
%  year={2022},  volume={60},  number={},  pages={1-16},  
%  doi={10.1109/TGRS.2021.3104998}}
%
%--------------------------------------------------------------------------

clc; 
close all;
warning off
clear all;

addpath(genpath('iForest'))
addpath(genpath('KPCA'))
addpath(genpath('PTA_LiWei'))
addpath(genpath('./Data'));
addpath(genpath('./Gabor'));
addpath(genpath('./3DROC')); 
addpath(genpath('./CRD')); 
addpath(genpath('./IIFD')); 

%% load data and instation

filename = 'SanDiegoI';   % SanDiegoI SanDiegoII TexasCoast Gulfport                            
load (filename)

%-------------------------imshow-------------------------------
f_show=data(:,:,[27,11,1]);
for i=1:3         
    max_f=max(max(f_show(:,:,i)));
    min_f=min(min(f_show(:,:,i)));
    f_show(:,:,i)=(f_show(:,:,i)-min_f)/(max_f-min_f);
end
figure,imshow(f_show,[]);
figure,imshow(f_show(:,:,1),[]);
figure('name','map'),imshow(map,[]);

%------------------------normalization-----------------------
DataTest =data;

mask = map;
mask2 = map_mx1;
[num_rows, num_cols, N] = size(DataTest);  
for i=1:N
    DataTest(:,:,i) = (DataTest(:,:,i)-min(min(DataTest(:,:,i)))) / (max(max(DataTest(:,:,i))-min(min(DataTest(:,:,i)))));
end
M = num_rows * num_cols;

%% Start testing
mask_reshape = reshape(mask, num_rows*num_cols, 1);
anomaly_map = logical(double(mask_reshape)>0  );
normal_map = logical(double(mask_reshape)==0 );
Data = reshape(DataTest,M,N); 

%% Run iForest imitation
% KPCA Processing
tic;
% X = Data;  % data_mxn
label = map_mx1;
% set kernel function
kernel = Kernel('type', 'gauss', 'width', 2); % 0.5 --> 2
% parameter setting
parameter = struct('application', 'dr', 'dim', 300, 'kernel', kernel); % 2 --> N  300
% build a KPCA object
kpca = KernelPCA(parameter);
% train KPCA model using given data
X_map = kpca.train(Data);
toc;
Data_KPCA = reshape(X_map,M,300);
Data_KPCA = Data_KPCA + abs(min(min(Data_KPCA)));
Data_KPCA_3D = reshape(X_map,num_rows, num_cols, 300);
% %------------------------------iForest---------------------------------------
% parameters for iForest
Methods = 4
NumTree = 512; % number of isolation trees 512  1024  1000
NumSub = ceil(num_rows*num_cols*0.03); % subsample size 256; 
NumDim = size(Data, 2); % do not perform dimension sampling  

tic 
score_KIF = RuniForest(Data_KPCA, NumTree, NumSub, NumDim);  %Data_KPCA   
toc 

score_KIF_show = reshape(score_KIF,[num_rows,num_cols]);
figure('name','iForest'), imshow(score_KIF_show,[]);
NBOXPLOT(:,4)=score_KIF_show(:);

%% Improved-iForest Detector
tic
Methods = 5

tic
[Data_PCA_3D] = pca(Data, 10);
Data_PCA_3D = Data_PCA_3D';
Data_PCA = reshape(Data_PCA_3D(:,1),num_rows, num_cols);
gaborArray = gaborFilterBank(5,8,199,199); % Generates Gabor filter bank; 
[ ~,featureVector] = gaborFeatures(Data_PCA,gaborArray,num_rows,num_cols);   % Extracts Gabor feature vector, 'featureVector', from the image,
Gabor_feature = reshape(featureVector,num_rows*num_cols,40) ;
toc
%-------------------------------------------------------------------------
NumTree = 32;  
K = 5;  

tic
[PMass_spectral,QMass_spectral,score_LIIF, ~, ~] = Local_IIFD(Data, Data, map_mx1, NumTree, num_rows, num_cols);
toc

tic
[PMass_spatial,QMass_spatial,score_GabIIF, AUC, ~] = IIFD_Gabor(Gabor_feature, map_mx1, NumTree, K);
toc

PMass_final = 0.618*PMass_spectral + 0.382*PMass_spatial;  
QMass_final = 0.618*QMass_spectral + 0.382*QMass_spatial;
ratio = PMass_final ./ QMass_final;
scoreraw = mean(ratio, 2);
value = 2 * (log(M-1)+0.5772156649) - 2*(M-1) / M;  
score_SSIIFD = 1 - 2.^(-scoreraw./value);

LIIFD_show = reshape(score_LIIF, num_rows, num_cols);
figure('name','LIIFD');imshow(LIIFD_show,[]);
IFD_spatial_show = reshape(score_GabIIF, num_rows, num_cols);
figure('name','IIF_Gabor');imshow(IFD_spatial_show,[]);
SSIIFD_show = reshape(score_SSIIFD, num_rows, num_cols);
figure('name','SSIIFD'); imshow(SSIIFD_show,[]);
NBOXPLOT(:,6)=SSIIFD_show(:);

%% Traditional RX detector
Methods = 1

tic;
r_RX = RX(Data');  % input: num_dim x num_sam    rx
r_RX_norm = (r_RX-min(r_RX)) / (max(r_RX)-min(r_RX));  %X*=£¨X-Xmin£©/(Xmax-Xmin)
toc;

r_RX_show = reshape(r_RX_norm,[num_rows,num_cols]);
figure('name','RX'), imshow(r_RX_show,[]);
NBOXPLOT(:,1)=r_RX_show(:);

 %% Perform CRD & PTA
% %-----------------------------CRD----------------------------------------------------
Methods = 2

tic
r_CRD_show = Unsupervised_NRS_Detect(data, 11, 7, 1e-6); %13  7  SanDiegoI&II:[31 3] airport4:[27 11] Texas Coast:[41 21] HYDICE:[15 7]  
r_CRD_show = (r_CRD_show-min(r_CRD_show(:)))/(max(r_CRD_show(:))-min(r_CRD_show(:))); %X*=£¨X-Xmin£©/(Xmax-Xmin)
toc   

r_CRD = reshape(r_CRD_show, M, 1);
figure('name','CRD'), imshow(r_CRD_show,[]);
NBOXPLOT(:,2)=r_CRD_show(:);

% %-------------PTA with LTV-norm----------------------------------------------------
Methods = 3

tol1 = 1e-4;
tol2 = 1e-6;
maxiter = 400;
truncate_rank = 1;
alphia = 1.7;  % 1.7
beta = 0.069;   %0.069
tau =0.1;       %0.1
tic;
 [X,S,area_PTA] = AD_Tensor_LILU1(DataTest,alphia,beta,tau,truncate_rank,maxiter,tol1,tol2,normal_map,anomaly_map);
toc

r_PTA_show = sqrt(sum(S.^2,3));
r_PTA_show = (r_PTA_show-min(r_PTA_show(:)))/(max(r_PTA_show(:))-min(r_PTA_show(:)));
r_PTA = reshape(r_PTA_show,M,1);

figure('name','PTA');
imshow(r_PTA_show,[]);
NBOXPLOT(:,3)=r_PTA_show(:);

%% Plot 3DROC and Calculate AUC values
results = [score_SSIIFD, score_KIF, r_PTA, r_CRD, r_RX_norm'];
legend_3DROC = {'SSIIFD', 'KIFD', 'PTA', 'CRD', 'RXD'};
Plot_3DROC_S(results,map_mx1,legend_3DROC,1);

[AUC,AUCnor] = cal_AUC(results,map_mx1,1,1);
AUCod = AUCnor.PFPD + AUCnor.tauPD - AUCnor.tauPF;
AUCsnpr = AUCnor.tauPD./AUCnor.tauPF;
