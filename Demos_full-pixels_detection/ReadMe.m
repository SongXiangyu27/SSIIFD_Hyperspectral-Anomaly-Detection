%% This set of files contains the Matlab code for AD experiments 
%% in the following paper:
%
%  [1] Xiangyu Song, Sunil Aryal, KaiMing Ting, Zhen Liu and Bin He, 
%  Spectral-Spatial Anomaly Detection of Hyperspectral Data Based on Improved Isolation Forest. 
%  IEEE Transactions on Geoscience and Remote Sensing, doi: 10.1109/TGRS.2021.3104998., Aug. 2021.
%
%  Note that: 
%  (1)the demo_ global _IIFD.m is a basic demo of using the proposed improved isolation forest detector for
%  anomaly detection in spectral-domain of original HSI data without pre-segmentation or any other preprocessing. 
%  (2)the demo_ local _IIFD.m is a demo of using IIFD for anomaly detection in spectral-domain with an ERS-based segmentation preprocessing
%  (3)the demo_SSIIFD.m is a demo of using improved isolation forest-based spectral-spatial anomaly detection framework for anomaly detection 
%  in both spectral and spatial domains
%
%  Since the global IIFD is easy to implement and not many parameters need to be adjusted, we recommend the demo_global_IIFD.m to readers
%  when detecting anomalies in other HSIs except ones listed in [1]
%
%% The KIFD are generated by:
%
%  [2] S. Li, K. Zhang, P. Duan, and X. Kang, 
%  Hyperspectral anomaly detection with kernel isolation forest, 
%  IEEE Transactions on Geoscience and Remote Sensing, vol. 58, no. 1, pp. 319�C329, Jan. 2020.
%
%% The MFIFD implimentation is by:
%
%  [3] R. Wang, F. Nie, Z. Wang, F. He and X. Li, 
%  Multiple Features and Isolation Forest-Based Fast Anomaly Detector for Hyperspectral Imagery, 
%  IEEE Transactions on Geoscience and Remote Sensing, vol. 58, no. 9, pp. 6664-6676, Sept. 2020.
%
%% the PTA implitation is by:
%
%  [4] L. Li, W. Li, Y. Qu, C. Zhao, R. Tao and Q. Du, 
%  Prior-Based Tensor Approximation for Anomaly Detection in Hyperspectral Imagery, 
%  IEEE Transactions on Neural Networks and Learning Systems,  doi: 10.1109/TNNLS.2020.3038659.
%
%% The CRD are built by:
%
%  [5] W. Li and Q. Du, 
%  Collaborative Representation for Hyperspectral Anomaly Detection,
%  IEEE Transactions on Geoscience and Remote Sensing, vol. 53, no. 3, pp. 1463-1474, March 2015.
%
%% The RXD is implemented by:
%
%  [6] I. S. Reed and X. Yu, 
%  Adaptive multiple-band CFAR detection of an optical pattern with unknown spectral distribution,
%  IEEE Transactions on Acoustics, Speech, and Signal Processing, vol. 38, no. 10, pp. 1760-1770, Oct. 1990.
%
%% The 3DROC implitation is by:
%
%  [7] C. -I. Chang, 
%  An Effective Evaluation Tool for Hyperspectral Target Detection: 3D Receiver Operating Characteristic Curve Analysis,
%  IEEE Transactions on Geoscience and Remote Sensing, vol. 59, no. 6, pp. 5131-5153, June 2021.
%
%%
%%  If you use this demo, please cite these references 
%%  --- [1][2][3][4][5][6][7] ---
%%
%%  For their BibTeX formats, please see 
%%  --- Reference.bib ---
%%
%%  For their RIS formats, please see 
%%  --- Reference.enl ---
%%  
%%
%
%  Any suggestions and comments are appreciated, and please send them to
%  the authors: 
%  songxiangyu17@mails.ucas.edu.cn; songxiangyu@crdc.com; sunil.aryal@deakin.edu.au; tingkm@nju.edu.cn