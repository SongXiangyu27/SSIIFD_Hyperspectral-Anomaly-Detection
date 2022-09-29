function [featureVector,feature_all] = gaborFeatures(img,gaborArray,m,n)

% GABORFEATURES extracts the Gabor features of an input image.
% It creates a column vector, consisting of the Gabor features of the input
% image. The feature vectors are normalized to zero mean and unit variance.
%
%
% Inputs:
%       img         :	Matrix of the input image 
%       gaborArray	:	Gabor filters bank created by the function gaborFilterBank
%       d1          :	The factor of downsampling along rows.
%       d2          :	The factor of downsampling along columns.
%               
% Output:
%       featureVector	:   A column vector with length (m*n*u*v)/(d1*d2). 
%                           This vector is the Gabor feature vector of an 
%                           m by n image. u is the number of scales and
%                           v is the number of orientations in 'gaborArray'.
%
%
% Sample use:
% 
% img = imread('cameraman.tif');
% gaborArray = gaborFilterBank(5,8,39,39);  % Generates the Gabor filter bank
% featureVector = gaborFeatures(img,gaborArray,4,4);   % Extracts Gabor feature vector, 'featureVector', from the image, 'img'.
% 
% 
% 
%   Details can be found in:
%   
%   M. Haghighat, S. Zonouz, M. Abdel-Mottaleb, "CloudID: Trustworthy 
%   cloud-based and cross-enterprise biometric identification," 
%   Expert Systems with Applications, vol. 42, no. 21, pp. 7905-7916, 2015.
% 
% 
% 
% (C)	Mohammad Haghighat, University of Miami
%       haghighat@ieee.org
%       PLEASE CITE THE ABOVE PAPER IF YOU USE THIS CODE.


if (nargin ~= 4)        % Check correct number of arguments
    error('Please use the correct number of input arguments!')
end

if size(img,3) == 3     % Check if the input image is grayscale
    warning('The input RGB image is converted to grayscale!')
    img = rgb2gray(img);
end

img = double(img);

bands=size(img,3);

%% Filter the image using the Gabor filter bank

% Filter input image by each Gabor filter
%  gabor_all=[];
% persistent p;
% p=1;
gaborResult_all=[];
for k=1:bands
[u,v] = size(gaborArray);
gaborResult_k = cell(u,v);
for i = 1:u
    for j = 1:v
        gaborResult_k{i,j} = imfilter(img(:,:,k), gaborArray{i,j});  
%         gaborAbs = abs(gaborResult{i,j});
%         gaborAbs = gaborAbs(:);
%         gaborAbs = (gaborAbs-mean(gaborAbs))/std(gaborAbs,1);
%         gabor_all(:,:,p)= reshape(gaborAbs,m,n); 
% %         gabor_all(:,:,p)= gaborResult{i,j};   
%        p= p+1;
    end
end
        gaborResult_all=[gaborResult_all gaborResult_k];
end

%% Create feature vector
% Extract feature vector from input image
featureVector = [];
feature_all=[];
persistent q;
q=1;
for k=1:bands
    gaborResult=gaborResult_all(:,1+v*(k-1):v*k);
for i = 1:u
    for j = 1:v        
        gaborAbs = abs(gaborResult{i,j});
%         gaborAbs = gaborAbs(:);       
%         % Normalized to zero mean and unit variance. (if not applicable, please comment this line)
%         gaborAbs = (gaborAbs-mean(gaborAbs))/std(gaborAbs,1);
        gaborAbs=gaborAbs./max(gaborAbs(:));
        featureVector =  [featureVector; gaborAbs];
        feature_all(:,:,q)= reshape(gaborAbs,m,n);  
       q= q+1;
    end
end
end

%----------------------------original----------------------------------------------%

% %% Filter the image using the Gabor filter bank  
% 
% % Filter input image by each Gabor filter
% [u,v] = size(gaborArray);
% gaborResult = cell(u,v);
% for i = 1:u
%     for j = 1:v
%         gaborResult{i,j} = imfilter(img, gaborArray{i,j});
%     end
% end
% 
% 
% %% Create feature vector
% 
% % Extract feature vector from input image
% featureVector = [];
% for i = 1:u
%     for j = 1:v
%         
%         gaborAbs = abs(gaborResult{i,j});
%         gaborAbs = downsample(gaborAbs,d1);
%         gaborAbs = downsample(gaborAbs.',d2);
%         gaborAbs = gaborAbs(:);
%         
%         % Normalized to zero mean and unit variance. (if not applicable, please comment this line)
%         gaborAbs = (gaborAbs-mean(gaborAbs))/std(gaborAbs,1);
%         
%         featureVector =  [featureVector; gaborAbs];
%         
%     end
% end


%% Show filtered images (Please comment this section if not needed!)

% % Show real parts of Gabor-filtered images
% figure('NumberTitle','Off','Name','Real parts of Gabor filters');
% for i = 1:u
%     for j = 1:v        
%         subplot(u,v,(i-1)*v+j)    
%         imshow(real(gaborResult{i,j}),[]);
%     end
% end
% 
% % Show magnitudes of Gabor-filtered images
% figure('NumberTitle','Off','Name','Magnitudes of Gabor filters');
% for i = 1:u
%     for j = 1:v        
%         subplot(u,v,(i-1)*v+j)    
%         imshow(abs(gaborResult{i,j}),[]);
%     end
% end
