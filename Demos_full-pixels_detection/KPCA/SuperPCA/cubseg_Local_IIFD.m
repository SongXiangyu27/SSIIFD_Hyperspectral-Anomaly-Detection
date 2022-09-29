function labels = cubseg_Local_IIFD(data_seg, num_Pixel, num_rows, num_cols)

[MN,B]=size(data_seg);
M = num_rows;   %sqrt(MN);       
N = num_cols;   %sqrt(MN);       
Y_scale=scaleForSVM(data_seg); % 
% Y=reshape(Y_scale, M, N, B);
% grey_img = im2uint8(mat2gray(Y(:, :, 30)));

p = 1;
[Y_pca] = pca(Y_scale, p);
img = im2uint8(mat2gray(reshape(Y_pca', M, N, p)));

% img = data_seg;
% sigma=0.05;
K=num_Pixel;
labels = mex_ers(double(img), K);

[height, width] = size(img);
[bmap] = seg2bmap(labels,width,height);
bmapOnImg = img;
idx = find(bmap>0);
timg = img;
timg(idx) = 255;
bmapOnImg(:,:,2) = timg;
bmapOnImg(:,:,1) = img;
bmapOnImg(:,:,3) = img;
figure;
imshow(bmapOnImg,[]);
title('superpixel boundary map');
