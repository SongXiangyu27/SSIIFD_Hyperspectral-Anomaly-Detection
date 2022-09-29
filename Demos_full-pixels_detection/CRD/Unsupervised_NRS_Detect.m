function result = Unsupervised_NRS_Detect(Data, win_out, win_in, lambda)
%Unsupervised NRS anomaly detector
%   Unsupervised NRS performs the RX anomaly detector
%
% Usage
%   [result] = Unsupervised_NRS_Detect(Data, window, lambda)
% Inputs
%   Data - 3D data matrix (num_row x num_col x num_dim)
%   window - spatial size window (e.g., 3, 5, 7, 9,...)
%   lambda - regularization parameter
% Outputs
%   result - Detector output (num_row x num_col)
%  

[a b c] = size(Data);        
result = zeros(a, b);
t = fix(win_out/2);
t1 = fix(win_in/2);
M = win_out^2;

% padding avoid edges
DataTest = zeros(3*a, 3*b, c);
DataTest(a+1:2*a, b+1:2*b, :) = Data;
DataTest(a+1:2*a, 1:b, :) = Data(:, b:-1:1, :);
DataTest(a+1:2*a, 2*b+1:3*b, :) = Data(:, b:-1:1, :);
DataTest(1:a, :, :) = DataTest(2*a:-1:(a+1), :, :);
DataTest(2*a+1:3*a, :, :) = DataTest(2*a:-1:(a+1), :, :);

for i = 1+b: 2*b 
    for j = 1+a: 2*a
        block = DataTest(j-t: j+t, i-t: i+t, :);
        y = squeeze(DataTest(j, i, :)).';
        block(t-t1+1:t+t1+1, t-t1+1:t+t1+1, :) = NaN;
        block = reshape(block, M, c);
        block(isnan(block(:, 1)), :) = [];
        H = block';  % num_dim x num_sam
        Psi = PCA_Train(H, 50);
        Sigma = (Psi * Psi');
        
        % implementation two
        z = H - repmat(y',1,size(H,2));
        C = z'*Sigma*z;
        tmp = diag(C);
        R1 = diag(lambda*(tmp./max(tmp(:))));   % optimal lambda = 0.1
        % R2 = diag(lambda*ones(size(H,2),1));    % optimal lambda = 0.1
        C = C + R1; % REGULARIZATION
        invC = pinv(C);
        weights = sum(invC)'/sum(sum(invC));  % K *1
        
        
        y_hat = (H*weights(:))';  % 1 x num_dim
        result(j-a, i-b) = norm(y - y_hat, 1);
    end
end

