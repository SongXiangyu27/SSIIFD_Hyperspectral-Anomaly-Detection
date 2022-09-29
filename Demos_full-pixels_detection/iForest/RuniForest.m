function score = RuniForest(Data, NumTree, NumSub, NumDim)
%
% Xiangyu Song
% 
% function RuniForest: run iForest anomaly detection 
% use the path length to calculate the anomaly score
%
% Input:
%     Data: n x d matrix; n: # of instance; d: dimension;
%     NumTree: # of isolation trees;
%     NumSub: # of sub-sample;
%     NumDim: # of sub-dimension;
%     rseed: random seed;
%
% Output:
%
%     score: anomaly score i.e. detection result
%

function value=adjustment(n)
%The average path length of an unsuccessful BST search
    if n<=1
        value = 0;
    else
        value = 2 * (log(n-1)+0.5772156649) - 2*(n-1) / n;
    end
end

[Num_pixels,~] = size(Data);
rseed = sum(100 * clock);
Forest = IsolationForest(Data, NumTree, NumSub, NumDim, rseed);
[Mass, ~] = IsolationEstimation(Data, Forest);
Emass =  mean(Mass, 2);
c = adjustment(Num_pixels);
score = 2.^(-(Emass/c));  % 

end
