function [Forest] = Improved_iForest(Data, NumTree, NumSub, NumDim, K, rseed, Exlevel)
    %======================================================================
    % Function Rel_IForest_K: build isolation forest
    %
    % Inputs:
    %     Data: n x d matrix; n: # of instance; d: dimension;
    %     NumTree: # of isolation trees;
    %     NumSub: # of sub-sample;
    %     NumDim: # of sub-dimension;
    %     K: Minimum number of instances in the leaf node.
    %     rseed: random seed;
    %     Exlevel: Exention level to be used in the creating splitting critera.
    %
    % Outputs:
    %     Forest: structure; an isolation forest model
    %     Forest.Trees: a half space forest model;
    %     Forest.NumTree: NumTree;
    %     Forest.NumSub: NumSub;
    %     Forest.NumDim: NumDim;
    %     Forest.K: K;
    %     Forest.HeightLimit: height limitation;
    %     Forest.ElapseTime: elapsed time;
    %     Forest.rseed: rseed;
    %====================================================================== 
    
    [NumInst, DimInst] = size(Data);
    
    Forest.Trees = cell(NumTree, 1);
    Forest.NumTree = NumTree;
    Forest.NumSub = NumSub;
    Forest.NumDim = NumDim;
    Forest.NumInst = NumInst;
    % Forest.HeightLimit = NumSub;
    Forest.HeightLimit = ceil(log2(NumSub));
    Forest.K = K;
    Forest.rseed = rseed;
    Forest.Exlevel = Exlevel; % Exention level to be used in the creating splitting critera.
%     rand('state', rseed);
    
    % data transformation
    Forest.DataMass = zeros(NumInst, NumTree);
    Forest.ParentMass = zeros(NumInst, NumTree);
    Forest.DataPL = zeros(NumInst, NumTree);

    % parameters for function IsolationTree
    Paras.HeightLimit = Forest.HeightLimit;
    Paras.NumDim = NumDim;
    Paras.K = K;
    
    et = cputime;

    for i = 1:NumTree

        if NumSub < NumInst % randomly selected sub-samples
            [~, SubRand] = sort(rand(1, NumInst));
            IndexSub = SubRand(1:NumSub);
        else
            IndexSub = 1:NumInst;
        end
        
        if ((NumDim < DimInst) && (NumDim > 0)) % randomly selected sub-dimensions
            [~, DimRand] = sort(rand(1, DimInst));
            IndexDim = DimRand(1:NumDim);
        else
            IndexDim = 1:DimInst;
            Paras.NumDim = DimInst;
        end

        Paras.IndexDim = IndexDim;

        DataIndex = 1:1:NumInst;
        
        [Forest.Trees{i}, DataMass, ParentMass, DataPL] = Improved_iTree(Data, IndexSub, DataIndex, 0, Paras,...
             (NumInst .* zeros(NumInst,1)), (NumInst .* zeros(NumInst,1)), (NumInst .* zeros(NumInst,1)) ,Forest.Exlevel); 
        
        Forest.DataMass(:,i) = DataMass;
        Forest.ParentMass(:,i) = ParentMass;
        Forest.DataPL(:,i) = DataPL;
        
    end

    Forest.ElapseTime = cputime - et;

end
