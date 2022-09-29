function [Tree, DataMass, ParentMass, DataPL] = Improved_iTree(Data, CurtIndex, DataIndex, CurtHeight, Paras, DMass, PMass, DPL, Exlevel)
    %======================================================================
    % Function Rel_ITree_K: build an isolation tree
    % 
    % Inputs:
    %     Data: n x d matrix; n: # of instance; d: dimension (the whole data set);
    %     CurtIndex: vector; indices of the current sample instances;
    %     DataIndex: vector; indices of the data instances;
    %     CurtHeight: current height;
    %     Paras.IndexDim: sub-dimension index;
    %     Paras.NumDim:  # of sub-dimension;
    %     Paras.HeightLimit: maximun height limitation;
    %     Paras.K: # minimum points in the leaf node;
    %     DMass: n x 1 - Data mass of the leaf node for each instance
    %     PMass: n x 1 - Mass of the parent node for each instance
    %     DPL: n x 1 - Pathlength of each instance
    %     Exlevel : int  Specifies degree of freedom in choosing the hyperplanes for dividing up data. 
    % 
    % Outputs:
    %     Tree: an isolation tree model
    %     Tree.NodeStatus: 1: inNode, 0: exNode;
    %     Tree.SplitAttribute: splitting attribute;
    %     Tree.SplitPoint: splitting point;
    %     Tree.LeftChild: left child tree which is also a half space tree model;
    %     Tree.RightChild: right child tree which is also a half space tree model;
    %     Tree.DataMass: node data mass;
    %     Tree.Size: node subsample mass;
    %     Tree.Height: node height;
    %     Tree.DataIndex: vector of data indices falling in that node
    %======================================================================
    
    DataMass = DMass;
    ParentMass = PMass;
    DataPL = DPL;
    
    Tree.Height = CurtHeight;
    Tree.DataIndex = DataIndex;
    Tree.DataMass = size(DataIndex,2); % Tree.DataMass is the # of instance (the whole data set)
    Tree.Size = size(CurtIndex,2);
    NumInst = length(CurtIndex);
    NumInst_Data = length(DataIndex);
    
    %DataMass(DataIndex) = Tree.DataMass;
    DataMass(DataIndex) = Tree.Size;
    DataPL(DataIndex) = Tree.Height;
    % avoid acept 0 in DataMass
    dmidx = DataMass == 0;
    DataMass(dmidx) = 1; 
%     dmidx = find(DataMass == 0);
%     DataMass(dmidx) = 1; 
    exlevel = Exlevel;
    
    % check if terminating condition reached
    if CurtHeight >= Paras.HeightLimit || Tree.Size <= Paras.K
        Tree.NodeStatus = 0;
        Tree.SplitAttribute = [];
        Tree.SplitPoint = [];
        Tree.LeftChild = [];
        Tree.RightChild = [];
        Tree.Nt = [];
        Tree.Pt = [];
        if (Tree.Size > 1)
            c = 2 * (log(Tree.Size - 1) + 0.5772156649) - 2 * (Tree.Size - 1) / Tree.Size;
            DataPL(DataIndex) = DataPL(DataIndex) + c;
        end
        return;
    else
        % get the samples
        Samples = Data(CurtIndex,Paras.IndexDim);
        % check for the candidate attribute list for split
        Sammax = max(Samples);
        Sammin = min(Samples);
        AttRange = Sammax - Sammin;
        AttList = find(AttRange > 0);
        SizeAttList = size(AttList,2);
        % return if no candiate attribute is found
        if (SizeAttList == 0)
            Tree.NodeStatus = 0;
            Tree.SplitAttribute = [];
            Tree.SplitPoint = [];
            Tree.LeftChild = [];
            Tree.RightChild = [];
            if (Tree.Size > 1)
                c = 2 * (log(Tree.Size - 1) + 0.5772156649) - 2 * (Tree.Size - 1) / Tree.Size;
                DataPL(DataIndex) = DataPL(DataIndex) + c;
            end
            return;
        else 
            Tree.NodeStatus = 1;
            
            % update parent mass
            % ParentMass(DataIndex) = Tree.DataMass;
             ParentMass(DataIndex) = Tree.Size;
            [~,d] = size(Data);
          
            best_sep = zeros(1,d);
            p = zeros(1,d);
            for i = 1 : d
                 [best_sep(i), p(i)] = find_nsplits_value_IIF( Samples(:,i) );
            end
            [~, index_bands] = sort(best_sep);

            index_zero = index_bands(1: d - exlevel);  %[1: d - exlevel, d - 2: d]
           
            Tree.Nt = normrnd(0,1,1,d);% A random normal vector£¨1*dim£© picked form a uniform n-sphere.
            Tree.Nt = Tree.Nt/norm(Tree.Nt); % A random normal vector picked form a uniform n-sphere. 
                                             % Note that in order to pick uniformly from n-sphere, we need to pick a random normal for each component of this vector.
            Tree.Nt(index_zero) = 0;
            Tree.Pt = p;  

            % select a split point
%            CurtSamples = Data(CurtIndex, Tree.SplitAttribute);
%            CurtData = Data(DataIndex, Tree.SplitAttribute);
            CurtSamples = Data(CurtIndex, :);
            CurtData = Data(DataIndex, :);
%           Tree.SplitPoint = min(CurtSamples) + (max(CurtSamples) - min(CurtSamples)) * rand(1);
            Tree.SplitAttribute = repmat(Tree.Nt, NumInst, 1);
            Tree.SplitAttribute_Data = repmat(Tree.Nt, NumInst_Data, 1);
            CurtSam = CurtSamples - Tree.Pt;
            CurtDat = CurtData -  Tree.Pt;
            SamW = dot(CurtSam' , Tree.SplitAttribute');
            DatW = dot(CurtDat' , Tree.SplitAttribute_Data');
            Tree.SplitAttribute = SamW';
            Tree.SplitAttribute_Data = DatW';
            Tree.SplitPoint = 0;
            % sample indices for the left child and the right child
            LeftCurtIndex = CurtIndex(Tree.SplitAttribute <= Tree.SplitPoint);
            RightCurtIndex = setdiff(CurtIndex, LeftCurtIndex);
            % Data indices
            LeftDataIndex = DataIndex(Tree.SplitAttribute_Data' <= Tree.SplitPoint);
            RightDataIndex = setdiff(DataIndex, LeftDataIndex);
            % built right and left trees
            [Tree.LeftChild, DataMass, ParentMass, DataPL] = Improved_iTree(Data, LeftCurtIndex, LeftDataIndex, CurtHeight + 1, Paras, DataMass, ParentMass, DataPL, Exlevel);
            [Tree.RightChild, DataMass, ParentMass, DataPL] = Improved_iTree(Data, RightCurtIndex, RightDataIndex, CurtHeight + 1, Paras, DataMass, ParentMass, DataPL, Exlevel);
        end
    end
end
