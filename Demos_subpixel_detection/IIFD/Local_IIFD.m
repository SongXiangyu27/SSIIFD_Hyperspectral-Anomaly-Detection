function [PMass_spatial,QMass_spatial,score,AUC, time] =  Local_IIFD(Data_seg, data, map, NumTree, num_rows, num_cols)

    %--------------------------------------------------------------------------
    % Run an Anomaly Detection algorithm on the given data set.
    % (C) Xiangyu Song, University of Chinese Academy of Sciences & Sunil Aryal, Monash University
    %
    % Input:
    %     Data_seg  : .mat file with num_rows*num_cols matrix; image for segmentation using ERS;
    %     data      : .mat file with n x d data matrix (n: # of instance d: dimension)
    %     map       : n x 1 matrix, Ground truth (1-Anomaly, 0-Normal);
    %     NumTree   : # Number of iTrees;
    %
    % Output:
    %     AUC                 : the AUC result (NumItx x 6 AUCs at different iterations with 6 sample size values);
    %     time.buildTime      : time elapsed in building template and sampling;
    %     time.estimationTime : time elapsed in estimating mass
    %     time.evaluationTime : evaluation time  
    %     PMass_spatial       : data mass of the immediate parent of point x_i 
    %     QMass_spatial       : data mass of the lesf node in which point x_i falls  
    %--------------------------------------------------------------------------
    data_seg = Data_seg;
 
    % get the information about data
    [NumInst, NumDim] = size(data);
    
%     NumTree = NumTree;   %512; %100  256  512  1024
%     NumSub = ceil(NumInst * 0.03); %1000
    num_Pixel = 3; 
    
    RanSeed = cputime;
    rand('state', RanSeed);

    % start the clock
    sTime = cputime;

 %------------------------------segment data with superpixel------------------------------------------
    labels = cubseg_Local_IIFD(data_seg, num_Pixel, num_rows, num_cols);
%     Results_segment = seg_im_class_for_RMiF_3D(data3D, labels);
    Results_segment = seg_im_class_IIF(data, labels); 
    Num=size(Results_segment.Y,2);
    score_raw = zeros(NumInst, 1);
    for i=1:Num
        [NumInst, NumD] = size(Results_segment.Y{1,i});
        NumSub = ceil(NumInst * 0.2); %0.03
        k_local = ceil(NumSub * 0.05); %ceil(NumInst * 0.025); 
        Exlevel = ceil(NumD * 0.5);  %7; 
        Forest = Improved_iForest(Results_segment.Y{1,i}, NumTree, NumSub, NumDim, k_local, RanSeed, Exlevel); % use the extended method
        % DataMass and Parent mass
        Q_Mass = Forest.DataMass;
        P_Mass = Forest.ParentMass;
        ratio = P_Mass ./ Q_Mass;
        subscore = mean(ratio, 2);

        score_raw(Results_segment.index{1,i},:) = subscore;
        P_Massraw(Results_segment.index{1,i},:) = P_Mass;
        Q_Massraw(Results_segment.index{1,i},:) = Q_Mass;
        
    end
    score = score_raw;
    
    PMass_spatial = P_Massraw;
    QMass_spatial = Q_Massraw;
    
    % get the template building time
    time.buildTime = cputime - sTime;
    
    % restart the clock
    sTime = cputime;
    
    % get the mass estimation time
    time.estimationTime = cputime - sTime;
    
    sTime = cputime;
    
    auc = Measure_AUC(score, map);
    
    time.evaluationTime = cputime - sTime;
    
    AUC = auc;
    RTIME = time;
   
end
