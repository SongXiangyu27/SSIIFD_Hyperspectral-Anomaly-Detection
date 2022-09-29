function [PMass_spectral,QMass_spectral,score,AUC, time] =  IIFD_Gabor(Data, map, NumTree, K)

    %--------------------------------------------------------------------------
    % Run an Anomaly Detection algorithm on the given data set.
    % (C) Xiangyu Song, University of Chinese Academy of Sciences & Sunil Aryal, Monash University
    %
    % Input:
    %     Data      : .mat file with n x d data matrix (n: # of instance d: dimension)
    %     map       : n x 1 matrix, Ground truth (1-Anomaly, 0-Normal);
    %     NumTree   : # Number of iTrees;
    %     K         : Minimum number of instances in the leaf node.
    %
    % Output:
    %     auc                 : the AUC result (NumItx x 6 AUCs at different iterations with 6 sample size values);
    %     time.buildTime      : time elapsed in building template and sampling;
    %     time.estimationTime : time elapsed in estimating mass
    %     time.evaluationTime : evaluation time  
    %     PMass_spatial       : data mass of the immediate parent of point x_i 
    %     QMass_spatial       : data mass of the lesf node in which point x_i falls  
    %--------------------------------------------------------------------------
    data_mxn = Data;
    
    % get the information about data
    [NumInst, NumDim] = size(data_mxn);
    
    NumSub = ceil(NumInst * 0.03); %1000
    Exlevel = ceil(0.5*NumDim);  %50;   %ceil(N/2);
    RanSeed = cputime;
    rand('state', RanSeed);
    
    % start the clock
    sTime = cputime;
    
    % buuild the template
    Forest = Improved_iForest(data_mxn, NumTree, NumSub, NumDim, K, RanSeed, Exlevel); % use the extended method
    
    % get the template building time
    time.buildTime = cputime - sTime;
    
    % restart the clock
    sTime = cputime;
    
    % DataMass and Parent mass
    Q_Mass = Forest.DataMass;
    P_Mass = Forest.ParentMass;
    ratio = P_Mass ./ Q_Mass;
    score_raw = mean(ratio, 2);
    
    value = 2 * (log(NumInst-1)+0.5772156649) - 2*(NumInst-1) / NumInst;
    score = 1 - 2.^(-score_raw ./ value);   % 2^(-scores(i) / length(forest) / adjustment(m));
    %     score = score';
    
    PMass_spectral = P_Mass;
    QMass_spectral = Q_Mass;
    
    % get the mass estimation time
    time.estimationTime = cputime - sTime;
    sTime = cputime;
    auc = Measure_AUC(score, map);
    time.evaluationTime = cputime - sTime;
    AUC = auc;
    RTIME = time;
   
end
