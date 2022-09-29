function [AUC,AUCnor] = cal_AUC(det_map,GT,mode_tau,mode_eq)
%  det_map; the detection result N*k, N is the number of pixels in the detection result, 
%          k is the number of detection results with different detector
%  GT: the ground truth vector N*1
%  mode_tau: if mode_tau==1,the tau is sample value; else the tau is uniform step 0.01
%  mode_eq: if mode_eq==1. the equation (7) in the paper is used;the equation (9)is used;
%%%%% 
num_map = size(det_map,2);
for i = 1:num_map
    det_map(:,i) = (det_map(:,i) - min(det_map(:,i))) /(max(det_map(:,i))-min(det_map(:,i)));
end

% PD and PF
if mode_tau == 1
    tau = sort(det_map,'descend');
else
    for k = 1:num_map
        tau(:,k) = [0:0.01:1]';
    end
    tau = sort(tau,'descend');
end
for k = 1:num_map
    for i = 1: length(tau)
        AD_bw =det_map(:,k);
        if mode_eq==1
           AD_bw(det_map(:,k)>=tau(i,k))=1;
           AD_bw(det_map(:,k)<tau(i,k))=0;
        else
           AD_bw(det_map(:,k)>tau(i,k))=1;
           AD_bw(det_map(:,k)<=tau(i,k))=0;
        end
        [PD(i,k),PF(i,k)] = cal_pdpf(AD_bw,GT);
    end
end

% Compute a1, a0, b1 and b0 and normalized PD PF
 a1 = min(PD(1,:));
 a0 = max(max(PD)); 
 b1=min(PF(1,:));
 b0 = max(max(PF));
 

AUCnor.a1 = a1;
AUCnor.a0 = a0;
AUCnor.b1 = b1;
AUCnor.b0 = b0;


% AUC(PF, PD) and Nor_AUC(PF, PD)
for i = 1:num_map
    AUC0 = trapz(PF(:,i),PD(:,i));    
    AUC.PFPD(i,1) = round(AUC0,6);
    AUCnor.PFPD(i,1) =round((AUC0-a1)/((a0-a1)*(b0-b1)),6);
end

% AUC(PD, tau) and Nor_AUC(PD, tau)
for i = 1:num_map
    AUC0 = -trapz(tau(:,i),PD(:,i));
    AUC.tauPD(i,1) = round(AUC0,6);
    AUCnor.tauPD(i,1) = round((AUC0-a1)/(a0-a1),6);
end

% AUC(PF£¬tau£©and Nor_AUC(PF£¬tau£©
for i = 1:num_map
    AUC0 = abs(trapz(tau(:,i),PF(:,i)));
    AUC.tauPF(i,1) =round(AUC0,6);
    AUCnor.tauPF(i,1) =round((AUC0-b1)/(b0-b1),6);
end


