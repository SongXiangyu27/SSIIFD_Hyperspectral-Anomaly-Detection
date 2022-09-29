function plot_3DROC(det_map,GT,detec_label,mode_eq)
% det_map: the detection result N*k, N is the number of pixels in the detection result, 
% GT: Ground truth
% detec_label: the name of detector for legend
% mode_eq: if mode_eq==1. the equation (7) in the paper is used;the equation (9)is used;

num_map = size(det_map,2);
for i = 1:num_map
    det_map(:,i) = (det_map(:,i) - min(det_map(:,i))) /(max(det_map(:,i))-min(det_map(:,i)));
end

%PD and PF based on uniform step and sample value
for k = 1:num_map
    tau1(:,k) = [0:0.01:1]';
end
 tau1 = sort(tau1,'descend');

for k = 1:num_map
    for i = 1: length(tau1)
        map =det_map(:,k);
        if mode_eq==1
           map(det_map(:,k)>=tau1(i,k))=1;
           map(det_map(:,k)<tau1(i,k))=0;
        else
           map(det_map(:,k)>tau1(i,k))=1;
           map(det_map(:,k)<=tau1(i,k))=0;
        end
        [PD1(i,k),PF1(i,k)] = cal_pdpf(map,GT);
    end
end
map = [];

tau2 = sort(det_map,'descend');
 
 for k = 1:num_map
    for i = 1: length(tau2)
        map =det_map(:,k);
        if mode_eq==1
           map(det_map(:,k)>=tau2(i,k))=1;
           map(det_map(:,k)<tau2(i,k))=0;
        else
           map(det_map(:,k)>tau2(i,k))=1;
           map(det_map(:,k)<=tau2(i,k))=0;
        end
        [PD2(i,k),PF2(i,k)] = cal_pdpf(map,GT);
    end
 end
 
% 
 a11 = min(PD1(1,:));
 a10 = max(max(PD1)); 
 b11=min(PF1(1,:));
 b10 = max(max(PF1));
 
 a21 = min(PD2(1,:));
 a20 = max(max(PD2)); 
 b21=min(PF2(1,:));
 b20 = max(max(PF2)); 
PD1nor = (PD1-a11)/(a10-a11);
PF1nor = (PF1-b11)/(b10-b11);
PD2nor = (PD2-a21)/(a20-a21);
PF2nor = (PF2-b21)/(b20-b21);

%% 2D ROC

% show ROC (PF, PD)
figure,plot(PF1,PD1,'LineWidth',2)
hold on 
plot(PF2,PD2,'LineWidth',2,'linestyle','--')
axis([0,1,0,1])

set(gca,'XTick',(0:0.2:1),'fontsize',16)
set(gca,'YTick',(0:0.2:1),'fontsize',16)
set(gca, 'XDir','reverse')
xlabel('P_F','fontsize',18) ; ylabel('P_D','fontsize',18) 
grid on

for i = 1:num_map
    name1(i) =strcat(detec_label(i),',','{\Delta}','=0.01'); 
    name2(i) = strcat(detec_label(i),',','sample values');
end
legend([name1,name2])
legend boxoff

% show normalized ROC (PFnor, PDnor)

figure,plot(PF1nor,PD1nor,'LineWidth',2)
hold on 
plot(PF2nor,PD2nor,'LineWidth',2,'linestyle','--')
axis([0,1,0,1])

set(gca,'XTick',(0:0.2:1),'fontsize',16)
set(gca,'YTick',(0:0.2:1),'fontsize',16)
set(gca, 'XDir','reverse')

xlabel('P_F','fontsize',18) ; ylabel('P_D','fontsize',18) 
grid on
legend([name1,name2])
hold off
legend boxoff


% show ROC (PD, Tau) based on uniform step and sample value
figure,plot(tau1,PD1,'LineWidth',2)
hold on 
plot(tau2,PD2,'LineWidth',2,'linestyle','--')
axis([0,1,0,1])

set(gca,'XTick',(0:0.2:1),'fontsize',16)
set(gca,'YTick',(0:0.2:1),'fontsize',16)
xlabel('\tau','fontsize',18) ; ylabel('P_D','fontsize',18) 
grid on

hold off
legend([name1,name2])
legend boxoff

% show normalized ROC (PDnor, Tau)

figure,plot(tau1,PD1nor,'LineWidth',2)
hold on 
plot(tau2,PD2nor,'LineWidth',2,'linestyle','--')
axis([0,1,0,1])

set(gca,'XTick',(0:0.2:1),'fontsize',16)
set(gca,'YTick',(0:0.2:1),'fontsize',16)
xlabel('\tau','fontsize',18) ; ylabel('P_D','fontsize',18) 
grid on
legend([name1,name2])
hold off
legend boxoff

 

% show ROC (PF, Tau) based on uniform step and sample value
figure,plot(tau1,PF1,'LineWidth',2)
hold on 
plot(tau2,PF2,'LineWidth',2,'linestyle','--')
axis([0,1,0,1])

set(gca,'XTick',(0:0.2:1),'fontsize',16)
set(gca,'YTick',(0:0.2:1),'fontsize',16)
xlabel('\tau','fontsize',18) ; ylabel('P_F','fontsize',18) 
grid on

hold off
legend([name1,name2])
legend boxoff

% show normalized ROC (PFnor, Tau)

figure,plot(tau1,PF1nor,'LineWidth',2)
hold on 
plot(tau2,PF2nor,'LineWidth',2,'linestyle','--')
axis([0,1,0,1])

set(gca,'XTick',(0:0.2:1),'fontsize',16)
set(gca,'YTick',(0:0.2:1),'fontsize',16)
xlabel('\tau','fontsize',18) ; ylabel('P_F','fontsize',18) 
grid on

legend([name1,name2])
hold off
legend boxoff

%% 3D ROC 

% 3D ROC

figure,plot3(PF1,tau1,PD1,'LineWidth',2)
hold on 
plot3(PF2,tau2,PD2,'LineWidth',2,'linestyle','--')
hold off

axis([0, 1, 0, 1, 0, 1])
set(gca,'XTick',(0:0.2:1),'fontsize',16)
set(gca,'YTick',(0:0.2:1),'fontsize',16)
set(gca,'ZTick',(0:0.2:1),'fontsize',16)
xlabel('P_F','fontsize',18) ; ylabel('\tau','fontsize',18); zlabel('P_D','fontsize',18) 
grid on
ax = gca;
ax.BoxStyle = 'full';
box on
legend([name1,name2],'fontsize',12)
legend boxoff

% Normalized 3D ROC 

figure,plot3(PF1nor,tau1,PD1nor,'LineWidth',2)
hold on 
plot3(PF2nor,tau2,PD2nor,'LineWidth',2,'linestyle','--')
hold off
legend([name1,name2],'fontsize',12)

axis([0, 1, 0, 1, 0, 1])
set(gca,'XTick',(0:0.2:1),'fontsize',16)
set(gca,'YTick',(0:0.2:1),'fontsize',16)
set(gca,'ZTick',(0:0.2:1),'fontsize',16)
xlabel('P_F','fontsize',18) ; ylabel('\tau','fontsize',18); zlabel('P_D','fontsize',18) 
grid on
ax = gca;
ax.BoxStyle = 'full';
box on
legend boxoff




