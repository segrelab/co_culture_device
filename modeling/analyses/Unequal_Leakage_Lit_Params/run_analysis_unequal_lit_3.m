%David Bernstein
%4/7/20
%simulate model for range of parameters

% Updated 4/19/20
% Save ratio of positive control growth at tspan(end) hrs to growth at 0 hours,
% normalize based on blank well

% Updated 3/30/22
% Narrow diffusion rate uncertainty based on lysine and isoleucine
% diffusion experiment

clear all
close all
clc

% Add path for model
addpath('../../model')
% Add path for data
addpath('../../../raw_data_and_plots')

%% Calculate Experimental Values of Growth Ratio

% Load Data
file = 'Figure 5 and 7 -- Syntrophic CoC.xlsx';
range = 'D28:BK412';
[num,~,~] = xlsread(file,range);

ind_blank = [32,34,36,38,40,41,43,45,47,49];
ind_initial = [31,33,35,37,39,42,44,46,48,50];
ind_PC = [53,54,55,56,57,58];
%ind_I = [4,14,...
%         6,16,26,...
%         18,28];
ind_I = [6,16,26];
%ind_K = [3,13,...
%         5,15,25,...
%         17,27];
ind_K = [5,15,25];

% Blank wells and initial amounts
m_blank = mean(mean(num(:,ind_blank)));
m_initial = mean(mean(num(:,ind_initial)));

% Fold Growth (for each time point)
Exp_r_K = (num(:,ind_K)-m_blank)./(m_initial-m_blank);
Exp_r_I = (num(:,ind_I)-m_blank)./(m_initial-m_blank);
Exp_r_PC = (num(:,ind_PC)-m_blank)./(m_initial-m_blank);

% Ratios @ time point
tspan = [0:0.25:48];

m_Exp_r_PC = mean(Exp_r_PC(length(tspan),:));
Exp_RT_K = m_Exp_r_PC./Exp_r_K(length(tspan),:);
Exp_RT_I = m_Exp_r_PC./Exp_r_I(length(tspan),:);

m_Exp_log_RT_K = mean([log10(Exp_RT_K)]);
s_Exp_log_RT_K = std([log10(Exp_RT_K)]);
m_Exp_log_RT_I = mean([log10(Exp_RT_I)]);
s_Exp_log_RT_I = std([log10(Exp_RT_I)]);

%% Simulate Co-culture experiments (Set Parameter Ranges)
% Add path for model
addpath('C:\Users\David\GitHub\co-culture_modeling\model')

n = 10000;
% Define Parameter Ranges
% Kinetics (Harcombe et al. cell systems. 2014; Gosset, Microb. Cell Fact. 2005)
r.vmaxG = 10.^unifrnd(1,1,n,1); % 10 %[mmol/(hr*g)]
r.vmaxA1 = 10.^unifrnd(1,1,n,1); % 10 %[mmol/(hr*g)]
r.vmaxA2 = 10.^unifrnd(1,1,n,1); % 10 %[mmol/(hr*g)]
r.kG = 10.^unifrnd(-2,-2,n,1); % 10e-2 %[mmol/L], 1.75e-3
r.kA1 = 10.^unifrnd(-2,-2,n,1); % 10e-2 %[mmol/L]
r.kA2 = 10.^unifrnd(-2,-2,n,1); % 10e-2 %[mmol/L]
% Biomass Stoichiometry
% grams e coli / grams glucose = 0.5/1 (BioNumbers 105318, Shiloach et al. Biotechnol Adv. 2005)
r.zG = 10.^unifrnd(log10(0.0901),log10(0.0901),n,1); % 0.5*180.156e-3=0.0901 %[g/mmol]
% cell / mol amino acid (Mee et al. PNAS. 2014)
% gram / cell = 500*10^-15 (cell biology by the numbers)
r.zA1 = 10.^unifrnd(log10(2.7374),log10(2.7374),n,1); %(5.4747e15)*(500e-15)*10^-3=2.7374 %[g/mmol] = [cell/mol]*[g/cell]*[mol/mmol] Lysine
r.zA2 = 10.^unifrnd(log10(4.0148),log10(4.0148),n,1); %(8.0295e15)*(500e-15)*10^-3=4.0148 %[g/mmol] = [cell/mol]*[g/cell]*[mol/mmol] Isoleucine
% Volume (From experiment)
r.v = 10.^unifrnd(log10(250*10^-6),log10(250*10^-6),n,1); %250*10^-6 [L]
% Diffusion (fit from chemical diffusion experiment: geometric mean of Phenol Red and Bromocresol Purple left and right side estimates)
r.d = 10.^unifrnd(log10(3.88e-5),log10(3.88e-5),n,1);

% Secretion Stoichiometry
y1 = -2;
y2 = 1;
r.yA1 = 10.^unifrnd(y1,y2,n,1);
r.yA2 = 10.^unifrnd(y1,y2,n,1);

% Initial conditions
initial_G = 5.55e-3; %[mmol] % Glucose [mmol]: 20% stock (20g/100mL) / MW * Conc. (2mL/100mL) * Vol (250uL)
initial_A1 = 1e-9; %[mmol] % AA1 [mmol]: aa/cell (1.1e8) * cell yield (10^9) / N_A (6.022e23) = 1.83e-4 mmol
initial_A2 = 1e-9; %[mmol] % AA2 [mmol]: aa/cell (7.5e7) * cell yield (10^9) / N_A (6.022e23) = 1.25e-4 mmol
initial_Bio1 = 1.23e-5; %[g]
initial_Bio2 = 1.23e-5; %[g] Biomass [g]:

MaxStep = 0.0025;

%% Simulate Co-culture experiments (Run model)
% % ================== COMMENT OUT TO RE-RUN ============= %
% % NOTE: This is a time consuming step so the results should be saved and
% % then loaded in the subsequent step for any re-analysis of data.
% 
% Sim_K = zeros(n,length(tspan));
% Sim_I = zeros(n,length(tspan));
% Sim_PC = zeros(n,length(tspan));
% for I1 = 1:n
%     I1
%     
%     p.vmaxG = r.vmaxG(I1);
%     p.vmaxA1 = r.vmaxA1(I1);
%     p.vmaxA2 = r.vmaxA2(I1);
%     p.kG = r.kG(I1);
%     p.kA1 = r.kA1(I1);
%     p.kA2 = r.kA2(I1);
%     p.zG = r.zG(I1);
%     p.zA1 = r.zA1(I1);
%     p.zA2 = r.zA2(I1);
%     p.v = r.v(I1);
%     p.d = r.d(I1);
%     p.yA1 = r.yA1(I1);
%     p.yA2 = r.yA2(I1);
%     % Run Model
%     % opposite wells (experiment)
%     % Initial Conditions
%     x0 = [initial_G,initial_G,... % Glucose [L,R]
%         initial_A1,initial_A1,... % AA1 [L,R]
%         initial_A2,initial_A2,... % AA2 [L,R]
%         initial_Bio1,0,... % Biomass1 [L,R]
%         0,initial_Bio2]; % Biomass2 [L,R]
%     options = odeset('NonNegative',[1:10],'MaxStep',MaxStep);
%     [to,xo] = ode23tb(@(t,x) f_co_culture_model(t,x,p), tspan, x0, options);
%     % same wells (positive control)
%     x0 = [initial_G,initial_G,... % Glucose [L,R]
%         initial_A1,initial_A1,... % AA1 [L,R]
%         initial_A2,initial_A2,... % AA2 [L,R]
%         initial_Bio1/2,initial_Bio1/2,... % Biomass1 [L,R]
%         initial_Bio2/2,initial_Bio2/2]; % Biomass2 [L,R]
%     options = odeset('NonNegative',[1:10],'MaxStep',MaxStep);
%     [ts,xs] = ode23tb(@(t,x) f_co_culture_model(t,x,p), tspan, x0, options);
%     
%     
%     Sim_K(I1,:) = xo(:,7);
%     Sim_I(I1,:) = xo(:,10);
%     Sim_PC(I1,:) = mean([xs(:,7)+xs(:,9),xs(:,8)+xs(:,10)],2);
% end
% save(strcat('r_',num2str(n),'_2'),'r')
% save(strcat('S_K_',num2str(n),'_2'),'Sim_K')
% save(strcat('S_I_',num2str(n),'_2'),'Sim_I')
% save(strcat('S_PC_',num2str(n),'_2'),'Sim_PC')
% % ================================= %

%% Evaluate and plot
% Load Data
load(strcat('r_',num2str(n),'_2'))
load(strcat('S_I_',num2str(n),'_2'))
load(strcat('S_K_',num2str(n),'_2'))
load(strcat('S_PC_',num2str(n),'_2'))

% Calculate Ratios
Sim_r_K = Sim_K(:,:)./Sim_K(:,1);
Sim_r_I = Sim_I(:,:)./Sim_I(:,1);
Sim_r_PC = Sim_PC(:,:)./Sim_PC(:,1);

% Check for numerical errors (all)
% rows = zeros(size(Sim_r_PC,1),1);
% for I = 1:size(Sim_r_PC,1)
%     for J = 1:size(Sim_r_PC,2)
%         if Sim_r_PC(I,J) < 1
%             I
%             rows(I) = 1;
%         end
%     end
% end
% rows = find(rows);


% Other ratios
Ratio_K = Sim_PC(:,end)./Sim_K(:,end);
Ratio_I = Sim_PC(:,end)./Sim_I(:,end);
Ratio_PC = Sim_PC(:,end)./Sim_PC(:,1);

% Check for numerical errors (end)
rows = find(Ratio_PC<1);
Ratio_K(rows) = [];
Ratio_I(rows) = [];
Ratio_PC(rows) = [];
n = n - length(rows);
r.yA1(rows) = [];
r.yA2(rows) = [];

% Rescale Ratio_PC between 0 and 1
Log_Ratio_PC = log10(Ratio_PC);
Ratio_PC_Scaled = (Log_Ratio_PC-min(Log_Ratio_PC))/(max(Log_Ratio_PC)-min(Log_Ratio_PC))/1.25+0.1;

lw = 1;

% addpath('C:\Users\David\Dropbox\My_Stuff\BU\Lab_Work\Segre_Lab\MATLAB\Misc_NonFBA\xkcd_rgb_v1.5\XKCD_RGB')
% addpath('C:\Users\David\Dropbox\My_Stuff\BU\Lab_Work\Segre_Lab\MATLAB\Misc_NonFBA\rgbmap_v2\')
% cmap = rgbmap('black','blue','light green','pale yellow',240);
% colormap(cmap)

colormap(parula)

figure(1)
set(gcf,'renderer','painters')
hold on
scatter(log10(r.yA1),log10(r.yA2),10,log10(Ratio_K),'filled','markerfacealpha',0.6)
hold on
x = log10(r.yA1);
y =log10(r.yA2);
z = log10(Ratio_K);
[xi,yi] = meshgrid(y1:0.1:y2, y1:0.1:y2);
zi = griddata(x,y,z,xi,yi,'linear');

thresh_K = 10*s_Exp_log_RT_K;
thresh_I = 10*s_Exp_log_RT_I;
lines = [m_Exp_log_RT_K+thresh_K,m_Exp_log_RT_K-thresh_K];
[C,h] = contour(xi,yi,zi,lines,'linewidth',lw);
%clabel(C,h)
cbh=colorbar;
cbh.Ticks = [0:0.2:1.6];
xlabel('log10(Lys Leakage [mmol/g])')
ylabel('log10(Ile Leakage [mmol/g])')
axis([y1 y2 y1 y2])

% figure(2)
% Ratio_I = 1./R_I(:,end);
% set(gcf,'renderer','painters')
% scatter(log10(r.d),log10(r.y),5,log10(Ratio_I),'filled','markerfacealpha',0.5)
% hold on
% colormap(cmap)
% x = log10(r.d);
% y =log10(r.y);
% z = log10(Ratio_I);
% [xi,yi] = meshgrid(-6.5:0.1:-3, -2:0.1:2);
% zi = griddata(x,y,z,xi,yi,'linear');
% lines = [0.5:2:5];
% [C,h] = contour(xi,yi,zi,lines,'linewidth',lw);
% %clabel(C,h)
% cbh=colorbar;
% cbh.Ticks = [0:0.5:5];
% xlabel('log10(Diffusion Rate [L/hr])')
% ylabel('log10(Leakage [mmol/g])')
% axis([-6.4 -3.1 -1.9 1.9])

%% Simulate Dynamics for select points
% leakage
leakages = [-1.5,0.9;-0.25,-0.25;0.9,-1.5;0.9,0.9;-1.5,-1.5];

blue1 = [0 0 0.7];
green1 = [0 0.7 0];
c_pc = [0.7 0.7 0.7];
fntsz = 12;

%MaxStep = 0.0025;

for I = 1:size(leakages,1)
    I
    figure(1)
    plot(leakages(I,1),leakages(I,2),'o','markersize',7,'markeredgecolor',[0.5 0 0],'markerfacecolor',[1 0.5 0.5])
    
%     figure(2)
%     plot(diffusion(I),leakage(I),'r.','markersize',15)
    
    p.yA1 = 10^leakages(I,1);
    p.yA2 = 10^leakages(I,2);
    
    p.d = r.d(1);
    p.vmaxG = r.vmaxG(1);
    p.vmaxA1 = r.vmaxA1(1);
    p.vmaxA2 = r.vmaxA2(1);
    p.kG = r.kG(1);
    p.kA1 = r.kA1(1);
    p.kA2 = r.kA2(1);
    p.zG = r.zG(1);
    p.zA1 = r.zA1(1);
    p.zA2 = r.zA2(1);
    p.v = r.v(1);

    % opposite wells
    % Initial Conditions
    x0 = [initial_G,initial_G,... % Glucose [L,R]
        initial_A1,initial_A1,... % AA1 [L,R]
        initial_A2,initial_A2,... % AA2 [L,R]
        initial_Bio1,0,... % Biomass1 [L,R]
        0,initial_Bio2]; % Biomass2 [L,R]
    % ODE
    options = odeset('NonNegative',[1:10],'MaxStep',MaxStep);
    [t,x] = ode23tb(@(t,x) f_co_culture_model(t,x,p), tspan, x0, options);
    % Record Results
    t3_1 = t;
    x3_1 = x;
    
    % same wells
    % Initial Conditions
    x0 = [initial_G,initial_G,... % Glucose [L,R]
        initial_A1,initial_A1,... % AA1 [L,R]
        initial_A2,initial_A2,... % AA2 [L,R]
        initial_Bio1/2,initial_Bio1/2,... % Biomass1 [L,R]
        initial_Bio2/2,initial_Bio2/2]; % Biomass2 [L,R]
    % ODE
    options = odeset('NonNegative',[1:10],'MaxStep',MaxStep);
    [t,x] = ode23tb(@(t,x) f_co_culture_model(t,x,p), tspan, x0, options);
    % Record Results
    t3_2 = t;
    x3_2 = x;
    
    figure(3)
    ind = (I-1)*2+1;
    subplot(size(leakages,1),2,ind)
    hold on
    plot(t,x3_1(:,7),'-','linewidth',1,'color',blue1)
    plot(t,x3_1(:,9),'-','linewidth',1,'color',green1)
    plot(t,x3_2(:,7)+x3_2(:,9),'linewidth',1,'color',c_pc,'linestyle',':')
    xlabel('Time [hrs]')
    ylabel('Biomass [g]')
    set(gca,'fontsize',fntsz)
    set(gca,'yscale','log')
    axis([0 tspan(end) 0 6*10^-4])
    set(gca,'xtick',0:24:tspan(end))
    set(gca,'ytick',[10^-4.5,10^-4,10^-3.5])
    set(gca,'yticklabels',{'10^{-4.5}','10^{-4}','10^{-3.5}'})
    set(gca, 'YMinorTick','off')
    
    subplot(size(leakages,1),2,ind+1)
    hold on
    plot(t,x3_1(:,8),'-','linewidth',1,'color',blue1)
    plot(t,x3_1(:,10),'-','linewidth',1,'color',green1)
    plot(t,x3_2(:,8)+x3_2(:,10),'linewidth',1,'color',c_pc,'linestyle',':')
    xlabel('Time [hrs]')
    set(gca,'fontsize',fntsz)
    set(gca,'yscale','log')
    axis([0 tspan(end) 0 6*10^-4])
    set(gca,'xtick',0:24:tspan(end))
    set(gca,'ytick',[10^-4.5,10^-4,10^-3.5])
    set(gca,'yticklabels',{'10^{-4.5}','10^{-4}','10^{-3.5}'})
    set(gca, 'YMinorTick','off')
    
    set(gcf,'Renderer', 'painters', 'Position', [0 0 500 800])
end


%% ABC Distribution Plot (End Point Ratio Data)
figure(4)
% Plot histogram of all points
Xedges = (y1:0.2:y2);
Yedges = (y1:0.2:y2);
histogram2(log10(r.yA1),log10(r.yA2),'XBinEdges',Xedges,'YBinEdges',Yedges,'normalization','pdf','facealpha',0.7,'facecolor','k');
% Plot histogram of all points for which ratio is within threshold
hold on
pass = zeros(n,1);
for I = 1:n
    if abs(log10(Ratio_K(I))-m_Exp_log_RT_K)<thresh_K
        if abs(log10(Ratio_I(I))-m_Exp_log_RT_I)<thresh_I
            pass(I) = 1;
        end
    end
end
x1 = log10(r.yA1(pass==1));
x2 = log10(r.yA2(pass==1));
histogram2(x1,x2,'XBinEdges',Xedges,'YBinEdges',Yedges,'normalization','pdf','facealpha',0.7,'facecolor',[0.9,0.75,0])

xlabel('log10(Lys Leakage [mmol/g])')
ylabel('log10(Ile Leakage [mmol/g])')
zlabel('Probability Density Function')

figure(5)
% Plot histogram of all points
histogram(log10(r.yA1),'BinEdges',Xedges,'normalization','pdf','facealpha',0.5,'facecolor','k');
% Plot histogram of all points for which ratio is within threshold
hold on
histogram(log10(r.yA1(pass==1)),'BinEdges',Xedges,'normalization','pdf','facealpha',0.5,'facecolor',[0.9,0.75,0])
xlabel('log10(Lys Leakage [mmol/g])')
ylabel('Probability Density Function')

figure(6)
% Plot histogram of all points
histogram(log10(r.yA1),'BinEdges',Yedges,'normalization','pdf','facealpha',0.5,'facecolor','k','Orientation','horizontal');
% Plot histogram of all points for which ratio is within threshold
hold on
histogram(log10(r.yA2(pass==1)),'BinEdges',Yedges,'normalization','pdf','facealpha',0.5,'facecolor',[0.9,0.75,0],'Orientation','horizontal')
ylabel('log10(Ile Leakage [mmol/g])')
xlabel('Probability Density Function')
axis([0 1.6 -2 1])

%% Format Plots
fs = 10;

figure(1)
set(gcf,'renderer','painters','position',[0 0 500 400])

figure(3)
set(gcf,'renderer','painters','position',[0 0 380 700])

figure(5)
set(gcf,'renderer','painters','position',[0 0 500 230])
set(gca,'fontsize',fs)

figure(6)
set(gcf,'renderer','painters','position',[0 0 250 400])
set(gca,'fontsize',fs)
