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
ind_I = [4,14,...
         6,16,26,...
         18,28];
ind_K = [3,13,...
         5,15,25,...
         17,27];

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

m_Exp_log_RT_I_K = mean([log10(Exp_RT_K),log10(Exp_RT_I)]);
s_Exp_log_RT_I_K = std([log10(Exp_RT_K),log10(Exp_RT_I)]);


%% Simulate Co-culture experiments (Set Parameter Ranges)

n = 10000;
dist = 0.5;
% Define Parameter Ranges
% Kinetics (Harcombe et al. cell systems. 2014; Gosset, Microb. Cell Fact. 2005)
r.vmaxG = 10.^unifrnd(1-dist,1+dist,n,1); % 10 %[mmol/(hr*g)]
r.vmaxA1 = 10.^unifrnd(1-dist,1+dist,n,1); % 10 %[mmol/(hr*g)]
r.vmaxA2 = 10.^unifrnd(1-dist,1+dist,n,1); % 10 %[mmol/(hr*g)]
r.kG = 10.^unifrnd(-2-dist,-2+dist,n,1); % 10e-2 %[mmol/L], 1.75e-3
r.kA1 = 10.^unifrnd(-2-dist,-2+dist,n,1); % 10e-2 %[mmol/L]
r.kA2 = 10.^unifrnd(-2-dist,-2+dist,n,1); % 10e-2 %[mmol/L]
% Biomass Stoichiometry
% grams e coli / grams glucose = 0.5/1 (BioNumbers 105318, Shiloach et al. Biotechnol Adv. 2005)
r.zG = 10.^unifrnd(log10(0.0901)-dist,log10(0.0901)+dist,n,1); % 0.5*180.156e-3=0.0901 %[g/mmol]
% cell / mol amino acid (Mee et al. PNAS. 2014)
% gram / cell = 500*10^-15 (cell biology by the numbers)
r.zA1 = 10.^unifrnd(log10(2.7374)-dist,log10(2.7374)+dist,n,1); %(5.4747e15)*(500e-15)*10^-3=2.7374 %[g/mmol] = [cell/mol]*[g/cell]*[mol/mmol] Lysine
r.zA2 = 10.^unifrnd(log10(4.0148)-dist,log10(4.0148)+dist,n,1); %(8.0295e15)*(500e-15)*10^-3=4.0148 %[g/mmol] = [cell/mol]*[g/cell]*[mol/mmol] Isoleucine
% Volume (From experiment)
r.v = 10.^unifrnd(log10(250*10^-6),log10(250*10^-6),n,1); %250*10^-6 [L]

% Diffusion
r1 = log10(3.88*10^-5)-0.5;
r2 = log10(3.88*10^-5)+0.5;
r.d = 10.^unifrnd(r1,r2,n,1);
% Secretion Stoichiometry
y1 = -2;
y2 = 1;
r.y = 10.^unifrnd(y1,y2,n,1);

% Initial conditions
initial_G = 5.55e-3; %[mmol] % Glucose [mmol]: 20% stock (20g/100mL) / MW * Conc. (2mL/100mL) * Vol (250uL)
initial_A1 = 1e-9; %[mmol] % AA1 [mmol]: aa/cell (1.1e8) * cell yield (10^9) / N_A (6.022e23) = 1.83e-4 mmol
initial_A2 = 1e-9; %[mmol] % AA2 [mmol]: aa/cell (7.5e7) * cell yield (10^9) / N_A (6.022e23) = 1.25e-4 mmol
% initial_Bio1 = 1.23e-5; %[g]
% initial_Bio2 = 1.23e-5; %[g] Biomass [g]:

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
%     p.zG = r.zG(I1);
%     % Calculate Theroetical Initial Bio based on Glucos Stoichiometry
%     initial_Bio1 = ((p.zG*initial_G)/(0.4960))*(0.0122); %((final biomass [grams])/(final OD - blank OD))*(initial OD - blank OD) % Numbers from calculate initial biomass script
%     initial_Bio2 = initial_Bio1;
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
%     p.yA1 = r.y(I1);
%     p.yA2 = r.y(I1);
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

% Other ratios
Ratio_K = Sim_PC(:,end)./Sim_K(:,end);
Ratio_I = Sim_PC(:,end)./Sim_I(:,end);
Ratio = mean([Ratio_K,Ratio_I],2);
Ratio_PC = Sim_PC(:,end)./Sim_PC(:,1);

% Check for numerical errors (end)
rows = find(log10(Ratio)<-0.2);
Ratio_K(rows) = [];
Ratio_I(rows) = [];
Ratio(rows) = [];
Ratio_PC(rows) = [];
n = n - length(rows);
r.y(rows) = [];
r.d(rows) = [];

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
scatter(log10(r.d),log10(r.y),10,log10(Ratio),'filled','markerfacealpha',0.6)
hold on

% x = log10(r.d);
% y =log10(r.y);
% z = log10(Ratio);
% [xi,yi] = meshgrid(r1:0.1:r2, y1:0.1:y2);
% zi = griddata(x,y,z,xi,yi,'linear');
% lines = [m_Exp_log_RT_I_K+thresh,m_Exp_log_RT_I_K-thresh];
% [C,h] = contour(xi,yi,zi,lines,'linewidth',lw);
%clabel(C,h)
cbh=colorbar;
%cbh.Ticks = [0:0.5:5];
xlabel('log10(Diffusion Rate [L/hr])')
ylabel('log10(Leakage [mmol/g])')
yticks([-2,-1,0,1])
xticks([r1,log10(3.88e-5),r2])
yticklabels({'-2','-1','0','1'})
xticklabels({num2str(round(r1,3)),num2str(round(log10(3.88e-5),3)),num2str(round(r2,3))})
axis([r1 r2 y1 y2])
set(gca,'fontsize',14)

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
% % leakage rates
% leakage = [0.4,-0.3,-0.7];
% % 0.1 pore size
% diffusion = [log10(2.39*10^-5),log10(2.39*10^-5),log10(2.39*10^-5)];
% 
% blue1 = [0 0 0.7];
% green1 = [0 0.7 0];
% c_pc = [0.7 0.7 0.7];
% fntsz = 12;
% 
% %MaxStep = 0.003;
% 
% for I = 1:length(diffusion)
%     I
%     figure(1)
%     plot(diffusion(I),leakage(I),'r.','markersize',15)
%     
% %     figure(2)
% %     plot(diffusion(I),leakage(I),'r.','markersize',15)
%     
%     p.d = 10^diffusion(I);
%     p.yA1 = 10^leakage(I);
%     p.yA2 = p.yA1;
%     
%     p.vmaxG = r.vmaxG(1);
%     p.vmaxA1 = r.vmaxA1(1);
%     p.vmaxA2 = r.vmaxA2(1);
%     p.kG = r.kG(1);
%     p.kA1 = r.kA1(1);
%     p.kA2 = r.kA2(1);
%     p.zG = r.zG(1);
%     p.zA1 = r.zA1(1);
%     p.zA2 = r.zA2(1);
%     p.v = r.v(1);
% 
%     % opposite wells
%     % Initial Conditions
%     x0 = [initial_G,initial_G,... % Glucose [L,R]
%         initial_A1,initial_A1,... % AA1 [L,R]
%         initial_A2,initial_A2,... % AA2 [L,R]
%         initial_Bio1,0,... % Biomass1 [L,R]
%         0,initial_Bio2]; % Biomass2 [L,R]
%     % ODE
%     options = odeset('NonNegative',[1:10],'MaxStep',MaxStep);
%     [t,x] = ode23tb(@(t,x) f_co_culture_model(t,x,p), tspan, x0, options);
%     % Record Results
%     t3_1 = t;
%     x3_1 = x;
%     
%     % same wells
%     % Initial Conditions
%     x0 = [initial_G,initial_G,... % Glucose [L,R]
%         initial_A1,initial_A1,... % AA1 [L,R]
%         initial_A2,initial_A2,... % AA2 [L,R]
%         initial_Bio1/2,initial_Bio1/2,... % Biomass1 [L,R]
%         initial_Bio2/2,initial_Bio2/2]; % Biomass2 [L,R]
%     % ODE
%     options = odeset('NonNegative',[1:10],'MaxStep',MaxStep);
%     [t,x] = ode23tb(@(t,x) f_co_culture_model(t,x,p), tspan, x0, options);
%     % Record Results
%     t3_2 = t;
%     x3_2 = x;
%     
%     figure(3)
%     ind = (I-1)*2+1;
%     subplot(length(diffusion),2,ind)
%     hold on
%     plot(t,x3_1(:,7),'-','linewidth',1,'color',blue1)
%     plot(t,x3_1(:,9),'-','linewidth',1,'color',green1)
%     plot(t,x3_2(:,7)+x3_2(:,9),'linewidth',1,'color',c_pc,'linestyle',':')
%     xlabel('Time [hrs]')
%     ylabel('Biomass [g]')
%     set(gca,'fontsize',fntsz)
%     set(gca,'yscale','log')
%     axis([0 tspan(end) 0 6*10^-4])
%     set(gca,'xtick',0:24:tspan(end))
%     set(gca,'ytick',[10^-4.5,10^-4,10^-3.5])
%     set(gca,'yticklabels',{'10^{-4.5}','10^{-4}','10^{-3.5}'})
%     set(gca, 'YMinorTick','off')
%     
%     subplot(length(diffusion),2,ind+1)
%     hold on
%     plot(t,x3_1(:,8),'-','linewidth',1,'color',blue1)
%     plot(t,x3_1(:,10),'-','linewidth',1,'color',green1)
%     plot(t,x3_2(:,8)+x3_2(:,10),'linewidth',1,'color',c_pc,'linestyle',':')
%     xlabel('Time [hrs]')
%     set(gca,'fontsize',fntsz)
%     set(gca,'yscale','log')
%     axis([0 tspan(end) 0 6*10^-4])
%     set(gca,'xtick',0:24:tspan(end))
%     set(gca,'ytick',[10^-4.5,10^-4,10^-3.5])
%     set(gca,'yticklabels',{'10^{-4.5}','10^{-4}','10^{-3.5}'})
%     set(gca, 'YMinorTick','off')
%     
%     set(gcf,'Renderer', 'painters', 'Position', [0 0 500 800])
% end

%% ABC Distribution Plot (End Point Ratio Data)
figure(4)
% Plot histogram of all points
h = histogram(log10(r.y),20,'normalization','pdf','facealpha',0.5,'facecolor','k','Orientation','horizontal');
% Plot histogram of all points for which ratio is within threshold
hold on
pass = zeros(length(r.y),1);

thresh = 2*s_Exp_log_RT_I_K;
for I = 1:length(r.y)
    if abs(log10(Ratio(I))-m_Exp_log_RT_I_K)<thresh
    	pass(I) = 1;
    end
end
histogram(log10(r.y(pass==1)),'BinEdges',h.BinEdges,'normalization','pdf','facealpha',0.5,'facecolor','b','Orientation','horizontal')

ylabel('log10(Leakage [mmol/g])')
xlabel('Probability Density Function')


%% format Plots
figure(1)
set(gcf,'renderer','painters','position',[0 0 500 400])

figure(4)
set(gcf,'renderer','painters','position',[0 0 250 400])
