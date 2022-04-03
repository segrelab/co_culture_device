% David Bernstein
% y19m11d11
% Read Data (Chemical Diffusion)

clear all
close all
clc

% Load Data
file = 'Figure 2 -- Chemical Dye Diffusion - Phenol Red & Bromocresol Purple.xlsx';
range = 'D28:BK316';
[A478_1,txt1,raw1] = xlsread(file,range);
range = 'D321:BK609';
[A490_1,txt2,raw2] = xlsread(file,range);
t_v1 = [0:0.25:72];

%% Plot OD
red = [.9 0 0];
purp = [.7 0 .7];
lw = 1;
% 60 well plot 
ymax = max(max([A478_1,A490_1]));
xmax = max(t_v1);
figure(1)
%478nm (Phenol Red)
hold on
for I = 1:3
    for J = 1:10
        ind = (I-1)*10+J;
        subplot(6,10,ind)
        plot(t_v1,A478_1(1:end,ind),'color',red,'linewidth',lw);
        axis([0 xmax 0 ymax])
    end
end
%490nm (Bromocresol Purple)
for I = 4:6
    for J = 1:10
        ind = (I-1)*10+J;
        subplot(6,10,ind)
        plot(t_v1,A490_1(1:end,ind),'color',purp,'linewidth',lw);
        axis([0 xmax 0 ymax])
    end
end

%% Load Calibration Data
file = 'Figure 2 -- Chemical Dye Diffusion - Calibration Curves.xlsx';
range = 'D28:BK52';
[A478_2,~,~] = xlsread(file,range);
range = 'D57:BK81';
[A490_2,~,~] = xlsread(file,range);
t_v2 = [0:0.25:6];
ymax = 3.3;
xmax = max(t_v2);
figure(2)
hold on
% Phenol Red Abs478
for I = 1:3
    for J = 1:10
        ind = (I-1)*10+J;
        subplot(6,10,ind)
        plot(t_v2,A478_2(1:end,ind),'color',red,'linewidth',lw);
        axis([0 xmax 0 ymax])
    end
end
figure(2)
% Bromocresol Purple Abs490
for I = 4:6
    for J = 1:10
        ind = (I-1)*10+J;
        subplot(6,10,ind)
        plot(t_v2,A490_2(1:end,ind),'color',purp,'linewidth',lw);
        axis([0 xmax 0 ymax])
    end
end

% Plot Abs vs Conc.
conc = [450:-50:0];
figure(3)
subplot(1,2,1)
hold on
title('Phenol Red')
DATA_PR = [];
for I = 1:3
    for J = 1:10
        ind = (I-1)*10+J;
        for K = 6:length(t_v2)
            plot(A478_2(K,ind),conc(J),'k.')
            DATA_PR = [DATA_PR;[conc(J),A478_2(K,ind)]];
        end
    end
end
xlabel('OD 478')
ylabel('Conc. uM')
set(gca,'fontsize',16)

figure(3)
subplot(1,2,2)
hold on
title('Bromocresol Purple')
DATA_BP = [];
for I = 4:6
    for J = 1:10
        ind = (I-1)*10+J;
        for K = 6:length(t_v2)
            plot(A490_2(K,ind),conc(J),'k.')
            DATA_BP = [DATA_BP;[conc(J),A490_2(K,ind)]];
        end
    end
end
xlabel('OD 490')
ylabel('Conc. uM')
set(gca,'fontsize',16)

% Fit linear regression model
T_PR = array2table(DATA_PR);
T_PR.Properties.VariableNames(1:2) = {'C','A'};
mdl_PR = fitlm(T_PR,'C~A');
T_BP = array2table(DATA_BP);
T_BP.Properties.VariableNames(1:2) = {'C','A'};
mdl_BP = fitlm(T_BP,'C~A');
subplot(1,2,1)
X = [min(DATA_PR(:,2)),max(DATA_PR(:,2))];
Y = mdl_PR.Coefficients{1,1} + mdl_PR.Coefficients{2,1}.*X;
plot(X,Y,'color',red,'linewidth',lw)
text(1,50,strcat('R-squared:',num2str(round(mdl_PR.Rsquared.Ordinary,3))),'fontsize',14)
subplot(1,2,2)
X = [min(DATA_BP(:,2)),max(DATA_BP(:,2))];
Y = mdl_BP.Coefficients{1,1} + mdl_BP.Coefficients{2,1}.*X;
plot(X,Y,'color',purp,'linewidth',lw)
text(0.8,50,strcat('R-squared:',num2str(round(mdl_BP.Rsquared.Ordinary,3))),'fontsize',14)

set(gcf,'renderer','painters')

% Convert data to mmol
% OD to uM
uM_PR = mdl_PR.Coefficients{1,1} + mdl_PR.Coefficients{2,1}.*A478_1;
uM_BP = mdl_BP.Coefficients{1,1} + mdl_BP.Coefficients{2,1}.*A490_1;
% uM to mmol
vol = 250e-6; %[L]
m_PR = (uM_PR/1000)*vol;
m_BP = (uM_PR/1000)*vol;
K = (400/1000)*vol;

%% Plot mmol
figure(4)
% (Phenol Red)
ymax = max(max([m_PR,m_BP]));
ymin = min(min([m_PR,m_BP]));
xmax = max(t_v1);
hold on
for I = 1:3
    for J = 1:10
        ind = (I-1)*10+J;
        subplot(6,10,ind)
        plot(t_v1,m_PR(1:end,ind),'r');
        axis([0 xmax ymin ymax])
    end
end
% (Bromocresol Purple)
for I = 4:6
    for J = 1:10
        ind = (I-1)*10+J;
        subplot(6,10,ind)
        plot(t_v1,m_BP(1:end,ind),'m');
        axis([0 xmax ymin ymax])
    end
end

% Phenol Red
figure(5)
c = [linspace(0.9,0.9,5)',linspace(0,0,5)',linspace(0,0,5)'];
lw = 1.5;
lw2 = 0.5;
for J = 1:5
    for K = 1:2
        ind = K+(J-1)*2;
        subplot(5,2,K+(J-1)*2);
        hold on
        inds = [ind+(1-1)*10,ind+(2-1)*10,ind+(3-1)*10];
        % Calculate mean
        m = mean(m_PR(:,inds),2);
        % Calculate standard error
        s = (std(m_PR(:,inds),'',2)./sqrt(3));
        m_ns = m-s;
        m_ps = m+s;
%         patch([t_v1,fliplr(t_v1)],[m_ns',fliplr(m_ps')],c(J,:),'facealpha',0.5,'edgecolor','none');
%         plot(t_v1,m,'color',c(J,:),'linewidth',lw);
%         plot(t_v1,m_ns,'color',c(J,:),'linewidth',lw2);
%         plot(t_v1,m_ps,'color',c(J,:),'linewidth',lw2);
        for I = 1:length(inds)
            plot(t_v1,m_PR(:,inds(I)),'color',c(J,:),'linewidth',lw2)
        end
        axis([0 xmax ymin ymax])
    end
end


% Bromocresol Purple
figure(6)
c = [linspace(0.7,0.7,5)',linspace(0,0,5)',linspace(0.7,0.7,5)'];
 lw = 1.5;
 lw2 = 0.5;
for J = 1:5
    for K = 1:2
        ind = K+(J-1)*2+30;
        subplot(5,2,K+(J-1)*2);
        hold on
        inds = [ind+(1-1)*10,ind+(2-1)*10,ind+(3-1)*10];
        % Calculate mean
        m = mean(m_BP(:,inds),2);
        % Calculate standard error
        s = (std(m_BP(:,inds),'',2)./sqrt(3));
        m_ns = m-s;
        m_ps = m+s;
%         patch([t_v1,fliplr(t_v1)],[m_ns',fliplr(m_ps')],c(J,:),'facealpha',0.5,'edgecolor','none');
%         plot(t_v1,m,'color',c(J,:),'linewidth',lw);
%         plot(t_v1,m_ns,'color',c(J,:),'linewidth',lw2);
%         plot(t_v1,m_ps,'color',c(J,:),'linewidth',lw2);
        for I = 1:length(inds)
            plot(t_v1,m_BP(:,inds(I)),'color',c(J,:),'linewidth',lw2)
        end
        axis([0 xmax ymin ymax])
    end
end

%% Lysine and Isoleucine Data

% Load Data
% Lysine
sheet = 1;
file = '/Users/david/Documents/GitHub/co_culture_device/raw_data_and_plots/Figure 2 -- AssayAnalysis-BioMeDiffusion.xlsx';
% Standard curve
range = 'C29:H30';
[L_s,~,~] = xlsread(file,sheet,range);
% Standard curve control
range = 'C35:H36';
[L_s_c,~,~] = xlsread(file,sheet,range);
% Data
range = 'C46:N48';
[L_d,~,~] = xlsread(file,sheet,range);
% Data control
range = 'C53:N55';
[L_d_c,~,~] = xlsread(file,sheet,range);

% Isoleucine
sheet = 2;
% Standard curve
range = 'C31:J33';
[I_s,~,~] = xlsread(file,sheet,range);
% Standard curve control
range = 'C38:J39';
[I_s_c,~,~] = xlsread(file,sheet,range);
% Data
range = 'C49:N51';
[I_d,~,~] = xlsread(file,sheet,range);
% Data control
range = 'C56:N58';
[I_d_c,~,~] = xlsread(file,sheet,range);


% Plot Standard Curve
figure(7)
% Lysine
L_st = [100,50,25,12.5,6.25,0]; %[micro molar]
subplot(1,2,1)
hold on
% Average standard
%a_L_s_c = mean(L_s_c,'all'); %blanks
a_L_s_c = mean([L_s(1,6),L_s(2,6)]); %0 conc samples
% Keep track of data for linear regression
DATA_L = [];
for I = 1:6 %concentrations
    for J = 1:2 %replicates
        plot(L_s(J,I)-a_L_s_c,L_st(I),'k.','markersize',10)
        DATA_L = [DATA_L;[L_st(I),L_s(J,I)-a_L_s_c]];
    end
end
% Linear regression model
T_L = array2table(DATA_L);
T_L.Properties.VariableNames(1:2) = {'C','A'};
mdl_L = fitlm(T_L,'C~A');
subplot(1,2,1)
X = [min(DATA_L(:,2)),max(DATA_L(:,2))];
Y = mdl_L.Coefficients{1,1} + mdl_L.Coefficients{2,1}.*X;
plot(X,Y,'color',red,'linewidth',lw)
text(0.1,20,strcat('R-squared:',num2str(round(mdl_L.Rsquared.Ordinary,3))),'fontsize',14)
% Format Plot
ylim([-1 120])
title('Lysine')
xlabel('OD 550')
ylabel('Conc. uM')
set(gca,'fontsize',16)
set(gcf,'renderer','painters')


% Isoleucine
I_st = [1000,500,250,125,62.5,31.25,15.625,0]; %[micro molar]
subplot(1,2,2)
hold on
% Average standard
%a_I_s_c = mean(I_s_c,'all'); %blanks
a_I_s_c = mean([I_s(1,8),I_s(2,8),I_s(3,8)]); %0 conc samples
% Keep track of data for linear regression
DATA_I = [];
for I = 1:8 %concentrations
    for J = 1:2 %replicates
        plot(I_s(J,I)-a_I_s_c,I_st(I),'k.','markersize',10)
        DATA_I = [DATA_I;[I_st(I),I_s(J,I)-a_I_s_c]];
    end
end
% Linear regression model
T_I = array2table(DATA_I);
T_I.Properties.VariableNames(1:2) = {'C','A'};
mdl_I = fitlm(T_I,'C~A');
subplot(1,2,2)
X = [min(DATA_I(:,2)),max(DATA_I(:,2))];
Y = mdl_I.Coefficients{1,1} + mdl_I.Coefficients{2,1}.*X;
plot(X,Y,'color',red,'linewidth',lw)
text(0.5,100,strcat('R-squared:',num2str(round(mdl_I.Rsquared.Ordinary,3))),'fontsize',14)
% Format Plot
ylim([-20 1000])
title('Isoleucine')
xlabel('OD 450')
ylabel('Conc. uM')
set(gca,'fontsize',16)
set(gcf,'renderer','painters')

% Plot Amino Acid Concentration Curves
figure(8)
time = [0,1.5,3,6,12,24]; %hours
% Lysine
% Convert data to micro molar
% Normalize with blank
a_L_c = a_L_s_c; %Use blank from standard curve
%a_L_C = mean(L_d_c,'all'); %Use blank from data
L_d_n = L_d - a_L_c;
% Convert to micro molar with linear regression
L_d_um = (L_d_n.*mdl_L.Coefficients{2,1}+mdl_L.Coefficients{1,1})*4;
% Plot data
subplot(1,2,1)
hold on
for I = 1:6 %time points
    l_ind = I*2-1;
    r_ind = I*2;
    for J = 1:3 %replicates
        % left well
        plot(time(I),L_d_um(J,l_ind),'b.','markersize',10)
        % right well
        plot(time(I),L_d_um(J,r_ind),'r.','markersize',10)
        % total
        plot(time(I),L_d_um(J,l_ind)+L_d_um(J,r_ind),'k.','markersize',10)
    end
end
% Plot means
m_L_d_um = mean(L_d_um,1);
m_L_d_um_left = m_L_d_um([1,3,5,7,9,11]);
m_L_d_um_right = m_L_d_um([2,4,6,8,10,12]);
m_L_d_um_total = m_L_d_um_left+m_L_d_um_right;
plot(time,m_L_d_um_left,'b:','linewidth',2)
plot(time,m_L_d_um_right,'r:','linewidth',2)
plot(time,m_L_d_um_total,'k:','linewidth',2)
% Format Plot
title('Lysine')
xlabel('Time [hours]')
ylabel('Conc. [uM]')
set(gca,'fontsize',16)
set(gcf,'renderer','painters')

% Isoleucine
% Convert data to micro molar
% Normalize with blank
a_I_c = a_I_s_c; %Use blank from standard curve
%a_L_C = mean(L_d_c,'all'); %Use blank from data
I_d_n = I_d - a_I_c;
% Convert to micro molar with linear regression
I_d_um = I_d_n.*mdl_I.Coefficients{2,1}+mdl_I.Coefficients{1,1};
% Plot data
subplot(1,2,2)
hold on
for I = 1:6 %time points
    l_ind = I*2-1;
    r_ind = I*2;
    for J = 1:3 %replicates
        % left well
        plot(time(I),I_d_um(J,l_ind),'b.','markersize',10)
        % right well
        plot(time(I),I_d_um(J,r_ind),'r.','markersize',10)
        % total
        plot(time(I),I_d_um(J,l_ind)+I_d_um(J,r_ind),'k.','markersize',10)
    end
end
% Plot means
m_I_d_um = mean(I_d_um,1);
m_I_d_um_left = m_I_d_um([1,3,5,7,9,11]);
m_I_d_um_right = m_I_d_um([2,4,6,8,10,12]);
m_I_d_um_total = m_I_d_um_left+m_I_d_um_right;
plot(time,m_I_d_um_left,'b:','linewidth',2)
plot(time,m_I_d_um_right,'r:','linewidth',2)
plot(time,m_I_d_um_total,'k:','linewidth',2)
% Format Plot
title('Isoleucine')
xlabel('Time [hours]')
ylabel('Conc. [uM]')
set(gca,'fontsize',16)
set(gcf,'renderer','painters')


