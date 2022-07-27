% David Bernstein
% y19m11d11
% Read Data (Chemical Diffusion)

% Note: the "fit" function used here is from the curve fitting toolbox,
% this toolbox must be installed or else MATLAB will try to use the stats
% function and give an error.

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

%% Load Calibration Data
file = 'Figure 2 -- Chemical Dye Diffusion - Calibration Curves.xlsx';
range = 'D28:BK52';
[A478_2,~,~] = xlsread(file,range);
range = 'D57:BK81';
[A490_2,~,~] = xlsread(file,range);
t_v2 = [0:0.25:6];

conc = [450:-50:0];
DATA_PR = [];
for I = 1:3
    for J = 1:10
        ind = (I-1)*10+J;
        for K = 6:length(t_v2)
            DATA_PR = [DATA_PR;[conc(J),A478_2(K,ind)]];
        end
    end
end

DATA_BP = [];
for I = 4:6
    for J = 1:10
        ind = (I-1)*10+J;
        for K = 6:length(t_v2)
            DATA_BP = [DATA_BP;[conc(J),A490_2(K,ind)]];
        end
    end
end

% Fit linear regression model
T_PR = array2table(DATA_PR);
T_PR.Properties.VariableNames(1:2) = {'C','A'};
mdl_PR = fitlm(T_PR,'C~A');
T_BP = array2table(DATA_BP);
T_BP.Properties.VariableNames(1:2) = {'C','A'};
mdl_BP = fitlm(T_BP,'C~A');

%% Calculate Diffusion Constants
% Convert data to mmol
% OD to uM
uM_PR = mdl_PR.Coefficients{1,1} + mdl_PR.Coefficients{2,1}.*A478_1;
uM_BP = mdl_BP.Coefficients{1,1} + mdl_BP.Coefficients{2,1}.*A490_1;
% uM to mmol
vol = 250e-6; %[L]
m_PR = (uM_PR/1000)*vol;
m_BP = (uM_PR/1000)*vol;
K = (400/1000)*vol;

% Transform
% mmol-K/2
T_m_PR = m_PR-(K/2);
T_m_BP = m_BP-(K/2);
figure(1)
red = [.9 0 0];
purp = [.7 0 .7];
% Plot transformed data
% (Phenol Red)
ymax = max(max([T_m_PR,T_m_BP]));
ymin = min(min([T_m_PR,T_m_BP]));
xmax = max(t_v1);
hold on
lw = 1;

for I = 1:3
    for J = 1:10
        ind = (I-1)*10+J;
        subplot(6,10,ind)
        plot(t_v1,T_m_PR(1:end,ind),'color',red,'linewidth',lw);
        axis([0 xmax ymin ymax])
    end
end
% (Bromocresol Purple)
for I = 4:6
    for J = 1:10
        ind = (I-1)*10+J;
        subplot(6,10,ind)
        plot(t_v1,T_m_BP(1:end,ind),'color',purp,'linewidth',lw);
        axis([0 xmax ymin ymax])
    end
end

% Fit Exponential
D = zeros(6,10);
figure(1)
% (Phenol Red)
for I = 1:3
    for J = 1:10
        ind = (I-1)*10+J;
        %fit exponential
        x = t_v1';
        y = T_m_PR(:,ind);
        f = fit(x,y,'exp1','lower',[(-1)^(ind+1)*K/2,-Inf],'upper',[(-1)^(ind+1)*K/2,Inf],'TolFun',1e-11);
        subplot(6,10,ind)
        hold on
        X = t_v1';
        Y = f.a*exp(f.b*X);
        plot(X,Y,'k:','linewidth',lw)
        D(I,J) = (f.b*vol)/-2;
    end
end
% (Bromocresol Purple)
for I = 4:6
    for J = 1:10
        ind = (I-1)*10+J;
        %fit exponential
        x = t_v1';
        y = T_m_BP(:,ind);
        f = fit(x,y,'exp1','lower',[(-1)^(ind+1)*K/2,-Inf],'upper',[(-1)^(ind+1)*K/2,Inf],'TolFun',1e-11);
        subplot(6,10,ind)
        hold on
        X = t_v1';
        Y = f.a*exp(f.b*X);
        plot(X,Y,'k:','linewidth',lw)
        D(I,J) = (f.b*vol)/-2;
    end
end
set(gcf,'renderer','painters')
%saveas(gcf,'fit.svg')


%% Plot Diffusion Constants
PS = [0,0.03,0.1,0.2,0.4]; % Pore Size
figure(2)
hold on
for I = 1:size(D,1)
    for J = 1:size(D,2)
        if I <= 3 % phenol red
            c = red;
        else % bromocresol purple
            c = purp;
        end
        if rem(J,2) == 1 % left side
            m = '.';
            ms = 20;
        else % right side
            m = '>';
            ms = 10;
        end
        X = PS(floor((J-1)/2)+1);
        plot(X,D(I,J),'color',c,'marker',m,'markersize',ms)
    end
end
ylabel('Diffusion Constant [L/hr]')
xlabel('Pore Size [{\mu}m]')

%% Average Diffusion Constants
m_D(1) = mean([D(:,1);D(:,2)]);
m_D(2) = mean([D(:,3);D(:,4)]);
m_D(3) = mean([D(:,5);D(:,6)]);
m_D(4) = mean([D(:,7);D(:,8)]);
m_D(5) = mean([D(:,9);D(:,10)]);

s_D(1) = std([D(:,1);D(:,2)]);
s_D(2) = std([D(:,3);D(:,4)]);
s_D(3) = std([D(:,5);D(:,6)]);
s_D(4) = std([D(:,7);D(:,8)]);
s_D(5) = std([D(:,9);D(:,10)]);

MAT = [PS',m_D',s_D'];

% Log Distribution
l_m_D(1) = mean(log10([D(:,1);D(:,2)]));
l_m_D(2) = mean(log10([D(:,3);D(:,4)]));
l_m_D(3) = mean(log10([D(:,5);D(:,6)]));
l_m_D(4) = mean(log10([D(:,7);D(:,8)]));
l_m_D(5) = mean(log10([D(:,9);D(:,10)]));

l_s_D(1) = std(log10([D(:,1);D(:,2)]));
l_s_D(2) = std(log10([D(:,3);D(:,4)]));
l_s_D(3) = std(log10([D(:,5);D(:,6)]));
l_s_D(4) = std(log10([D(:,7);D(:,8)]));
l_s_D(5) = std(log10([D(:,9);D(:,10)]));

LOG_MAT = [PS',l_m_D',l_s_D'];

%% Save Figure
figure(2)
set(gcf,'renderer','painters')
set(gca,'fontsize',16)
xticks([0,0.03,0.1,0.2,0.4])
xticklabels({'0','0.03','0.1','0.2','0.4'})
axis([0,0.4,0,10e-5])

%% Lysine and Isoleucine
% Load data and calculate concentrations

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
% Lysine
L_st = [100,50,25,12.5,6.25,0]; %[micro molar]
% Average standard
%a_L_s_c = mean(L_s_c,'all'); %blanks
a_L_s_c = mean([L_s(1,6),L_s(2,6)]); %0 conc samples
% Keep track of data for linear regression
DATA_L = [];
for I = 1:6 %concentrations
    for J = 1:2 %replicates
        DATA_L = [DATA_L;[L_st(I),L_s(J,I)-a_L_s_c]];
    end
end
% Linear regression model
T_L = array2table(DATA_L);
T_L.Properties.VariableNames(1:2) = {'C','A'};
mdl_L = fitlm(T_L,'C~A');
X = [min(DATA_L(:,2)),max(DATA_L(:,2))];
Y = mdl_L.Coefficients{1,1} + mdl_L.Coefficients{2,1}.*X;

% Isoleucine
I_st = [1000,500,250,125,62.5,31.25,15.625,0]; %[micro molar]
% Average standard
%a_I_s_c = mean(I_s_c,'all'); %blanks
a_I_s_c = mean([I_s(1,8),I_s(2,8),I_s(3,8)]); %0 conc samples
% Keep track of data for linear regression
DATA_I = [];
for I = 1:8 %concentrations
    for J = 1:2 %replicates
        DATA_I = [DATA_I;[I_st(I),I_s(J,I)-a_I_s_c]];
    end
end
% Linear regression model
T_I = array2table(DATA_I);
T_I.Properties.VariableNames(1:2) = {'C','A'};
mdl_I = fitlm(T_I,'C~A');
X = [min(DATA_I(:,2)),max(DATA_I(:,2))];
Y = mdl_I.Coefficients{1,1} + mdl_I.Coefficients{2,1}.*X;

% Plot Amino Acid Concentration Curves
time = [0,1.5,3,6,12,24]; %hours
% Lysine
% Convert data to micro molar
% Normalize with blank
a_L_c = a_L_s_c; %Use blank from standard curve
%a_L_c = mean(L_d_c,'all'); %Use blank from data
L_d_n = L_d - a_L_c;
% Convert to micro molar with linear regression
L_d_um = (L_d_n.*mdl_L.Coefficients{2,1}+mdl_L.Coefficients{1,1})*4;

% Isoleucine
% Convert data to micro molar
% Normalize with blank
a_I_c = a_I_s_c; %Use blank from standard curve
%a_I_c = mean(I_d_c,'all'); %Use blank from data
I_d_n = I_d - a_I_c;
% Convert to micro molar with linear regression
I_d_um = I_d_n.*mdl_I.Coefficients{2,1}+mdl_I.Coefficients{1,1};


% Fit to exponential
% uM to mmol
vol = 250e-6; %[L]
m_L = (L_d_um/1000)*vol;
m_I = (I_d_um/1000)*vol;
K_L = (mean([L_d_um(1,1),L_d_um(2,1),L_d_um(3,1)])/1000)*vol;
K_I = (mean([I_d_um(1,1),I_d_um(2,1),I_d_um(3,1)])/1000)*vol;

% Transform
% mmol-K/2
T_m_L = m_L-(K_L/2);
T_m_I = m_I-(K_I/2);

% Plot transformed data
figure(3)
% Lysine
subplot(1,2,1)
hold on
for I = 1:6 %time points
    l_ind = I*2-1;
    r_ind = I*2;
    for J = 1:3 %replicates
        % left well
        plot(time(I),T_m_L(J,l_ind),'b.','markersize',10)
        % right well
        plot(time(I),T_m_L(J,r_ind),'r.','markersize',10)
        % total
        plot(time(I),T_m_L(J,l_ind)+T_m_L(J,r_ind),'k.','markersize',10)
    end
end
subplot(1,2,2)
hold on
for I = 1:6 %time points
    l_ind = I*2-1;
    r_ind = I*2;
    for J = 1:3 %replicates
        % left well
        plot(time(I),T_m_I(J,l_ind),'b.','markersize',10)
        % right well
        plot(time(I),T_m_I(J,r_ind),'r.','markersize',10)
        % total
        plot(time(I),T_m_I(J,l_ind)+T_m_I(J,r_ind),'k.','markersize',10)
    end
end

% Fit exponential
time = [0,1.5,3,6,12,24];
x = [time,time,time];
X = linspace(min(time),max(time),100);

% Lysine
% Left
y = [T_m_L(1,[1,3,5,7,9,11]),T_m_L(2,[1,3,5,7,9,11]),T_m_L(3,[1,3,5,7,9,11])];
f = fit(x',y','exp1','lower',[K_L/2,-Inf],'upper',[K_L/2,Inf],'TolFun',1e-11);
D_L_l = (f.b*vol)/-2;
% plot exponential fit
Y = f.a*exp(f.b*X);
subplot(1,2,1)
plot(X,Y,'b:','linewidth',2)
% Right
y = [T_m_L(1,[2,4,6,8,10,12]),T_m_L(2,[2,4,6,8,10,12]),T_m_L(3,[2,4,6,8,10,12])];
f = fit(x',y','exp1','lower',[-1*K_L/2,-Inf],'upper',[-1*K_L/2,Inf],'TolFun',1e-11);
D_L_r = (f.b*vol)/-2;
% plot exponential fit
Y = f.a*exp(f.b*X);
subplot(1,2,1)
plot(X,Y,'r:','linewidth',2)

% Isoleucine
% Left
y = [T_m_I(1,[1,3,5,7,9,11]),T_m_I(2,[1,3,5,7,9,11]),T_m_I(3,[1,3,5,7,9,11])];
f = fit(x',y','exp1','lower',[K_I/2,-Inf],'upper',[K_I/2,Inf],'TolFun',1e-11);
D_I_l = (f.b*vol)/-2;
% plot exponential fit
Y = f.a*exp(f.b*X);
subplot(1,2,2)
plot(X,Y,'b:','linewidth',2)
% Right
y = [T_m_I(1,[2,4,6,8,10,12]),T_m_I(2,[2,4,6,8,10,12]),T_m_I(3,[2,4,6,8,10,12])];
f = fit(x',y','exp1','lower',[-1*K_I/2,-Inf],'upper',[-1*K_I/2,Inf],'TolFun',1e-11);
D_I_r = (f.b*vol)/-2;
% plot exponential fit
Y = f.a*exp(f.b*X);
subplot(1,2,2)
plot(X,Y,'r:','linewidth',2)

% Plote diffusion coefficients on figure
figure(2)
plot(0.1,D_L_l,'b.','markersize',20)
plot(0.1,D_L_r,'b>','markersize',10)
plot(0.1,D_I_l,'g.','markersize',20)
plot(0.1,D_I_r,'g>','markersize',10)
