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
            ms = 15;
        else % right side
            m = '>';
            ms = 7;
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
set(gca,'fontsize',14)
xticks([0.03,0.1,0.2,0.4])
xticklabels({'0.03','0.1','0.2','0.4'})


