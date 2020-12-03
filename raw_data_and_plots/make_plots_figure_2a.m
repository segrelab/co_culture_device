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