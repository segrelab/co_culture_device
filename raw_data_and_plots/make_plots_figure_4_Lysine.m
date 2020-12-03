% David Bernstein
% y19m11d11
% Read Data (Lysine Diffusion)

clear all
close all
clc

% Load Data
file = 'Figure 4 -- Amino Acid Diffusion - Lysine.xlsx';
range = 'D29:BK413';
[num1,txt,raw] = xlsread(file,range);
time = [0:0.25:96];

blue = [0 0 0.7];
c_nc = [1 0 0];
c_pc = [0.7 0.7 0.7];
lw = 1;
lw_c = 1.5;

%% 60 well plot
num = num1;
ymax = 1.2;
figure(1)
for I = 1:6
    for J = 1:10
        ind1 = (I-1)*10+J;
        subplot(6,10,ind1)
        if ind1 > 40 && ind1 < 51
            plot(time,num(1:end,ind1),'color',c_nc,'linewidth',lw);
        elseif ind1 > 50
            plot(time,num(1:end,ind1),'color',c_pc,'linewidth',lw);
        else
            plot(time,num(1:end,ind1),'color',blue,'linewidth',lw);
        end
        axis([0 96 0 ymax]) 
    end
end
set(gcf,'renderer','painters')
%saveas(gcf,'Lys_all.svg')

%% Replicates
figure(2)
for I = 1:5
    for J = 1:2
        ind1 = (I-1)*2+J;
        subplot(5,2,ind1)
        hold on
        for K = 1:6
            ind2 = ind1 + (K-1)*10;
            if K == 6
                plot(time,num(1:end,ind2),':','color',c_pc,'linewidth',lw_c)  
            elseif K == 5
                plot(time,num(1:end,ind2),':','color',c_nc,'linewidth',lw_c)
            else
                plot(time,num(1:end,ind2),'color',blue,'linewidth',lw)
            end
        end
        axis([0 96 0 ymax])
        set(gca,'YScale','log')
        set(gca,'fontsize',12)
        xticks(0:24:96)
        yticks([0.25,0.5,1])
        yticklabels({'2^{-2}','2^{-1}','2{^0}'});
    end
end

set(gcf,'renderer','painters','Position', [0 0 400 800])
%saveas(gcf,'Lys.svg')

%% One Pore Size
lw = 1;
% 0.1 micron pores size
BI = 5;
% Experimental
figure(3)
for I = 1:2
    subplot(1,2,I)
    hold on
    axis([0 96 0 ymax])
    set(gca,'YScale','log')
    set(gca,'fontsize',12)
    xticks(0:24:96)
    yticks([0.25,0.5,1])
    yticklabels({'2^{-2}','2^{-1}','2{^0}'});
    for ind2 = [BI,BI+10,BI+20,BI+30]+I-1
        plot(time,num(1:end,ind2),'color',blue,'linewidth',lw)
    end
end
set(gcf,'renderer','painters','Position', [0 0 400 150])
%saveas(gcf,'Lys_1.svg')

% No Pores
figure(4)
for I = 1:2
    subplot(1,2,I)
    hold on
    axis([0 96 0 ymax])
    set(gca,'YScale','log')
    set(gca,'fontsize',12)
    xticks(0:24:96)
    yticks([0.25,0.5,1])
    yticklabels({'2^{-2}','2^{-1}','2{^0}'});
    for ind2 = [1,11,21,31]+I-1
        plot(time,num(1:end,ind2),'color',blue,'linewidth',lw)
    end
end
set(gcf,'renderer','painters','Position', [0 0 400 150])
%saveas(gcf,'Lys_2.svg')

% Positive Control
figure(5)
for I = 1:2
    subplot(1,2,I)
    hold on
    axis([0 96 0 ymax])
    set(gca,'YScale','log')
    set(gca,'fontsize',12)
    xticks(0:24:96)
    yticks([0.25,0.5,1])
    yticklabels({'2^{-2}','2^{-1}','2{^0}'});
    for ind2 = [BI+50]+I-1
        plot(time,num(1:end,ind2),'color',blue,'linewidth',lw)
    end
end
set(gcf,'renderer','painters','Position', [0 0 400 150])
%saveas(gcf,'Lys_3.svg')

% Negative Control
figure(6)
for I = 1:2
    subplot(1,2,I)
    hold on
    axis([0 96 0 ymax])
    set(gca,'YScale','log')
    set(gca,'fontsize',12)
    xticks(0:24:96)
    yticks([0.25,0.5,1])
    yticklabels({'2^{-2}','2^{-1}','2{^0}'});
    for ind2 = [BI+40]+I-1
        plot(time,num(1:end,ind2),'color',blue,'linewidth',lw)
    end
end
set(gcf,'renderer','painters','Position', [0 0 400 150])
%saveas(gcf,'Lys_4.svg')
