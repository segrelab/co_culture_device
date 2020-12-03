% David Bernstein
% 11/12/19
% Plot Data (Co-culture)

clear all
close all
clc

% Load Data
file = 'Figure 5 and 7 -- Syntrophic CoC.xlsx';
range = 'D28:BK412';
[num,txt,raw] = xlsread(file,range);
t_v = [0:0.25:96];

green = [0 0.7 0];
blue = [0 0 0.7];
c_nc = [1 0 0];
c_pc = [0.7 0.7 0.7];
lw = 1;
lw_c = 1.5;

%% 60 well plot
ymax = 1;
xmax = max(t_v);
figure(1) %478nm (Phenol Red)
for I = 1:6
    for J = 1:10
        ind = (I-1)*10+J;
        subplot(6,10,ind)
        
        % Cross-over
        if ind == 30 || ind == 29 || ind == 20 || ind == 19 || ind == 10 || ind == 9 || ind == 8 || ind == 7 || ind == 24 || ind == 23
            ls = '--';
        else
            ls = '-';
        end
                
        if ind > 30 && ind < 51
            plot(t_v,num(1:end,ind),'color',c_nc,'linewidth',lw,'linestyle',ls);
        elseif ind > 50
            plot(t_v,num(1:end,ind),'color',c_pc,'linewidth',lw,'linestyle',ls);
        else
            if mod(ind,2) == 1
                plot(t_v,num(1:end,ind),'color',blue,'linewidth',lw,'linestyle',ls);
            else
                plot(t_v,num(1:end,ind),'color',green,'linewidth',lw,'linestyle',ls);
            end
        end
        axis([0 xmax 0 ymax])
    end
end
set(gcf,'renderer','painters')
%saveas(gcf,'CC_all.svg')

%% Replicates
figure(2)
ls = '-';
for I = 1:5
    for J = 1:2
        ind1 = (I-1)*2+J;
        subplot(5,2,ind1)
        hold on
        for K = 1:6
            ind2 = ind1 + (K-1)*10;
            if ind2 == 30 || ind2 == 29 || ind2 == 20 || ind2 == 19 || ind2 == 10 || ind2 == 9 || ind2 == 8 || ind2 == 7 || ind2 == 24 || ind2 == 23
                ls = '--';
            else
                ls = '-';
            end
            if K == 6
                plot(t_v,num(1:end,ind2),'color',c_pc,'linewidth',lw_c,'linestyle',':')  
            elseif K == 5
                plot(t_v,num(1:end,ind2),'color',c_nc,'linewidth',lw_c,'linestyle',':')
             elseif K == 4
                plot(t_v,num(1:end,ind2),'color',c_nc,'linewidth',lw_c,'linestyle',':')
            else
                if mod(ind1,2) == 0
                    plot(t_v,num(1:end,ind2),'color',green,'linewidth',lw,'linestyle',ls)
                else
                    plot(t_v,num(1:end,ind2),'color',blue,'linewidth',lw,'linestyle',ls)                    
                end
            end
        end
        axis([0 96 0 ymax])
        set(gca,'fontsize',12)
        xticks(0:24:96)
        yticks([0.25,0.5,1])
        set(gca,'YScale','log')
        yticklabels({'2^{-2}','2^{-1}','2{^0}'});
    end
end

set(gcf,'renderer','painters','Position', [0 0 400 800])
%saveas(gcf,'CoCulture.svg')

%% One Pore Size
lw = 1;
c = [blue;green];
time = t_v;
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
    for ind2 = [BI,BI+10,BI+20]+I-1
        plot(time,num(1:end,ind2),'color',c(I,:),'linewidth',lw)
    end
end
set(gcf,'renderer','painters','Position', [0 0 400 150])
%saveas(gcf,'CC_1.svg')

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
    for ind2 = [1,11,21]+I-1
        plot(time,num(1:end,ind2),'color',c(I,:),'linewidth',lw)
    end
end
set(gcf,'renderer','painters','Position', [0 0 400 150])
%saveas(gcf,'CC_2.svg')

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
        plot(time,num(1:end,ind2),'color','k','linewidth',lw)
    end
end
set(gcf,'renderer','painters','Position', [0 0 400 150])
%saveas(gcf,'CC_3.svg')

% Negative Control 1
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
    for ind2 = [BI+30]+I-1
        plot(time,num(1:end,ind2),'color',blue,'linewidth',lw)
    end
end
set(gcf,'renderer','painters','Position', [0 0 400 150])
%saveas(gcf,'CC_4.svg')

% Negative Control 2
figure(7)
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
        plot(time,num(1:end,ind2),'color',green,'linewidth',lw)
    end
end
set(gcf,'renderer','painters','Position', [0 0 400 150])
%saveas(gcf,'CC_5.svg')

