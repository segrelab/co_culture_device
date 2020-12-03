% David Bernstein
% y19m11d11
% Read Data (Isoleucine Diffusion)

clear all
close all
clc

% Load Data
file = 'Figure 4 -- Amino Acid Diffusion - Isoleucine.xlsx';
range = 'D29:BK413';
[num1,txt,raw] = xlsread(file,range);
time = [0:0.25:96];

green = [0 0.7 0];
c_nc = [1 0 0];
c_pc = [0.7 0.7 0.7];
lw = 1;
lw_c = 1.5;

%% 60 well plot
% NOTE: Plate arrangement was flipped (membrane pore sizes are flipped left
% to right). Data plotted here does not allign with the protocol

num = num1;
ymax = 1.2;
figure(1)

% Flip Data
num = zeros(size(num1));
i = [5:-1:1];
for I = 1:5
    for J = 1:6
        I1 = (I-1)*2+(J-1)*10+1;
        I2 = (i(I)-1)*2+(J-1)*10+1;
        num(:,I1) = num1(:,I2);
        num(:,I1+1) = num1(:,I2+1);
    end
end

for I = 1:6
    for J = 1:10
        ind1 = (I-1)*10+J;
        subplot(6,10,ind1)
        if ind1 > 40 && ind1 < 51
            plot(time,num(1:end,ind1),'color',c_nc,'linewidth',lw);
        elseif ind1 > 50
            plot(time,num(1:end,ind1),'color',c_pc,'linewidth',lw);
        else
            plot(time,num(1:end,ind1),'color',green,'linewidth',lw);
        end
        axis([0 96 0 ymax]) 
    end
end
set(gcf,'renderer','painters')
%saveas(gcf,'Ile_all.svg')

%% Replicates
% NOTE: Plate arrangement was flipped (membrane pore sizes are flipped left
% to right). Data is plotted here to align with the lysine experiment

% Flip Data
num = zeros(size(num1));
i = [5:-1:1];
for I = 1:5
    for J = 1:6
        I1 = (I-1)*2+(J-1)*10+1;
        I2 = (i(I)-1)*2+(J-1)*10+1;
        num(:,I1) = num1(:,I2);
        num(:,I1+1) = num1(:,I2+1);
    end
end

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
                plot(time,num(1:end,ind2),'color',green,'linewidth',lw)
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
%saveas(gcf,'Ile.svg')
