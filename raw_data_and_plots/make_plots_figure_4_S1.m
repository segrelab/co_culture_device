% David Bernstein
% 2/26/20
% growth curve data

clear all
close all
clc

Range = 'D25:O313';
[num,txt,raw] = xlsread('Figure 4 -- EC NR1 Auxotroph Growth Curves.xlsx',Range);

figure(1)
hold on
blue = [0 0 0.7];
green = [0 0.7 0];
lw = 1;
t_v = 0:0.25:72;
for I = 1:4
    plot(t_v,num(:,I),'-','color',blue,'linewidth',lw);
end
for I = 5:8
    plot(t_v,num(:,I),'-','color',green,'linewidth',lw);
end
% Plot dotted line at 6 hours
plot([6,6],[0.09,1.4],'k:','linewidth',1)
xticks(0:24:72)
ylabel('OD 600')
xlabel('Time [hours]')
set(gca,'YScale','log')
set(gca,'fontsize',14)

