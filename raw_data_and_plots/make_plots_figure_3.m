% David Bernstein
% 4/6/20
% Plot Data (Drosophila gut)

clear all
close all
clc

% Load Data
file = 'Figure 3 -- DM Gut Microbiome CoC.xlsx';
range = 'D28:BK412';
[num,txt,raw] = xlsread(file,range);
t_v = [0:0.25:96];

yellow = [0.8 0.8 0]; % LP
orange = [0.8 0.4 0]; % LB
purple = [0.7 0 0.7]; % AO

%% 60 well plot
lw = 1;
ls = '-';

ymax = 2;
xmax = max(t_v);
figure(1) %478nm (Phenol Red)
for I = 1:6
    for J = 1:10
        ind = (I-1)*10+J;
        subplot(6,10,ind) 
        % AO
        if ismember(ind,[1,11,21,31,41,51,13,33,53,7,17,37,47,57,38,48,58])
            plot(t_v,num(1:end,ind),'color',purple,'linewidth',lw,'linestyle',ls);            
        % LP   
        elseif ismember(ind,[12,32,52,5,15,25,35,45,55,27,8,9,19,29,10,20,30])
            plot(t_v,num(1:end,ind),'color',yellow,'linewidth',lw,'linestyle',ls);                        
        % LB  
        elseif ismember(ind,[4,14,24,34,44,54,16,36,56,18,28,39,49,59,40,50,60])
            plot(t_v,num(1:end,ind),'color',orange,'linewidth',lw,'linestyle',ls);                        
        % empty
        else
            plot(t_v,num(1:end,ind),'color','k','linewidth',lw,'linestyle',ls);                        
        end
        axis([0 xmax 0 ymax])
    end
end
set(gcf,'renderer','painters')
%saveas(gcf,'DM_all.svg')

%% Individual Plots
figure(2)

lw = 0.5;

% AO, LP
Inds = [11,31,51,7];
c = [purple;yellow];
for I = 1:2
    subplot(3,6,0+I)
    hold on
    for ind2 = Inds+I-1
        plot(t_v,num(1:end,ind2),'color',c(I,:),'linewidth',lw)
    end
    axis([0 96 0 ymax])
    set(gca,'YScale','log')
    set(gca,'fontsize',12)
    xticks(0:24:96)
    yticks([0.25,0.5,1,2])
    yticklabels({'2^{-2}','2^{-1}','2{^0}','2'});
    xlabel('Time [hrs]')
    if I == 1
        ylabel('OD 600')
    end
end

% AO, LB
Inds = [13,33,53,17];
c = [purple;orange];
for I = 1:2
    subplot(3,6,2+I)
    hold on
    axis([0 96 0 ymax])
    set(gca,'YScale','log')
    set(gca,'fontsize',12)
    xticks(0:24:96)
    yticks([0.25,0.5,1,2])
    yticklabels({'2^{-2}','2^{-1}','2{^0}','2'});
    for ind2 = Inds+I-1
        plot(t_v,num(1:end,ind2),'color',c(I,:),'linewidth',lw)
    end
    xlabel('Time [hrs]')
    if I == 1
        ylabel('OD 600')
    end
end

% LP, LB
Inds = [15,35,55,27];
c = [yellow;orange];
for I = 1:2
    subplot(3,6,4+I)
    hold on
    axis([0 96 0 ymax])
    set(gca,'YScale','log')
    set(gca,'fontsize',12)
    xticks(0:24:96)
    yticks([0.25,0.5,1,2])
    yticklabels({'2^{-2}','2^{-1}','2{^0}','2'});
    for ind2 = Inds+I-1
        plot(t_v,num(1:end,ind2),'color',c(I,:),'linewidth',lw)
    end
    xlabel('Time [hrs]')
    if I == 1
        ylabel('OD 600')
    end
end

% AO, AO
Inds = [37,47,57];
c = [purple;purple];
for I = 1:2
    subplot(3,6,6+I)
    hold on
    axis([0 96 0 ymax])
    set(gca,'YScale','log')
    set(gca,'fontsize',12)
    xticks(0:24:96)
    yticks([0.25,0.5,1,2])
    yticklabels({'2^{-2}','2^{-1}','2{^0}','2'});
    for ind2 = Inds+I-1
        plot(t_v,num(1:end,ind2),'color',c(I,:),'linewidth',lw)
    end
    xlabel('Time [hrs]')
    if I == 1
        ylabel('OD 600')
    end
end

% LP, LP
Inds = [9,19,29];
c = [yellow;yellow];
for I = 1:2
    subplot(3,6,8+I)
    hold on
    axis([0 96 0 ymax])
    set(gca,'YScale','log')
    set(gca,'fontsize',12)
    xticks(0:24:96)
    yticks([0.25,0.5,1,2])
    yticklabels({'2^{-2}','2^{-1}','2{^0}','2'});
    for ind2 = Inds+I-1
        plot(t_v,num(1:end,ind2),'color',c(I,:),'linewidth',lw)
    end
    xlabel('Time [hrs]')
    if I == 1
        ylabel('OD 600')
    end
end

% LB, LB
Inds = [39,49,59];
c = [orange;orange];
for I = 1:2
    subplot(3,6,10+I)
    hold on
    axis([0 96 0 ymax])
    set(gca,'YScale','log')
    set(gca,'fontsize',12)
    xticks(0:24:96)
    yticks([0.25,0.5,1,2])
    yticklabels({'2^{-2}','2^{-1}','2{^0}','2'});
    for ind2 = Inds+I-1
        plot(t_v,num(1:end,ind2),'color',c(I,:),'linewidth',lw)
    end
    xlabel('Time [hrs]')
    if I == 1
        ylabel('OD 600')
    end
end

% AO, -
Inds = [1,21,41];
c = [purple;[0 0 0]];
for I = 1:2
    subplot(3,6,12+I)
    hold on
    axis([0 96 0 ymax])
    set(gca,'YScale','log')
    set(gca,'fontsize',12)
    xticks(0:24:96)
    yticks([0.25,0.5,1,2])
    yticklabels({'2^{-2}','2^{-1}','2{^0}','2'});
    for ind2 = Inds+I-1
        plot(t_v,num(1:end,ind2),'color',c(I,:),'linewidth',lw)
    end
    xlabel('Time [hrs]')
    if I == 1
        ylabel('OD 600')
    end
end

% LP, -
Inds = [5,25,45];
c = [yellow;[0 0 0]];
for I = 1:2
    subplot(3,6,14+I)
    hold on
    axis([0 96 0 ymax])
    set(gca,'YScale','log')
    set(gca,'fontsize',12)
    xticks(0:24:96)
    yticks([0.25,0.5,1,2])
    yticklabels({'2^{-2}','2^{-1}','2{^0}','2'});
    for ind2 = Inds+I-1
        plot(t_v,num(1:end,ind2),'color',c(I,:),'linewidth',lw)
    end
    xlabel('Time [hrs]')
    if I == 1
        ylabel('OD 600')
    end
end

% LB, -
Inds = [3,23,43];
c = [[0 0 0];orange];
for I = 1:2
    I1 = 1; if I == 1; I1=2; end
    subplot(3,6,16+I)
    hold on
    axis([0 96 0 ymax])
    set(gca,'YScale','log')
    set(gca,'fontsize',12)
    xticks(0:24:96)
    yticks([0.25,0.5,1,2])
    yticklabels({'2^{-2}','2^{-1}','2{^0}','2'});
    for ind2 = Inds+I-1
        plot(t_v,num(1:end,ind2),'color',c(I,:),'linewidth',lw)
    end
    xlabel('Time [hrs]')
    if I == 1
        ylabel('OD 600')
    end
end

set(gcf,'renderer','painters','position',[0 0 1400 800])
%saveas(gcf,'DM_ind.svg')


%% AO Means
figure(3)

lw = 1;
fs = 16;
fa = 0.3;

subplot(1,3,1)
hold on
axis([0 96 0 ymax])
set(gca,'YScale','log')
set(gca,'fontsize',12)
xticks(0:24:96)
yticks([0.25,0.5,1,2])
yticklabels({'2^{-2}','2^{-1}','2{^0}','2'});

% AO, -
Inds = [1,21,41];
c = [[0 0 0]];
ls = '-';
me = mean(num(1:end,Inds),2);
plot(t_v,me,'color',c,'linewidth',lw,'linestyle',ls)
se = std(num(1:end,Inds),'',2)./sqrt(length(Inds));
%plot(t_v,me+se,'color',c,'linewidth',lw2,'linestyle',ls);
%plot(t_v,me-se,'color',c,'linewidth',lw2,'linestyle',ls);
patch([t_v,flipud(t_v')'],[me+se;flipud(me-se)],c,'facealpha',fa,'edgealpha',0);

% AO, AO
%Inds = [37,47,57,38,48,58];
Inds = [37,57,38,58];%remove outlier interaction
c = [purple];
%c = [0.5,0.1,0.5];
ls = '--';
me = mean(num(1:end,Inds),2);
plot(t_v,me,'color',c,'linewidth',lw,'linestyle',ls)
se = std(num(1:end,Inds),'',2)./sqrt(length(Inds));
%plot(t_v,me+se,'color',c,'linewidth',lw2,'linestyle',ls);
%plot(t_v,me-se,'color',c,'linewidth',lw2,'linestyle',ls);
patch([t_v,flipud(t_v')'],[me+se;flipud(me-se)],c,'facealpha',fa,'edgealpha',0);

% AO, LP
Inds = [11,31,51,7];
c = [yellow];
%c = [0.9,0.2,0.9];
ls = '-';
me = mean(num(1:end,Inds),2);
plot(t_v,me,'color',c,'linewidth',lw,'linestyle',ls)
se = std(num(1:end,Inds),'',2)./sqrt(length(Inds));
%plot(t_v,me+se,'color',c,'linewidth',lw2,'linestyle',ls);
%plot(t_v,me-se,'color',c,'linewidth',lw2,'linestyle',ls);
patch([t_v,flipud(t_v')'],[me+se;flipud(me-se)],c,'facealpha',fa,'edgealpha',0);

% AO, LB
Inds = [13,33,53,17];
c = [orange];
%c = [0.9,0.2,0.9];
ls = '-';
me = mean(num(1:end,Inds),2);
plot(t_v,me,'color',c,'linewidth',lw,'linestyle',ls)
se = std(num(1:end,Inds),'',2)./sqrt(length(Inds));
%plot(t_v,me+se,'color',c,'linewidth',lw2,'linestyle',ls);
%plot(t_v,me-se,'color',c,'linewidth',lw2,'linestyle',ls);
patch([t_v,flipud(t_v')'],[me+se;flipud(me-se)],c,'facealpha',fa,'edgealpha',0);

xlabel('Time [hrs]')
ylabel('OD 600')
set(gca,'fontsize',fs)


%% LP Means
subplot(1,3,2)
hold on
axis([0 96 0 ymax])
set(gca,'YScale','log')
set(gca,'fontsize',12)
xticks(0:24:96)
yticks([0.25,0.5,1,2])
yticklabels({'2^{-2}','2^{-1}','2{^0}','2'});

% LP, -
Inds = [5,25,45];
c = [[0 0 0]];
ls = '-';
me = mean(num(1:end,Inds),2);
plot(t_v,me,'color',c,'linewidth',lw,'linestyle',ls)
se = std(num(1:end,Inds),'',2)./sqrt(length(Inds));
%plot(t_v,me+se,'color',c,'linewidth',lw2,'linestyle',ls);
%plot(t_v,me-se,'color',c,'linewidth',lw2,'linestyle',ls);
patch([t_v,flipud(t_v')'],[me+se;flipud(me-se)],c,'facealpha',fa,'edgealpha',0);

% LP, LP
Inds = [9,19,29,10,20,30];
c = [yellow];
%c = [0.5,0.5,0.1];
ls = '--';
me = mean(num(1:end,Inds),2);
plot(t_v,me,'color',c,'linewidth',lw,'linestyle',ls)
se = std(num(1:end,Inds),'',2)./sqrt(length(Inds));
%plot(t_v,me+se,'color',c,'linewidth',lw2,'linestyle',ls);
%plot(t_v,me-se,'color',c,'linewidth',lw2,'linestyle',ls);
patch([t_v,flipud(t_v')'],[me+se;flipud(me-se)],c,'facealpha',fa,'edgealpha',0);

% LP, AO
Inds = [12,32,52,8];
c = [purple];
%c = [0.9,0.9,0.2];
ls = '-';
me = mean(num(1:end,Inds),2);
plot(t_v,me,'color',c,'linewidth',lw,'linestyle',ls)
se = std(num(1:end,Inds),'',2)./sqrt(length(Inds));
%plot(t_v,me+se,'color',c,'linewidth',lw2,'linestyle',ls);
%plot(t_v,me-se,'color',c,'linewidth',lw2,'linestyle',ls);
patch([t_v,flipud(t_v')'],[me+se;flipud(me-se)],c,'facealpha',fa,'edgealpha',0);

% LP, LB
Inds = [15,35,55,27];
c = [orange];
%c = [0.9,0.9,0.2];
ls = '-';
me = mean(num(1:end,Inds),2);
plot(t_v,me,'color',c,'linewidth',lw,'linestyle',ls)
se = std(num(1:end,Inds),'',2)./sqrt(length(Inds));
%plot(t_v,me+se,'color',c,'linewidth',lw2,'linestyle',ls);
%plot(t_v,me-se,'color',c,'linewidth',lw2,'linestyle',ls);
patch([t_v,flipud(t_v')'],[me+se;flipud(me-se)],c,'facealpha',fa,'edgealpha',0);

xlabel('Time [hrs]')
ylabel('OD 600')
set(gca,'fontsize',fs)


%% LB Means
subplot(1,3,3)
hold on
axis([0 96 0 ymax])
set(gca,'YScale','log')
set(gca,'fontsize',12)
xticks(0:24:96)
yticks([0.25,0.5,1,2])
yticklabels({'2^{-2}','2^{-1}','2{^0}','2'});

% LB, -
Inds = [4,24,44];
c = [[0 0 0]];
ls = '-';
me = mean(num(1:end,Inds),2);
plot(t_v,me,'color',c,'linewidth',lw,'linestyle',ls)
se = std(num(1:end,Inds),'',2)./sqrt(length(Inds));
%plot(t_v,me+se,'color',c,'linewidth',lw2,'linestyle',ls);
%plot(t_v,me-se,'color',c,'linewidth',lw2,'linestyle',ls);
patch([t_v,flipud(t_v')'],[me+se;flipud(me-se)],c,'facealpha',fa,'edgealpha',0);

% LB, LB
Inds = [39,49,59,40,50,60];
c = [orange];
%c = [0.5,0.3,0.1];
ls = '--';
me = mean(num(1:end,Inds),2);
plot(t_v,me,'color',c,'linewidth',lw,'linestyle',ls)
se = std(num(1:end,Inds),'',2)./sqrt(length(Inds));
%plot(t_v,me+se,'color',c,'linewidth',lw2,'linestyle',ls);
%plot(t_v,me-se,'color',c,'linewidth',lw2,'linestyle',ls);
patch([t_v,flipud(t_v')'],[me+se;flipud(me-se)],c,'facealpha',fa,'edgealpha',0);

% LB, AO
Inds = [14,34,54,18];
c = [purple];
%c = [0.9,0.5,0.2];
ls = '-';
me = mean(num(1:end,Inds),2);
plot(t_v,me,'color',c,'linewidth',lw,'linestyle',ls)
se = std(num(1:end,Inds),'',2)./sqrt(length(Inds));
%plot(t_v,me+se,'color',c,'linewidth',lw2,'linestyle',ls);
%plot(t_v,me-se,'color',c,'linewidth',lw2,'linestyle',ls);
patch([t_v,flipud(t_v')'],[me+se;flipud(me-se)],c,'facealpha',fa,'edgealpha',0);

% LB, LP
Inds = [16,36,56,28];
c = [yellow];
%c = [0.9,0.5,0.2];
ls = '-';
me = mean(num(1:end,Inds),2);
plot(t_v,me,'color',c,'linewidth',lw,'linestyle',ls)
se = std(num(1:end,Inds),'',2)./sqrt(length(Inds));
%plot(t_v,me+se,'color',c,'linewidth',lw2,'linestyle',ls);
%plot(t_v,me-se,'color',c,'linewidth',lw2,'linestyle',ls);
patch([t_v,flipud(t_v')'],[me+se;flipud(me-se)],c,'facealpha',fa,'edgealpha',0);

xlabel('Time [hrs]')
ylabel('OD 600')
set(gca,'fontsize',fs)

set(gcf,'renderer','painters','position',[0 0 1400 300])

%saveas(gcf,'DM_Means.svg')

%% STATS
% Calculate significance of difference between growth yeilds
% KS test and Student T test

Inds = [1,21,41];AO_B = num(end,Inds);

%Inds = [37,47,57,38,48,58];
Inds = [37,57,38,58];%remove outlier interaction
AO_AO = num(end,Inds);

Inds = [11,31,51,7];AO_LP = num(end,Inds);

Inds = [13,33,53,17];AO_LB = num(end,Inds);

Inds = [5,25,45];LP_B = num(end,Inds);

Inds = [9,19,29,10,20,30];LP_LP = num(end,Inds);

Inds = [12,32,52,8];LP_AO = num(end,Inds);

Inds = [15,35,55,27];LP_LB = num(end,Inds);

Inds = [4,24,44];LB_B = num(end,Inds);

Inds = [39,49,59,40,50,60];LB_LB = num(end,Inds);

Inds = [14,34,54,18];LB_AO = num(end,Inds);

Inds = [16,36,56,28];LB_LP = num(end,Inds);


% AO_LB/AO_LB vs AO_AO/AO_B
x1 = [AO_LP,AO_LB];
x2 = [AO_B,AO_AO];
[~,p_KS_1]=kstest2(x1,x2,'Tail','smaller') %one sided KS test
[~,p_t_1]=ttest2(x1,x2,'Vartype','unequal','Tail','right') %one sided t-test, don't assume equal variance

% LP_AO/LP_B vs LP_LP/LP_LB
x1 = [LP_AO,LP_B]
x2 = [LP_LP,LP_LB]
[~,p_KS_2]=kstest2(x1,x2,'Tail','smaller') %one sided KS test
[~,p_t_2]=ttest2(x1,x2,'Vartype','unequal','Tail','right') %one sided t-test, don't assume equal variance

% LB_AO/LB_B vs LB_LB/LB_LP
x1 = [LB_AO,LB_B]
x2 = [LB_LB,LB_LP]
[~,p_KS_3]=kstest2(x1,x2,'Tail','smaller') %one sided KS test
[~,p_t_3]=ttest2(x1,x2,'Vartype','unequal','Tail','right') %one sided t-test, don't assume equal variance


% LP_AO vs LP_B
x1 = [LP_AO];
x2 = [LP_B];
[~,p_KS_4]=kstest2(x1,x2,'Tail','smaller') %one sided KS test
[~,p_t_4]=ttest2(x1,x2,'Vartype','unequal','Tail','right') %one sided t-test, don't assume equal variance

% LB_AO vs LB_B
x1 = [LB_AO];
x2 = [LB_B];
[~,p_KS_5]=kstest2(x1,x2,'Tail','smaller') %one sided KS test
[~,p_t_5]=ttest2(x1,x2,'Vartype','unequal','Tail','right') %one sided t-test, don't assume equal variance

