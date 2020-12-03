% David Bernstein
% 6/27/19
% Co-culture Modeling
% Simulate Experiments

clear all
close all
clc

fntsz = 12;

% Add path to model
addpath('C:\Users\David\GitHub\co_culture_device\modeling\model')

%% Define Parameters
% Kinetics (Harcombe et al. cell systems. 2014; Gosset, Microb. Cell Fact. 2005)
p.vmaxG = 10; % 10 %[mmol/(hr*g)]
p.vmaxA1 = 10; % 10 %[mmol/(hr*g)]
p.vmaxA2 = 10; % 10 %[mmol/(hr*g)]
p.kG = 0.01; % 0.01 %[mmol/L], 1.75e-3
p.kA1 = 0.01; % 0.01 %[mmol/L]
p.kA2 = 0.01; % 0.01 %[mmol/L]
% Biomass Stoichiometry 
% grams e coli / grams glucose = 0.5/1 (BioNumbers 105318, Shiloach et al. Biotechnol Adv. 2005)
p.zG = 0.0901; % 0.5*180.156e-3=0.0901 %[g/mmol]
% cell / mol amino acid (Mee et al. PNAS. 2014)
% gram / cell = 500*10^-15 (cell biology by the numbers)
p.zA1 = 2.7374; %(5.4747e15)*(500e-15)*10^-3=2.7374 %[g/mmol] = [cell/mol]*[g/cell]*[mol/mmol] Lysine
p.zA2 = 4.0148; %(8.0295e15)*(500e-15)*10^-3=4.0148 %[g/mmol] = [cell/mol]*[g/cell]*[mol/mmol] Isoleucine
% Secretion Stoichiometry
p.yA1 = 10^-0.4; %[mmol/g]
p.yA2 = p.yA1; %[mmol/g]
% Volume (From experiment)
p.v = 250*10^-6; %[L]
% Diffusion
p.d = 2.39*10^-5; %[L/hr]

%% Experiment1: Simulate Chemical Diffusion
% Initial Conditions
IG = 5.55e-3;
initial_G = IG; %[mmol]
initial_A1 = 0; %[mmol]
initial_A2 = 0; %[mmol]
initial_Bio1 = 0; %[g]
initial_Bio2 = 0; %[g]
x0 = [initial_G,0,... % Glucose [L,R]
    initial_A1,initial_A1,... % AA1 [L,R]
    initial_A2,initial_A2,... % AA2 [L,R]
    initial_Bio1,initial_Bio1,... % Biomass1 [L,R]
    initial_Bio2,initial_Bio2]; % Biomass2 [L,R]

% ODE
tspan = [0:0.25:96];
options = odeset('NonNegative',[1:10],'MaxStep',0.01);
[t,x] = ode23tb(@(t,x) f_co_culture_model(t,x,p), tspan, x0, options);
% Record Results
t1 = t;
x1 = x;

% Plot
yellow = [0.9 0.6 0];
figure(1) %Glucose
subplot(1,2,1) %left
plot(t,x(:,1),'-','linewidth',2,'color',yellow)
xlabel('time [hrs]')
ylabel('Glucose (left) [mmol]')
axis([0 t(end) 0 initial_G])
set(gca,'fontsize',fntsz)
subplot(1,2,2) %right
plot(t,x(:,2),'-','linewidth',2,'color',yellow)
xlabel('time [hrs]')
ylabel('Glucose (right) [mmol]')
axis([0 t(end) 0 initial_G])
set(gca,'fontsize',fntsz)
set(gcf,'position',[100,100,500,200])
xticks(0:24:96)

set(gcf,'renderer','painters')
%saveas(gcf,'B.svg')

%% Experiment 2: Amino Acid Growth
% A1, B1
% Initial Conditions
initial_G = IG; %[mmol] % Glucose [mmol]: 20% stock (20g/100mL) / MW * Conc. (2mL/100mL) * Vol (250uL)
initial_A1 = 1.83e-4; %[mmol] % AA1 [mmol]: aa/cell (1.1e8) * cell yield (10^9) / N_A (6.022e23) = 1.83e-4 mmol
initial_A2 = 0; %[mmol] % AA2 [mmol]: aa/cell (7.5e7) * cell yield (10^9) / N_A (6.022e23) = 1.25e-4 mmol
initial_Bio1 = 6e-8; %[g]
initial_Bio2 = 0; %[g] Biomass [g]: (OD600 Culture (0.1) - OD M9 (0.04)) * OD600 E Coli Calib (8.0e8 cells/mL)  * Vol (0.250mL) * Dilution Factor (1/100) * Mass per cell (500e-15 g/cell) = 6.0e-8 g
x0 = [initial_G,initial_G,... % Glucose [L,R]
    0,initial_A1,... % AA1 [L,R]
    initial_A2,initial_A2,... % AA2 [L,R]
    initial_Bio1,0,... % Biomass1 [L,R]
    initial_Bio2,initial_Bio2]; % Biomass2 [L,R]

% ODE
%tspan = [0:0.01:48];
options = odeset('NonNegative',[1:10],'MaxStep',0.01);
[t,x] = ode23tb(@(t,x) f_co_culture_model(t,x,p), tspan, x0, options);
% Record Results
t2_1 = t;
x2_1 = x;

% Plot
green1 = [0 0.7 0];
blue1 = [0 0 0.7];
green2 = [0.4 1 0.4];
blue2 =  [0.4 0.4 1];
figure(2)
subplot(2,2,1) %AA1 Left
plot(t,x(:,3),'-','linewidth',2,'color',blue2)
xlabel('time [hrs]')
ylabel('Lysine (left) [mmol]')
axis([0 t(end) 0 initial_A1])
set(gca,'fontsize',fntsz)
subplot(2,2,2) %AA1 Right
plot(t,x(:,4),'-','linewidth',2,'color',blue2)
xlabel('time [hrs]')
ylabel('Lysine (right) [mmol]')
axis([0 t(end) 0 initial_A1])
set(gca,'fontsize',fntsz)
subplot(2,2,3) %Bio1 Left
plot(t,x(:,7),'-','linewidth',2,'color',blue1)
xlabel('time [hrs]')
ylabel('Bio \DeltaK (left) [g]')
axis([0 t(end) 0 max(x(:,7))])
set(gca,'fontsize',fntsz)
%set(gca,'YScale','log')
xticks(0:24:96)
%yticks([1.25*10^-4,2.5*10^-4,5*10^-4])
subplot(2,2,4) %Bio1 Right
plot(t,x(:,8),'-','linewidth',2,'color',blue1)
xlabel('time [hrs]')
ylabel('Bio \DeltaK (right) [g]')
axis([0 t(end) 0 max(x(:,7))])
set(gca,'fontsize',fntsz)
%set(gca,'YScale','log')
xticks(0:24:96)
%yticks([1.25*10^-4,2.5*10^-4,5*10^-4])
set(gcf,'position',[100,100,500,400])

set(gcf,'renderer','painters')
%saveas(gcf,'C.svg')

% A2, B2
% Initial Conditions
initial_G = IG; %[mmol] % Glucose [mmol]: 20% stock (20g/100mL) / MW * Conc. (2mL/100mL) * Vol (250uL)
initial_A1 = 0; %[mmol] % AA1 [mmol]: aa/cell (1.1e8) * cell yield (10^9) / N_A (6.022e23) = 1.83e-4 mmol
initial_A2 = 1.25e-4; %[mmol] % AA2 [mmol]: aa/cell (7.5e7) * cell yield (10^9) / N_A (6.022e23) = 1.25e-4 mmol
initial_Bio1 = 0; %[g]
initial_Bio2 = 6e-8; %[g] Biomass [g]: (OD600 Culture (0.1) - OD M9 (0.04)) * OD600 E Coli Calib (8.0e8 cells/mL)  * Vol (0.250mL) * Dilution Factor (1/100) * Mass per cell (500e-15 g/cell) = 6.0e-8 g
x0 = [initial_G,initial_G,... % Glucose [L,R]
    initial_A1,initial_A1,... % AA1 [L,R]
    0,initial_A2,... % AA2 [L,R]
    initial_Bio1,initial_Bio1,... % Biomass1 [L,R]
    initial_Bio2,0]; % Biomass2 [L,R]

% ODE
%tspan = [0:0.01:48];
options = odeset('NonNegative',[1:10],'MaxStep',0.01);
[t,x] = ode23tb(@(t,x) f_co_culture_model(t,x,p), tspan, x0, options);
% Record Results
t2_2 = t;
x2_2 = x;

% Plot
figure(3)
subplot(2,2,1) %AA2 Left
plot(t,x(:,5),'-','linewidth',2,'color',green2)
xlabel('time [hrs]')
ylabel('Isoleucine (left) [mmol]')
axis([0 t(end) 0 initial_A2])
set(gca,'fontsize',fntsz)
subplot(2,2,2) %AA2 Right
plot(t,x(:,6),'-','linewidth',2,'color',green2)
xlabel('time [hrs]')
ylabel('Isoleucine (right) [mmol]')
axis([0 t(end) 0 initial_A2])
set(gca,'fontsize',fntsz)
subplot(2,2,3) %Bio2 Left
plot(t,x(:,9),'-','linewidth',2,'color',green1)
xlabel('time [hrs]')
ylabel('Bio \DeltaI (left) [g]')
axis([0 t(end) 0 max(x(:,9))])
set(gca,'fontsize',fntsz)
subplot(2,2,4) %Bio2 Right
plot(t,x(:,10),'-','linewidth',2,'color',green1)
xlabel('time [hrs]')
ylabel('Bio \DeltaI (right) [g]')
axis([0 t(end) 0 max(x(:,9))])
set(gca,'fontsize',fntsz)
set(gcf,'position',[100,100,500,400])
xticks(0:24:96)

set(gcf,'renderer','painters')
%saveas(gcf,'D.svg')

%% Experiment 3: Co-culture
% opposite wells
% Initial Conditions
initial_G = IG; %[mmol] % Glucose [mmol]: 20% stock (20g/100mL) / MW * Conc. (2mL/100mL) * Vol (250uL)
initial_A1 = 1e-9; %[mmol] % AA1 [mmol]: aa/cell (1.1e8) * cell yield (10^9) / N_A (6.022e23) = 1.83e-4 mmol
initial_A2 = 1e-9; %[mmol] % AA2 [mmol]: aa/cell (7.5e7) * cell yield (10^9) / N_A (6.022e23) = 1.25e-4 mmol
initial_Bio1 = 6e-8; %[g]
initial_Bio2 = 6e-8; %[g] Biomass [g]: (OD600 Culture (0.1) - OD M9 (0.04)) * OD600 E Coli Calib (8.0e8 cells/mL)  * Vol (0.250mL) * Dilution Factor (1/100) * Mass per cell (500e-15 g/cell) = 6.0e-8 g
x0 = [initial_G,initial_G,... % Glucose [L,R]
    initial_A1,initial_A1,... % AA1 [L,R]
    initial_A2,initial_A2,... % AA2 [L,R]
    initial_Bio1,0,... % Biomass1 [L,R]
    0,initial_Bio2]; % Biomass2 [L,R]

% ODE
%tspan = [0:0.01:48];
options = odeset('NonNegative',[1:10],'MaxStep',0.01);
[t,x] = ode23tb(@(t,x) f_co_culture_model(t,x,p), tspan, x0, options);
% Record Results
t3_1 = t;
x3_1 = x;

% Plot
figure(4)
subplot(2,2,1) % Bio1 L
plot(t,x(:,7),'-','linewidth',2,'color',blue1)
xlabel('time [hrs]')
ylabel('Bio \DeltaK (left) [g]')
set(gca,'fontsize',fntsz)
%axis([0 t(end) 0 max([x(:,7);x(:,8);x(:,9);x(:,10)])])
subplot(2,2,2) % Bio1 R
plot(t,x(:,8),'-','linewidth',2,'color',blue1)
xlabel('time [hrs]')
ylabel('Bio \DeltaK (right) [g]')
%axis([0 t(end) 0 max([x(:,7);x(:,8);x(:,9);x(:,10)])])
set(gca,'fontsize',fntsz)
subplot(2,2,3) % Bio2 L
plot(t,x(:,9),'-','linewidth',2,'color',green1)
xlabel('time [hrs]')
ylabel('Bio \DeltaI (left) [g]')
%axis([0 t(end) 0 max([x(:,7);x(:,8);x(:,9);x(:,10)])])
set(gca,'fontsize',fntsz)
subplot(2,2,4) % Bio2 R
plot(t,x(:,10),'-','linewidth',2,'color',green1)
xlabel('time [hrs]')
ylabel('Bio \DeltaI (right) [g]')
%axis([0 t(end) 0 max([x(:,7);x(:,8);x(:,9);x(:,10)])])
set(gca,'fontsize',fntsz)
set(gcf,'position',[100,100,500,400])

% same wells
% Initial Conditions
initial_G = IG; %[mmol] % Glucose [mmol]: 20% stock (20g/100mL) / MW * Conc. (2mL/100mL) * Vol (250uL)
initial_A1 = 1e-9; %[mmol] % AA1 [mmol]: aa/cell (1.1e8) * cell yield (10^9) / N_A (6.022e23) = 1.83e-4 mmol
initial_A2 = 1e-9; %[mmol] % AA2 [mmol]: aa/cell (7.5e7) * cell yield (10^9) / N_A (6.022e23) = 1.25e-4 mmol
initial_Bio1 = 6e-8; %[g]
initial_Bio2 = 6e-8; %[g] Biomass [g]: (OD600 Culture (0.1) - OD M9 (0.04)) * OD600 E Coli Calib (8.0e8 cells/mL)  * Vol (0.250mL) * Dilution Factor (1/100) * Mass per cell (500e-15 g/cell) = 6.0e-8 g
x0 = [initial_G,initial_G,... % Glucose [L,R]
    initial_A1,initial_A1,... % AA1 [L,R]
    initial_A2,initial_A2,... % AA2 [L,R]
    initial_Bio1,initial_Bio1,... % Biomass1 [L,R]
    initial_Bio2,initial_Bio2]; % Biomass2 [L,R]

% ODE
%tspan = [0:0.01:48];
options = odeset('NonNegative',[1:10],'MaxStep',0.01);
[t,x] = ode23tb(@(t,x) f_co_culture_model(t,x,p), tspan, x0, options);
% Record Results
t3_2 = t;
x3_2 = x;

% Plot
figure(5)
subplot(2,2,1) % Bio1 L
plot(t,x(:,7),'-','linewidth',2,'color',blue1)
xlabel('time [hrs]')
ylabel('Bio \DeltaK (left) [g]')
%axis([0 t(end) 0 max([x(:,7);x(:,8);x(:,9);x(:,10)])])
set(gca,'fontsize',fntsz)
subplot(2,2,2) % Bio1 R
plot(t,x(:,8),'-','linewidth',2,'color',blue1)
xlabel('time [hrs]')
ylabel('Bio \DeltaK (right) [g]')
%axis([0 t(end) 0 max([x(:,7);x(:,8);x(:,9);x(:,10)])])
set(gca,'fontsize',fntsz)
subplot(2,2,3) % Bio2 L
plot(t,x(:,9),'-','linewidth',2,'color',green1)
xlabel('time [hrs]')
ylabel('Bio \DeltaI (left) [g]')
%axis([0 t(end) 0 max([x(:,7);x(:,8);x(:,9);x(:,10)])])
set(gca,'fontsize',fntsz)
subplot(2,2,4) % Bio2 R
plot(t,x(:,10),'-','linewidth',2,'color',green1)
xlabel('time [hrs]')
ylabel('Bio \DeltaI (right) [g]')
%axis([0 t(end) 0 max([x(:,7);x(:,8);x(:,9);x(:,10)])])
set(gca,'fontsize',fntsz)
set(gcf,'position',[100,100,500,400])

% Rescale figure 4 and 5
figure(4)
subplot(2,2,1) % Bio1 L
axis([0 t(end) 0 3*10^-4])
%set(gca,'YScale','log')
subplot(2,2,2) % Bio1 L
axis([0 t(end) 0 3*10^-4])
%set(gca,'YScale','log')
subplot(2,2,3) % Bio1 L
axis([0 t(end) 0 3*10^-4])
%set(gca,'YScale','log')
subplot(2,2,4) % Bio1 L
axis([0 t(end) 0 3*10^-4])
%set(gca,'YScale','log')
set(gcf,'renderer','painters')
%saveas(gcf,'E.svg')

figure(5)
subplot(2,2,1) % Bio1 L
axis([0 t(end) 0 3*10^-4])
%set(gca,'YScale','log')
subplot(2,2,2) % Bio1 L
axis([0 t(end) 0 3*10^-4])
%set(gca,'YScale','log')
subplot(2,2,3) % Bio1 L
axis([0 t(end) 0 3*10^-4])
%set(gca,'YScale','log')
subplot(2,2,4) % Bio1 L
axis([0 t(end) 0 3*10^-4])
%set(gca,'YScale','log')
set(gcf,'renderer','painters')
%saveas(gcf,'F.svg')

