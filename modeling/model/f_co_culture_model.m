% David Bernstein
% 5/14/19
% Co-culture Auxotroph Equations

function dxdt = f_co_culture_model(t,x,p)
% Variables
% x1 = GL [mmol]
% x2 = GR [mmol]
% x3 = A1L [mmol]
% x4 = A1R [mmol]
% x5 = A2L [mmol]
% x6 = A2R [mmol]
% x7 = B1L [g]
% x8 = B1R [g]
% x9 = B2L [g]
% x10 = B2R [g]

% Parameters
% % Kinetics 
% p.vmaxG [mmol/(hr*g)]
% p.vmaxA1 [mmol/(hr*g)]
% p.vmaxA2 [mmol/(hr*g)]
% p.kG [mmol/L]
% p.kA1 [mmol/L]
% p.kA2 [mmol/L]
% % Biomass Stoichiometry 
% p.zG [g/mmol]
% p.zA1 [g/mmol]
% p.zA2 [g/mmol]
% % Secretion Stoichiometry
% p.yA1 [mmol/g]
% p.yA2 [mmol/g]
% % Volume
% p.v [L]
% % Diffusion
% p.d [L/hr]

%% NonNegative (commented because applied in ODE parameters)
% for I = 1:length(x)
%     if x(I) < 0
%         x(I)= 0;
%     end
% end

%% Left side
% Bio 1
% uptake bounds 1
uGL1 = (p.vmaxG*(x(1)/p.v))/(p.kG+(x(1)/p.v)); %[mmol/(hr*g)]
uA1L1 = (p.vmaxA1*(x(3)/p.v))/(p.kA1+(x(3)/p.v));
% biomass flux 1
vBL1 = min(uGL1*p.zG,uA1L1*p.zA1);
% uptake flux 1
vGL1 = vBL1/p.zG;
vA1L1 = vBL1/p.zA1;
% production flux 1
vA2L1 = vBL1*p.yA2;

% Bio 2
% uptake bounds
uGL2 = (p.vmaxG*(x(1)/p.v))/(p.kG+(x(1)/p.v));
uA2L2 = (p.vmaxA2*(x(5)/p.v))/(p.kA2+(x(5)/p.v));
% biomass flux
vBL2 = min(uGL2*p.zG,uA2L2*p.zA2);
% uptake flux
vGL2 = vBL2/p.zG;
vA2L2 = vBL2/p.zA2;
% production flux
vA1L2 = vBL2*p.yA1;

%% Right side
% Bio 2
% uptake bounds
uGR2 = (p.vmaxG*(x(2)/p.v))/(p.kG+(x(2)/p.v));
uA2R2 = (p.vmaxA2*(x(6)/p.v))/(p.kA2+(x(6)/p.v));
% biomass flux
vBR2 = min(uGR2*p.zG,uA2R2*p.zA2);
% uptake flux
vGR2 = vBR2/p.zG;
vA2R2 = vBR2/p.zA2;
% production flux
vA1R2 = vBR2*p.yA1;

% Bio 1
% uptake bounds 1
uGR1 = (p.vmaxG*(x(2)/p.v))/(p.kG+(x(2)/p.v)); %[mmol/(hr*g)]
uA1R1 = (p.vmaxA1*(x(4)/p.v))/(p.kA1+(x(4)/p.v));
% biomass flux 1
vBR1 = min(uGR1*p.zG,uA1R1*p.zA1);
% uptake flux 1
vGR1 = vBR1/p.zG;
vA1R1 = vBR1/p.zA1;
% production flux 1
vA2R1 = vBR1*p.yA2;

%% Diffusion
dG = p.d*((x(1)/p.v)-(x(2)/p.v));
dA1 = p.d*((x(3)/p.v)-(x(4)/p.v));
dA2 = p.d*((x(5)/p.v)-(x(6)/p.v));

%% Dynamics
dxdt = zeros(length(x),1);
dxdt(1) = -vGL1*x(7)-vGL2*x(9)-dG; %GL
dxdt(2) = -vGR2*x(10)-vGR1*x(8)+dG; %GR
dxdt(3) = -vA1L1*x(7)+vA1L2*x(9)-dA1; %A1L
dxdt(4) = vA1R2*x(10)-vA1R1*x(8)+dA1; %A1R
dxdt(5) = vA2L1*x(7)-vA2L2*x(9)-dA2; %A2L
dxdt(6) = -vA2R2*x(10)+vA2R1*x(8)+dA2; %A2R
dxdt(7) = vBL1*x(7); %B1L
dxdt(8) = vBR1*x(8); %B1R
dxdt(9) = vBL2*x(9);%B2L
dxdt(10) = vBR2*x(10);%B2R

end



