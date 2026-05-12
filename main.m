clear
close all
clc
addpath(genpath('.\OrbitalDynamics'))
addpath(genpath('.\PoCIntegrals'))

%% Initialisation
[primary,secondary] = generateInitShort(243);

% Non-dimensionalization
Lsc   = primary.a;
musc  = 398600.4418;
Vsc   = sqrt(musc/Lsc);
Tsc   = Lsc/Vsc;
Asc   = Vsc/Tsc;
scale = diag([ones(3,1)./Lsc; ones(3,1)./Vsc]);

x_p     = scale*primary.x0;
x_s     = scale*secondary.x0;
r2ep    = rtn2eci(x_p(1:3),x_p(4:6));
r2es    = rtn2eci(x_s(1:3),x_s(4:6));
cov     = scale(1:3,1:3)*(r2ep*primary.C0(1:3,1:3)*r2ep'+r2es*secondary.C0(1:3,1:3)*r2es')*scale(1:3,1:3);
e2b     = eci2Bplane(primary.x0(4:6),secondary.x0(4:6)); 
e2b     = e2b([1 3],:);
Pb      = e2b*cov*e2b';
ctrlMax = 1e-7/Asc;
T       = primary.T/Tsc;
HBR     = (primary.HBR + secondary.HBR)/Lsc;
md_lim  = 1; %km
pocLim  = 1e-6; 
smdLim  = PoC2SMD(Pb,HBR,pocLim, 5, 1, 1e-3, 200); 
nx_orb  = 60;
n_orb   = .5;

% Write json input
input         = struct();
input.N       = nx_orb*n_orb+1; 
input.Lsc     = Lsc; 
input.et      = 659871.07119168108; 
input.t_back  = n_orb*2*pi; 
input.uMax    = ctrlMax;
input.scaling = Lsc;
input.xp_tCA  = x_p';
input.xs_tCA  = x_s';
input.P       = cov;
input.HBR     = HBR;
input.lim     = (md_lim/Lsc)^2;
input.metric_case = 2;
input.tCAHandling = 2;

fid = fopen('./input.json','w'); 
fwrite(fid,jsonencode(input),'char'); 
fclose(fid);

%% Optimisation
!wsl ./build/bin/backSweep

[control, rB, tca, md] = readBackSweepOutput('output.json');

%% Validation and postprocessing
[rB_val, miss_dist_val, smd, poc] = validateBackSweep('output.json');

figure
plot(tca.shift_s*Tsc)
% hold on
% plot(dtca*Tsc)
% hold off

figure
plot(control)
legend('R','T','N')

figure
plot(sqrt(md))
hold on
plot(miss_dist_val*Lsc)
hold off

figure
semilogy(poc)


figure
plot(smd)
hold on
plot(smdLim)

showEllipseBplane(Pb,smdLim,rB,rB_val,input.metric_case,Lsc);