clear
close all
clc
addpath(genpath('.\Functions'))

%% Initialisation

%Parameters
params = struct( ...
    'ctrlMax_dim', 1e-7, ...
    'md_lim_dim',  1, ...
    'pocLim', 1e-6, ... 
    'nx_orb', 60, ...
    'n_orb', 1 ...
    );

% Scenario definition
[primary,secondary] = generateInitShort(1);
scenario            = nondimensionalise(primary,secondary,params);

% Write json input
input = write_input(scenario, params);

%% Optimisation
!wsl ./build/bin/backSweep

[control, rB, tca, md] = readBackSweepOutput('output.json');

%% Validation and postprocessing
[rB_val, miss_dist_val, smd, poc] = validateBackSweep('output.json');
rB_val(rB_val==0) = nan;
rB(rB==0)         = nan;
% figure
% plot(tca.shift_s*Tsc)
% hold on
% plot(dtca*Tsc)
% hold off

figure
plot(control)
legend('R','T','N')

figure
plot(sqrt(md))
hold on
plot(miss_dist_val*input.Lsc)
hold off

figure
semilogy(poc)
ylim([1e-10,1e-2])


figure
plot(smd)
hold on
plot(scenario.smdLim)

showEllipseBplane(scenario.Pb,input.lim,rB,rB_val,input.metric_case,input.Lsc);