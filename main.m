clear
% close all
clc
addpath(genpath('.\Functions'))

%% Initialisation

%Parameters
params = struct( ...
    'mass',             800, ... % kg (Starlink v2)
    'thrust',           300, ... % mN (Starlink v2)
    'md_lim_dim',       2, ...   % km
    'pocLim',           1e-6, ... 
    'nx_orb',           360, ...  % Number of nodes per orbit
    'n_orb',            0.2, ...   % Number of orbits before TCA
    'breakOnThreshold', 1, ...
    'metric_case',      2 ...    % 1: miss distance, 2: SMD
    );
params.ctrlMax_dim = params.thrust/1e6/params.mass; % km/s^2

% Scenario definition
[primary,secondary] = generateInitShort(1);
scenario            = nondimensionalise(primary,secondary,params);

% Write json input
input = write_input(scenario, params);

%% Optimisation
!wsl ./build/bin/backSweep

%%
outSim = readBackSweepOutput('output.json',input,scenario,params.metric_case);

%% Validation and postprocessing
validation = validateBackSweep('output.json',scenario);

mainPostprocess(outSim, validation, scenario, input, params)