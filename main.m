clear
close all
clc
addpath(genpath('.\Functions'))

%% Initialisation

%Parameters
params = struct( ...
    'ctrlMax_dim', 1e-7, ...
    'md_lim_dim',  0.5, ...
    'pocLim', 1e-6, ... 
    'nx_orb', 30, ...
    'n_orb', 1, ...
    'breakOnThreshold', 1, ...
    'metric_case', 1 ...
    );

% Scenario definition
[primary,secondary] = generateInitShort(1);
scenario            = nondimensionalise(primary,secondary,params);

% Write json input
input = write_input(scenario, params);

%% Optimisation
!wsl ./build/bin/backSweep

%%
outSim = readBackSweepOutput('output.json',params.metric_case);

%% Validation and postprocessing
validation = validateBackSweep('output.json');

mainPostprocess(outSim, validation, scenario, input, params)