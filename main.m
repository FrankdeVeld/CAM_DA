clear
close all
clc
addpath(genpath('.\Functions'))
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultAxesFontName','Times');
set(0,'DefaultUicontrolFontName','Times', 'DefaultUicontrolFontSize', 16);
set(0,'DefaultUitableFontName','Times', 'DefaultUitableFontSize', 16);
set(0,'DefaultTextFontName','Times', 'DefaultTextFontSize', 16);
set(0,'DefaultUipanelFontName','Times', 'DefaultUipanelFontSize', 16);

set(0, 'DefaultLineLineWidth', 1);
set(0,'defaultfigurecolor',[1 1 1])
    
%% Initialisation

%Parameters
params = struct( ...
    'mass',             800, ... % kg (Starlink v2)
    'thrust',           300, ... % mN (Starlink v2)
    'md_lim_dim',       2, ...   % km
    'pocLim',           1e-6, ... 
    'nx_orb',           360, ...  % Number of nodes per orbit
    'n_orb',            0.2, ...   % Number of orbits before TCA
    'metric_case',      2, ...    % 1: miss distance, 2: SMD
    'order',            6, ...    % 1: DA order
    'tCAHandling',      2 ...     % 1: Fix, 2: changing
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
validation = validateBackSweep(input,outSim);

mainPostprocess(outSim, validation, scenario, input, params)