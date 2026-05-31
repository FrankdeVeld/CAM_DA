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
    'mass',               800, ...  % kg (Starlink v2)
    'thrust',             300, ...  % mN (Starlink v2)
    'md_lim_dim',         2, ...    % km MD limit
    'pocLim',             1e-6, ... % - PoC limit
    'nx_orb',             120, ...  % Number of nodes per orbit
    'n_orb',              0.7, ...    % Number of orbits before TCA
    'n_orb_start',        0, ...    % Number of orbits before TCA (redundant)
    'metric_case',        1, ...    % 1: miss distance, 2: SMD
    'order',              2, ...    % DA order (between 2 and 6)
    'tCAHandling',        2, ...    % 0: Fix (not working), 1: DA, 2: Picard-Lindelhof
    'refineLastInterval', 1 ...     % 0: no, 1: yes
    );
params.ctrlMax_dim = params.thrust/1e6/params.mass; % km/s^2 control magnitude 

% Scenario definition
% [primary,secondary] = generateInitShort(184);
[primary,secondary] = readCDM_CONGEN();
% Nondimensionalisation of the initial conditions and parameters
scenario            = nondimensionalise(primary,secondary,params);

% Write json input for C++ code
input = write_input(scenario, params);

%% Optimisation
!wsl ./build/bin/backSweep

%Read json output
outSim = readBackSweepOutput('output.json',input,scenario,params.metric_case);

%% Validation
validation = validateBackSweep(input,outSim);

%% Postprocessing
mainPostprocess(outSim, validation, scenario, input, params)