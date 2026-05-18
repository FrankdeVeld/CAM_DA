clear
close all
clc
addpath(genpath('.\Functions'))

%% Initialisation
case_ids = 1:2170;
numCases = numel(case_ids);

%Parameters
params = struct( ...
    'mass',             800, ... % kg (Starlink v2)
    'thrust',           300, ... % mN (Starlink v2)
    'md_lim_dim',       2, ...   % km
    'pocLim',           1e-6, ... 
    'nx_orb',           60, ...  % Number of nodes per orbit
    'n_orb',            0.5, ...   % Number of orbits before TCA
    'metric_case',      2, ...    % 1: miss distance, 2: SMD
    'order',            2, ...    % 1: DA order
    'tCAHandling',      2 ...     % 1: Fix, 2: changing
    );
params.ctrlMax_dim = params.thrust/1e6/params.mass; % km/s^2

%% Preallocation
scenario_id      = nan(numCases,1);
dvTot     = nan(numCases,1);
abs_err = nan(numCases,1);
rel_err = nan(numCases,1);
max_tca_shift    = nan(numCases,1);
solver_success   = false(numCases,1);
error_message    = strings(numCases,1);

%% Optional detailed storage
outSim_all(numCases,1)     = struct();
validation_all(numCases,1) = struct();
scenario_all(numCases,1)   = struct();
input_all(numCases,1)      = struct();

%% Loop over scenarios
for kk = 1:numCases
    icase = case_ids(kk);
    scenario_id(kk) = icase;
    fprintf('\n============================================================\n');
    fprintf('Running scenario %d / %d (ID = %d)\n', kk, numCases, icase);
    fprintf('============================================================\n');
    
    try
        % Initialisation
        [primary, secondary] = generateInitShort(icase);
        scenario = nondimensionalise(primary, secondary, params);
        input    = write_input(scenario, params);

        % Optimisation
        !wsl ./build/bin/backSweep
        
        % Read optimisation output
        outSim = readBackSweepOutput('output.json',input,scenario,params.metric_case);

        % Validation
        validation = validateBackSweep(input,outSim);

        %% Requested metrics
        % 1) total deltaV
        % Assume outSim.control is N x 3 and scenario.dt or input.dt may exist.
        % Fallback: reconstruct dt from time history if needed.
        u_hist    = outSim.control;
        t         = outSim.t_nodes*scenario.Tsc;
        dvTot(kk) = normOfVec(outSim.control(1:end-1,:)')*...
                 diff(-t)'*scenario.Asc*scenario.ctrlMax*1000; % m/s
        
        % 2) maximum validation error on SMD or miss distance
        % Prefer explicit error histories if validateBackSweep returns them.
        if params.metric_case == 1
            md_abs_err = sqrt(outSim.m_d)-validation.m_d;
            md_err     = 100*abs(1-sqrt(outSim.m_d)./validation.m_d);
            abs_err(kk) = max(abs(md_abs_err));
            rel_err(kk) = max(md_err);
        else
            poc_abs_err = outSim.poc-validation.poc;
            smd_abs_err = outSim.smd-validation.smd;
            smd_err     = 100*abs(1-outSim.smd./validation.smd);
            abs_err(kk) = max(abs(poc_abs_err));
            rel_err(kk) = max(smd_err);
        end

        % 3) maximum TCA shift
        max_tca_shift(kk) = max(abs(outSim.deltaTca))*scenario.Tsc;
       

        %% Save detailed structs if desired
        solver_success(kk) = true;

    catch ME
        solver_success(kk) = false;
        error_message(kk) = string(ME.message);
        warning('Scenario %d failed: %s', icase, ME.message);
    end
end

%% Summary table
results = table(scenario_id, solver_success, dvTot, abs_err, rel_err, max_tca_shift, error_message, ...
    'VariableNames', {'scenario_id','solver_success','dvTot [m/s]','max_abs_error','max_rel_error','max_tca_shift [s]','error_message'});

%% Save results
save('batch_results.mat', 'results', 'outSim_all', 'validation_all', 'scenario_all', 'input_all', 'params');