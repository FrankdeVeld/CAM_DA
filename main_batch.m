clear
close all
clc
addpath(genpath('.\Functions'))

%% Batch settings
case_ids = 1:2;
numCases = numel(case_ids);

params = struct( ...
    'ctrlMax_dim', 1e-8, ...
    'md_lim_dim',  0.5, ...
    'pocLim', 1e-6, ... 
    'nx_orb', 60, ...
    'n_orb', 1, ...
    'breakOnThreshold', 1, ...
    'metric_case', 2 ...
    );

%% Preallocation
scenario_id      = nan(numCases,1);
total_deltaV     = nan(numCases,1);
max_metric_error = nan(numCases,1);
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
        %% Initialisation
        [primary, secondary] = generateInitShort(icase);
        scenario = nondimensionalise(primary, secondary, params);
        input    = write_input(scenario, params);

        %% Optimisation
        status = system('wsl ./build/bin/backSweep');
        if status ~= 0
            error('Optimizer returned non-zero exit code: %d', status);
        end

        %% Read optimisation output
        outSim = readBackSweepOutput('output.json', params.metric_case);

        %% Validation
        validation = validateBackSweep('output.json', scenario);

        %% Requested metrics
        % 1) total deltaV
        % Assume outSim.control is N x 3 and scenario.dt or input.dt may exist.
        % Fallback: reconstruct dt from time history if needed.
        if isfield(outSim, 'control')
            u_hist = outSim.control;
        elseif isfield(outSim, 'u')
            u_hist = outSim.u;
        else
            error('Could not find control history in outSim.');
        end

        if isfield(outSim, 't')
            t_hist = outSim.t(:);
            if numel(t_hist) > 1
                dt_hist = abs(diff(t_hist));
                if size(u_hist,1) == numel(dt_hist)
                    total_deltaV(kk) = sum(vecnorm(u_hist,2,2) .* dt_hist);
                elseif size(u_hist,1) == numel(t_hist)
                    total_deltaV(kk) = sum(vecnorm(u_hist(1:end-1,:),2,2) .* dt_hist);
                else
                    error('Size mismatch between control history and time history in outSim.');
                end
            else
                error('Time history in outSim is too short to compute deltaV.');
            end
        elseif isfield(scenario, 'dt')
            total_deltaV(kk) = sum(vecnorm(u_hist,2,2)) * scenario.dt;
        elseif isfield(input, 'dt')
            total_deltaV(kk) = sum(vecnorm(u_hist,2,2)) * input.dt;
        elseif isfield(input, 't_back') && isfield(input, 'N')
            total_deltaV(kk) = sum(vecnorm(u_hist,2,2)) * (input.t_back / input.N);
        else
            error('Could not infer time step for deltaV computation.');
        end

        % 2) maximum validation error on SMD or miss distance
        % Prefer explicit error histories if validateBackSweep returns them.
        if params.metric_case == 1
            if isfield(validation, 'md_error')
                max_metric_error(kk) = max(abs(validation.md_error));
            elseif isfield(validation, 'miss_distance_error')
                max_metric_error(kk) = max(abs(validation.miss_distance_error));
            elseif isfield(validation, 'md_hist') && isfield(outSim, 'metric')
                max_metric_error(kk) = max(abs(validation.md_hist(:) - outSim.metric(:)));
            else
                error('Could not find miss-distance validation error/history.');
            end
        elseif params.metric_case == 2
            if isfield(validation, 'smd_error')
                max_metric_error(kk) = max(abs(validation.smd_error));
            elseif isfield(validation, 'smd_hist') && isfield(outSim, 'metric')
                max_metric_error(kk) = max(abs(validation.smd_hist(:) - outSim.metric(:)));
            else
                error('Could not find SMD validation error/history.');
            end
        else
            if isfield(validation, 'metric_error')
                max_metric_error(kk) = max(abs(validation.metric_error));
            else
                error('Metric case 3 selected, but no generic validation.metric_error found.');
            end
        end

        % 3) maximum TCA shift
        if isfield(validation, 'deltaTca')
            max_tca_shift(kk) = max(abs(validation.deltaTca));
        elseif isfield(validation, 'tca_shift')
            max_tca_shift(kk) = max(abs(validation.tca_shift));
        else
            error('Could not find TCA shift history in validation output.');
        end

        %% Save detailed structs if desired
        outSim_all(kk)     = outSim;
        validation_all(kk) = validation;
        scenario_all(kk)   = scenario;
        input_all(kk)      = input;

        solver_success(kk) = true;
        fprintf('Scenario %d completed. dV = %.6e, max err = %.6e, max TCA shift = %.6e\n', ...
            icase, total_deltaV(kk), max_metric_error(kk), max_tca_shift(kk));

    catch ME
        solver_success(kk) = false;
        error_message(kk) = string(ME.message);
        warning('Scenario %d failed: %s', icase, ME.message);
    end
end

%% Summary table
results = table(scenario_id, solver_success, total_deltaV, max_metric_error, max_tca_shift, error_message, ...
    'VariableNames', {'scenario_id','solver_success','total_deltaV','max_metric_error','max_tca_shift','error_message'});

%% Save results
save('batch_results.mat', 'results', 'outSim_all', 'validation_all', 'scenario_all', 'input_all', 'params');
writetable(results, 'batch_results.csv');

fprintf('\nBatch completed: %d / %d successful scenarios.\n', nnz(solver_success), numCases);