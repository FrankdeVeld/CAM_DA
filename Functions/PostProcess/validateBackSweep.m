function output = validateBackSweep(input,outSim)
% VALIDATEBACKSWEEP  Validates the backward-sweep CAM optimisation output.
%
%   Reads the JSON produced by DAGreedyBackSweepJson, back-propagates the
%   primary and secondary from TCA to t_alert, applies the optimised
%   piecewise-constant RTN control forward in time, and plots:
%     (1) Danger metric history (miss distance^2, SMD, or PoC) vs thrust time
%     (2) Relative position in the B-plane at each node
%
%   USAGE
%     validateBackSweep(input,outSim) 
%

mu = 1.0;

%% Input
xp_tca    = input.xp_tCA(:);   % 6x1 
xs_tca    = input.xs_tCA(:);
P         = input.P;            % 3x3 covariance 
uMax      = input.uMax;         % Thrust magnitude
HBR       = input.HBR;         % HBR
N         = input.N;

%% ── 2. Extract control and timing from sim ───────────────────────────────
t_nodes   = outSim.t_nodes';
u_control = outSim.control;

%% ── 4. Forward integration with piecewise-constant RTN control ───────────
%  Re-integrate primary step by step applying the optimised u at each node.
%  Secondary propagates freely.

xp_hist     = nan(6, N);   % primary state history
xp_hist_nom = nan(6, N);   
rB_hist     = nan(N,2);   
smd_hist    = nan(N,1);
poc_hist    = nan(N,1);
deltaTca    = nan(N,1);

dt      = -diff(t_nodes);   % step size [s]
tcaOff  = outSim.deltaTca;

% Nominal primary trajectory
xp_hist_nom(:,end) = xp_tca;
for k = N-1:-1:1
    xp_hist_nom(:,k) = propKepOde(xp_hist_nom(:,k+1), zeros(3,1), -dt(k), mu);
end

% Quantities at TCA without control
deltaTca(end) = 0;
[Pb, rB_hist(end,:), smd_hist(end)] = Bplane(xp_tca, xs_tca, P);   
poc_hist(end)  = poc_Chan(HBR,Pb,smd_hist(end));

% Backward sweep
for j = N-1:-1:1
    if norm(u_control(j,:)) == 0
        deltaTca(j) = deltaTca(j+1);
        rB_hist(j,:) = rB_hist(j+1,:);
        smd_hist(j)  = smd_hist(j+1);
        poc_hist(j)  = poc_hist(j+1);
        continue
    end
    xp_hist(:,j) = xp_hist_nom(:,j);
    ddt = 0;
    xs_tcaOff = propKepOde(xs_tca, zeros(3,1), tcaOff(j), mu);
    for k = j:N-1
        if k == N-1; ddt = tcaOff(j); end
        u_k     = u_control(k,:)'*uMax;
        xp_hist(:,k+1) = propKepOde(xp_hist(:,k), u_k, dt(k) + ddt, mu);
    end

    [Pb, rB_hist(j,:), smd_hist(j)] = Bplane(xp_hist(:,end), xs_tcaOff, P);   
    poc_hist(j) = poc_Chan(HBR,Pb,smd_hist(j));
end
e2b = eci2Bplane(xp_hist(4:6,N),xs_tcaOff(4:6));
e2b = e2b([1 3],:); 
miss_dist  = normOfVec(rB_hist')';   

output = struct( ...
    'e2b',      e2b, ...
    't',        t_nodes, ...
    'rB',       rB_hist, ...
    'm_d',      miss_dist, ...
    'smd',      smd_hist, ...
    'poc',      poc_hist, ...
    'deltaTca', deltaTca ...
    );
end
