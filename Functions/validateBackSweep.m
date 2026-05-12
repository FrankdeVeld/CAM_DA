function [rB_hist, miss_dist, smd_hist, poc_hist] = validateBackSweep(jsonFile)
% VALIDATEBACKSWEEP  Validates the backward-sweep CAM optimisation output.
%
%   Reads the JSON produced by DAGreedyBackSweepJson, back-propagates the
%   primary and secondary from TCA to t_alert, applies the optimised
%   piecewise-constant RTN control forward in time, and plots:
%     (1) Danger metric history (miss distance^2, SMD, or PoC) vs thrust time
%     (2) Relative position in the B-plane at each node
%
%   USAGE
%     validateBackSweep()                 % uses 'output.json' in current folder
%     validateBackSweep('path/out.json')
%
%   DEPENDENCIES (your existing helpers)
%     propKepOde, eci2Bplane, Bplane
%
%   Author:  auto-generated validation script

if nargin < 1, jsonFile = 'output.json'; end

mu = 1.0;

%% ── 1. Read JSON ──────────────────────────────────────────────────────────
raw = jsondecode(fileread(jsonFile));
nodes = raw.nodes;
N     = raw.N;

% Input states: TCA states in km / km/s
jin       = jsondecode(fileread('./input.json'));
xp_tca    = jin.xp_tCA(:);   % 6x1 [km, km/s]
xs_tca    = jin.xs_tCA(:);
P_in      = jin.P;            % 3x3 covariance [km^2]
tCA_Nom_s = jin.t_back;      % nominal TCA time from epoch [s]
uMax      = jin.uMax;         % physical thrust magnitude [km/s^2]
HBR       = jin.HBR;         % HBR

% Convert covariance to m^2 if needed – here we keep km^2 to match Bplane.m
P = P_in;

%% ── 2. Extract control and timing from JSON ───────────────────────────────
t_nodes    = zeros(N,1);
u_control  = zeros(N,3);   % RTN, already scaled (u_dir * uMax from JSON)
tCA_abs    = zeros(N,1);
rB_opt     = zeros(N,3);

for k = 1:N
    nd = nodes(k);
    t_nodes(k)   = nd.tNode_s;         % time of node from epoch [s]
    u_control(k,:) = nd.controlRTN(:)'; % RTN thrust [km/s^2]
    tCA_abs(k)   = nd.tcaAbsolute_s;
    rB_opt(k,:)  = nd.relativePositionBPlane_km(:)';
end

% Thrust duration = tCA_Nom - first active node time
t_alert = tCA_Nom_s;
dt      = tCA_Nom_s / (N-1);   % step size [s]

%% ── 3. Back-propagate both objects from TCA to t_alert ───────────────────

xp_backprop = propKepOde(xp_tca, zeros(3,1), -t_alert, mu);  % 6 x (N+1)

% Initial states at t_alert (last column of back-propagation)
xp0 = xp_backprop(:, end);   % 6x1 at t_alert

%% ── 4. Forward integration with piecewise-constant RTN control ───────────
%  Re-integrate primary step by step applying the optimised u at each node.
%  Secondary propagates freely.

xp_hist     = nan(6, N);   % primary state history
xp_ref       = nan(6, 100); % primary state history
xp_hist_nom = nan(6, N);   
miss_dist   = zeros(1, N);
rB_hist     = zeros(N, 2);   
smd_hist    = zeros(1, N);
poc_hist    = zeros(1, N);
deltaTca    = zeros(1, N);

xp_hist_nom(:,1) = xp0;

% secondary refined near tca
% N_ref = 500;
% ddt = dt/N_ref;
% xs_ref(:,1) = propKepOde(xs_tca, zeros(3,1), -dt, mu);
% for i = 1:2*N_ref
%     xs_ref(:,i+1) = propKepOde(xs_ref(:,i), zeros(3,1), ddt, mu);
% end

for k = 1:N-1
    xp_hist_nom(:,k+1) = propKepOde(xp_hist_nom(:,k), zeros(3,1), dt, mu);
end

for j = 1:N-1
    xp_hist(:,j) = xp_hist_nom(:,j);
    for k = j:N-1
        u_k     = u_control(k,:)'*uMax;
        xp_hist(:,k+1) = propKepOde(xp_hist(:,k), u_k, dt, mu);
    end
   
    % %find new TCA
    % xp_ref(:,1) = xp_hist(:,end-1);
    % for i = 1:2*N_ref
    %     xp_ref(:,i+1) = propKepOde(xp_ref(:,i), u_k, ddt, mu);
    % end
    % md_ref  = normOfVec(xp_ref - xs_ref);
    % [~,ind] = min(md_ref);
    % 
    % deltaTca(j) = dt - ddt*ind; 
    [Pb, rB_hist(j,:), smd_hist(j)] = Bplane(xp_hist(:,end), xs_tca, P);   
    miss_dist(j) = norm(rB_hist(j,:));   
    poc_hist(j) = poc_Chan(HBR,Pb,smd_hist(j));

end

[Pb, rB_hist(end,:), smd_hist(end)] = Bplane(xp_tca, xs_tca, P);   
poc_hist(end)  = poc_Chan(HBR,Pb,smd_hist(end));
miss_dist(end) = norm(rB_hist(end,:));
deltaTca(end)  = 0;
end
