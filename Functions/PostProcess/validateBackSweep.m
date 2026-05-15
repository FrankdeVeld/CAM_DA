function output = validateBackSweep(jsonFile,scenario)
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
tCA_Nom   = jin.t_back;      % nominal TCA time from epoch [s]
uMax      = jin.uMax;         % physical thrust magnitude [km/s^2]
HBR       = jin.HBR;         % HBR

% Convert covariance to m^2 if needed – here we keep km^2 to match Bplane.m
P = P_in;

%% ── 2. Extract control and timing from JSON ───────────────────────────────
t_nodes    = zeros(N,1);
u_control  = zeros(N,3);   % RTN, already scaled (u_dir * uMax from JSON)
rB_opt     = zeros(N,2);

for k = 1:N
    nd = nodes(k);
    t_nodes(k)   = nd.tNode;         % time of node from epoch [s]
    u_control(k,:) = nd.controlRTN(:)'; % RTN thrust [km/s^2]
    rB_opt(k,:)  = nd.relativePositionBPlane(:)';
end

% Thrust duration = tCA_Nom - first active node time
t_alert = tCA_Nom;
dt      = tCA_Nom / (N-1);   % step size [s]

%% ── 3. Back-propagate both objects from TCA to t_alert ───────────────────

xp_backprop = propKepOde(xp_tca, zeros(3,1), -t_alert, mu);  % 6 x (N+1)

% Initial states at t_alert (last column of back-propagation)
xp0 = xp_backprop(:, end);   % 6x1 at t_alert

%% ── 4. Forward integration with piecewise-constant RTN control ───────────
%  Re-integrate primary step by step applying the optimised u at each node.
%  Secondary propagates freely.

N_ref = 50;
dt1   = 0.02/scenario.Tsc;
ddt   = dt1/N_ref;

xp_hist     = nan(6, N);   % primary state history
xp_ref      = nan(6, N_ref); 
xs_ref      = nan(6, N_ref); 
xp_hist_nom = nan(6, N);   
rB_hist     = nan(N,2);   
smd_hist    = nan(N,1);
poc_hist    = nan(N,1);
deltaTca    = nan(N,1);

xp_hist_nom(:,1) = xp0;

% % secondary refined near tca
% xs_ref(:,1) = propKepOde(xs_tca, zeros(3,1), -dt1, mu);
% for i = 1:N_ref-1
%     xs_ref(:,i+1) = propKepOde(xs_ref(:,i), zeros(3,1), ddt, mu);
% end

% Nominal primary trajectory
for k = 1:N-1
    xp_hist_nom(:,k+1) = propKepOde(xp_hist_nom(:,k), zeros(3,1), dt, mu);
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
    for k = j:N-1
        u_k     = u_control(k,:)'*uMax;
        xp_hist(:,k+1) = propKepOde(xp_hist(:,k), u_k, dt, mu);
    end

    % find new TCA
    % xp_ref(:,1) = propKepOde(xp_hist(:,N-1), u_k, dt-dt1, mu);
    % for i = 1:N_ref-1
    %     x_rel_ref(:,i) = xp_ref(:,i) - xs_ref(:,i);
    %     ca_ref(i)      = norm(x_rel_ref(1:3,i)');
    %     % if i > 1 && ca_ref(i) > ca_ref(i-1); break; end
    %     xp_ref(:,i+1)  = propKepOde(xp_ref(:,i), u_k, ddt, mu);
    % end
    % [~,ind] = min(ca_ref);

    % deltaTca(j) =  ddt * (ind - (N_ref + 1)); 
    [Pb, rB_hist(j,:), smd_hist(j)] = Bplane(xp_hist(:,end), xs_tca, P);   
    % [Pb, rB_hist(j,:), smd_hist(j)] = Bplane(xp_ref(:,ind), xs_ref(:,ind), P);   
    poc_hist(j) = poc_Chan(HBR,Pb,smd_hist(j));
end
miss_dist  = normOfVec(rB_hist')';   

output = struct( ...
    'rB',       rB_hist, ...
    'm_d',      miss_dist, ...
    'smd',      smd_hist, ...
    'poc',      poc_hist, ...
    'deltaTca', deltaTca ...
    );
end
