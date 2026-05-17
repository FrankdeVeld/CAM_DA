function output = readBackSweepOutput(jsonFile,input,scenario,metric)
% READBACKSWEEPOUTPUT  Read the output JSON from DAGreedyBackSweepJson and
%                      extract nodewise control, B-plane relative position,
%                      and TCA history as plain MATLAB arrays.
%
%   [control, rB, tca] = readBackSweepOutput(jsonFile)
%
%   INPUT
%     jsonFile  : path to the output JSON file (default: 'output.json')
%
%   OUTPUT
%     control   : N x 3 array, optimal RTN thrust direction at each node
%                 [-, normalised; multiply by uMax to get physical units]
%     rB        : N x 3 array, relative position in B-plane at each node [km]
%     tca       : struct with fields
%                   .shift_s    N x 1, TCA shift w.r.t. nominal TCA [s]
%                   .absolute_s N x 1, absolute TCA from epoch [s]
%
%   EXAMPLE
%     [u, rB, tca] = readBackSweepOutput('output.json');
%     figure; plot(tca.absolute_s, vecnorm(rB,2,2)); ylabel('miss distance [km]')

if nargin < 1
    jsonFile = 'output.json';
end

raw = jsondecode(fileread(jsonFile));
nodes = raw.nodes;
N = numel(nodes);

control         = nan(N, 3);
rB              = nan(N, 2);
tca_shift       = nan(N, 1);
m_d             = nan(N, 1);
smd             = nan(N, 1);
poc             = nan(N, 1);
for k = N:-1:1
    nd = nodes(k);
    t_nodes(k)     = nd.tNode;         % time of node back from TCA
    control(k, :)  = nd.controlRTN(:)';
    if norm(control(k, :)) == 0 && k < N
        rB(k, :)     = rB(k+1, :);
        tca_shift(k) = tca_shift(k+1);
        m_d(k)       = m_d(k+1);
        smd(k)       = smd(k+1);
    else
        t_start        = nd.tNode;
        rB(k, :)       = nd.relativePositionBPlane(:)';
        tca_shift(k)   = nd.tcaShift;
        if metric == 1
            m_d(k) = nd.dangerMetric;
            smd(k) = nan;
            poc(k) = nan;
        else
            m_d(k) = norm(rB(k, :))^2;
            smd(k) = nd.dangerMetric;
            poc(k) = poc_Chan(input.HBR,scenario.Pb,smd(k));
        end
    end
end
control(end, :)  = control(end-1, :);
output = struct( ...
    't_start',  t_start, ...
    't_nodes',  t_nodes, ...
    'control',  control, ...
    'rB',       rB, ...
    'm_d',      m_d, ...
    'smd',      smd, ...
    'deltaTca', tca_shift, ...
    'poc',      poc ...
    );

end