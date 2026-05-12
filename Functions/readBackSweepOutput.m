function [control, rB, tca, m_d] = readBackSweepOutput(jsonFile)
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

control  = zeros(N, 3);
rB       = zeros(N, 3);
tca_shift    = zeros(N, 1);
tca_absolute = zeros(N, 1);
m_d             = zeros(N, 1);
for k = 1:N
    nd = nodes(k);
    control(k, :)  = nd.controlRTN(:)';
    rB(k, :)       = nd.relativePositionBPlane_km(:)';
    tca_shift(k)   = nd.tcaShift_s;
    tca_absolute(k)= nd.tcaAbsolute_s;
    m_d(k)         = nd.dangerMetric_km2;
end

tca.shift_s    = tca_shift;
tca.absolute_s = tca_absolute;

end