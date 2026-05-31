function [] = histograms(var,edges,varargin)
% Build histogram variable
num = size(var,1);
for k = 1:num
    counts(k,:) = histcounts(var(k,:),edges);
end
if nargin > 2
    normalizationType = varargin{1};
else
    normalizationType = 'none';
end
binWidth = edges(2)-edges(1);
switch lower(normalizationType)
    case 'none'
        counts = counts';
    case 'prob_density'
        counts = counts'/(size(var,2)*binWidth);
    case 'prob_mass'
        counts = counts'/size(var,2);
    case 'min-max'
        counts = (counts'-min(counts))/(max(counts)-min(counts));
    otherwise
        error('invalid normalization method')
end

% Plot histograms
bar(edges(2:end),counts,'group','FaceColor',0.6*ones(3,1),'EdgeColor','k');
ylabel('Probability')
grid on
box on
end