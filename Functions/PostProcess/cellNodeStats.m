function [mu, sigma] = cellNodeStats(C)
% cellNodeStats  Nodewise mean and std across a cell array of numeric arrays
%
%   [mu, sigma] = cellNodeStats(C)
%
%   C must be a cell array where each cell contains a numeric array of the
%   same size. The function computes the mean and standard deviation across
%   the cell entries, node by node.

    if ~iscell(C) || isempty(C)
        error('Input must be a non-empty cell array.');
    end
    
    % Check that all cells have the same size
    refSize = size(C{1});
    for k = 2:numel(C)
        if ~isequal(size(C{k}), refSize)
            error('All cells must contain arrays of the same size.');
        end
    end
    
    % for k = 2:numel(C)
    %     p = C{k};
    %     p(isnan(p)) = 
    % end

    % Stack along 3rd dimension
    X = cat(ndims(C{1}) + 1, C{:});

    % Nodewise statistics across cells
    mu    = mean(X, ndims(C{1}) + 1);
    sigma = std(X, 0, ndims(C{1}) + 1);
end