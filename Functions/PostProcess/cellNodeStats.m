function [mu, sigma] = cellNodeStats(C,yesnan)
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
    
    % Stack along 3rd dimension
    X = abs(cat(ndims(C{1}) + 1, C{:}));
    if yesnan
        for j = size(X,1)-1:-1:2
            if sum(squeeze(X(j,:,:))) == 0
                X = X(j+1:end,:,:);
                break
            end
        end
        X(X==0) = nan;
    end
    % Nodewise statistics across cells
    mu    = mean(X, ndims(C{1}) + 1,'omitnan');
    sigma = std(X, 0, ndims(C{1}) + 1,'omitnan');
end