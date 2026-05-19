function h = plotMeanSigma(t, mu, sigma)
% plotMean3Sigma  Plot mean with shaded +/- sigma band
%
%   h = plotMean3Sigma(t, mu, sigma)
%
%   Inputs:
%       t     - t-axis values
%       mu    - mean values
%       sigma - standard deviation values
%
%   Outputs:
%       h     - struct with handles to patch and mean line
%
%   Example:
%       t = 1:length(mu);
%       plotMean3Sigma(t, mu, sigma);

    t = t(:);
    mu = mu(:);
    sigma = sigma(:);

    if ~isequal(length(t), length(mu), length(sigma))
        error('t, mu, and sigma must have the same length.');
    end

    upper = mu + sigma;
    lower = mu - sigma;
    
    % Shaded region
    h.patch = fill([t; flipud(t)], [upper; flipud(lower)], ...
                   [0.7 0.7 0.7], ...
                   'EdgeColor', 'none', ...
                   'FaceAlpha', 0.35);
    hold on;

    % Mean line
    h.mean = plot(t, mu, 'k', 'LineWidth', 2);

    % Optional styling
    grid on;
    hold off;
end