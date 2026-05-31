function [] = plot_cdf(xdata,varargin)
%PLOT_CDF Plot the empirical cumulative distribution of a vector.
%   [f, x] = PLOT_CDF(xdata) returns the CDF values f and support x.
    if nargin > 1
        col = varargin{1};
    else
        col = 'k';
    end
    if nargin > 2
        line = varargin{2};
    else
        line = '-';
    end
    xdata = xdata(:);
    xdata = xdata(~isnan(xdata));

    [f, x] = ecdf(xdata);

    stairs(x, f, 'LineWidth', 1.5,'Color',col,'LineStyle',line)
    grid on
    ylabel('Cumulative probability')
end