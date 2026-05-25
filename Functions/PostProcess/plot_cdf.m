function [] = plot_cdf(xdata,col)
%PLOT_CDF Plot the empirical cumulative distribution of a vector.
%   [f, x] = PLOT_CDF(xdata) returns the CDF values f and support x.

    xdata = xdata(:);
    xdata = xdata(~isnan(xdata));

    [f, x] = ecdf(xdata);

    stairs(x, f, 'LineWidth', 1.5,'Color',col)
    grid on
    ylabel('Cumulative probability')
end