function [] = postProcessBatch()
%POSTPROCESSBATCH Summary of this function goes here
%   Detailed explanation goes here

load('batch_results.mat');
prc = 95;
abs_err = less_than(abs_err,prc);
rel_err = less_than(rel_err,prc);
dvTot = less_than(dvTot,prc);
max_tca_shift = less_than(max_tca_shift,prc);
comp_time = less_than(comp_time,prc);
start_time = less_than(start_time,prc);

t_vec = flip(linspace(0,params.n_orb,params.nx_orb*params.n_orb+1));

% figure
% violin(abs_err);
% figure
% violin(rel_err);
% figure
% violin(dvTot);
% figure
% violin(max_tca_shift);

figure
plot_cdf(abs_err)
xlabel('Absolute error [m]')

figure
plot_cdf(rel_err/100)
xlabel('SMD relative error')


figure
plot_cdf(dvTot)
xlabel('Total $\Delta v$ [m/s]')


figure
plot_cdf(max_tca_shift)
xlabel('TCA shift [s]')


figure
plot_cdf(comp_time)
xlabel('Computation time [ms]')


figure
plot_cdf(start_time)
xlabel('Start time [s]')

figure
[mu_R,sigma_R] = cellNodeStats(uR_all);
[mu_T,sigma_T] = cellNodeStats(uT_all);
[mu_N,sigma_N] = cellNodeStats(uN_all);

subplot(3,1,1)
plotMeanSigma(t_vec,mu_R,sigma_R);
subplot(3,1,2)
plotMeanSigma(t_vec,mu_T,sigma_T);
subplot(3,1,3)
plotMeanSigma(t_vec,mu_N,sigma_N);

figure
[mu_err,sigma_err] = cellNodeStats(smd_abs_err_all);

plotMeanSigma(t_vec,mu_err,sigma_err);

end

function out = less_than(in,prc)
    in  = in(:);
    in  = in(~isnan(in));
    out = in(in<prctile(in,prc));
end