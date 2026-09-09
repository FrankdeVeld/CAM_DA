function [] = postProcessBatch()
%POSTPROCESSBATCH Summary of this function goes here
%   Detailed explanation goes here

load('data/batch_results.mat');
prc = 95;
rel_err = less_than(rel_err,prc);
dvTot = less_than(dvTot,prc);
max_tca_shift = less_than(max_tca_shift,prc);
comp_time = less_than(comp_time,prc);
start_time = less_than(start_time,prc);
xmin = min(dvTot);

t_vec = flip(linspace(0,params.n_orb,params.nx_orb*params.n_orb+1));

col = [0, 0.4470, 0.7410];

f2 = figure;
plot_cdf(rel_err/100,col)
xlabel('SMD relative error')

f4 = figure;
plot_cdf(max_tca_shift,col)
xlabel('TCA shift [s]')

f5 = figure;
plot_cdf(comp_time,col)
xlabel('Computation time [ms]')

f6 = figure;
plot_cdf(start_time,col)
xlabel('Thrust duration [s]')

f7 = figure;
[mu_R,sigma_R] = cellNodeStats(uR_all,'true');
[mu_T,sigma_T] = cellNodeStats(uT_all,'true');
[mu_N,sigma_N] = cellNodeStats(uN_all,'true');

subplot(3,1,1)
plotMeanSigma(t_vec,mu_R,sigma_R);
axis tight
ylabel('$|u_R|$')
subplot(3,1,2)
plotMeanSigma(t_vec,mu_T,sigma_T);
axis tight
ylabel('$|u_T|$')
subplot(3,1,3)
plotMeanSigma(t_vec,mu_N,sigma_N);
axis tight
xlabel('Orbits before TCA')
ylabel('$|u_N|$')
ylim([0,1])

f8 = figure;
[mu_err,sigma_err] = cellNodeStats(smd_rel_err_all,'false');

plotMeanSigma(t_vec,mu_err,sigma_err);
xlabel('Orbits before TCA')
ylabel('SMD relative error')

%% Miss distance

load('data/batch_results_md.mat');
prc = 95;
rel_err = less_than(rel_err,prc);
dvTot = less_than(dvTot,prc);
max_tca_shift = less_than(max_tca_shift,prc);
comp_time = less_than(comp_time,prc);
start_time = less_than(start_time,prc);
xmax = max(dvTot);

t_vec = flip(linspace(0,params.n_orb,params.nx_orb*params.n_orb+1));
col = [0.8500, 0.3250, 0.0980];

f2 = figure(f2);
hold on
plot_cdf(rel_err/100,col)
hold off

f4 = figure(f4);
hold on
plot_cdf(max_tca_shift,col)
hold off
legend('SMD','MD','Location','southeast')

f5 = figure(f5);
hold on
plot_cdf(comp_time,col)
hold off
legend('SMD','MD','Location','southeast')

f6 = figure(f6);
hold on
plot_cdf(start_time,col)

f9 = figure;
[mu_R,sigma_R] = cellNodeStats(uR_all,'true');
[mu_T,sigma_T] = cellNodeStats(uT_all,'true');
[mu_N,sigma_N] = cellNodeStats(uN_all,'true');

subplot(3,1,1)
plotMeanSigma(t_vec(length(t_vec)-length(mu_R)+1:end),mu_R,sigma_R);
ylabel('$|u_R|$')
axis tight
xlim([0,0.65])
subplot(3,1,2)
plotMeanSigma(t_vec(length(t_vec)-length(mu_T)+1:end),mu_T,sigma_T);
ylabel('$|u_T|$')
axis tight
xlim([0,0.65])
subplot(3,1,3)
plotMeanSigma(t_vec(length(t_vec)-length(mu_N)+1:end),mu_N,sigma_N);
xlabel('Orbits before TCA')
ylabel('$|u_N|$')
axis tight
xlim([0,0.65])

f10 = figure;
[mu_err,sigma_err] = cellNodeStats(md_rel_err_all,'false');

plotMeanSigma(t_vec,mu_err,sigma_err);
axis tight
xlabel('Orbits before TCA')
ylabel('MD relative error')
xlim([0,0.65])
% 
% %% FO SMD
% load('data/FO_comparison_SMD.mat')
% col = [0, 0.4470, 0.7410];
% 
% f6 = figure(f6);
% plot_cdf(t_thrust,col,'--')
% 
% 
% %% FO MD
% load('data/FO_comparison_MD.mat')
% col = [0.8500, 0.3250, 0.0980];
% 
% f6 = figure(f6);
% plot_cdf(t_thrust,col,'--')

%% OCP time optimal
t_thrust = -load('data/Julia_OC_Results_t0.csv');
col = [0.8500, 0.3250, 0.0980];
f6 = figure(f6);
plot_cdf(t_thrust,col,'-.')
hold off
legend('SMD-greedy','MD-greedy','SMD-FO','MD-FO','MD-OC')
ax1 = gca;
box off
ax2 = axes('Position', ax1.Position, ...
           'Color', 'none', ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'right', ...
           'XLim', [xmin, xmax], ...
           'YTick', [], ...
           'YColor', 'black'...
           );
xlabel(ax2,'$\Delta v$ [m/s]')


t_0_diff = (start_time - t_thrust)./t_thrust*100;
t_0_diff = less_than(t_0_diff,prc);
f11 = figure();
histograms(t_0_diff',linspace(0,max(t_0_diff),30),'prob_density')
hold on
plot_cdf(t_0_diff)
m = mean(t_0_diff);
m1 = median(t_0_diff);
plot([m, m],[0,1],'k--')
plot([m1, m1],[0,1],'k-.')
xlabel('Relative error [%]')

%% Save all figures
saveFigurePDF(f2, 'rel_err', 15, 10, 'centimeters')
saveFigurePDF(f4, 'Figures/tca_shift', 15, 6, 'centimeters')
saveFigurePDF(f5, 'Figures/comp_time', 15, 6, 'centimeters')
saveFigurePDF(f6, 'Figures/start_time', 20, 10, 'centimeters')
saveFigurePDF(f7, 'Figures/control_smd', 15, 10, 'centimeters')
saveFigurePDF(f8, 'Figures/error_smd', 15, 6, 'centimeters')
saveFigurePDF(f9, 'Figures/control_md', 15, 10, 'centimeters')
saveFigurePDF(f10, 'Figures/error_md', 15, 6, 'centimeters')
saveFigurePDF(f11, 'Figures/error_t0', 15, 6, 'centimeters')
close all
end