close all
clc
load('data/batch_results_360.mat');
prc = 95;
rel_err = less_than(rel_err,prc);
dvTot = less_than(dvTot,prc);
max_tca_shift = less_than(max_tca_shift,prc);
comp_time = less_than(comp_time,prc);
start_time = less_than(start_time,prc);
xmin = min(dvTot);

t_vec = flip(linspace(0,params.n_orb,params.nx_orb*params.n_orb+1));

col = [0, 0.4470, 0.7410];

f3 = figure;
plot_cdf(comp_time,col)
xlabel('Computation time [ms]')

f4 = figure;
plot_cdf(start_time,col)
xlabel('Thrust duration [s]')

%%
load('data/batch_results.mat');
rel_err = less_than(rel_err,prc);
dvTot = less_than(dvTot,prc);
max_tca_shift = less_than(max_tca_shift,prc);
comp_time = less_than(comp_time,prc);
start_time = less_than(start_time,prc);
xmax = max(dvTot);

t_vec = flip(linspace(0,params.n_orb,params.nx_orb*params.n_orb+1));
col = [0.8500, 0.3250, 0.0980];

f3 = figure(f3);
hold on
plot_cdf(comp_time,col)
hold off

f4 = figure(f4);
hold on
plot_cdf(start_time,col)

%%
load('data/batch_results_60.mat');
rel_err = less_than(rel_err,prc);
dvTot = less_than(dvTot,prc);
max_tca_shift = less_than(max_tca_shift,prc);
comp_time = less_than(comp_time,prc);
start_time = less_than(start_time,prc);
xmax = max(dvTot);

t_vec = flip(linspace(0,params.n_orb,params.nx_orb*params.n_orb+1));
col = [0.9290 0.6940 0.1250];


f3 = figure(f3);
hold on
plot_cdf(comp_time,col)
hold off
legend('N/orbit = 360','N/orbit = 120','N/orbit = 60','Location','southeast')

f4 = figure(f4);
hold on
plot_cdf(start_time,col)
legend('N/orbit = 360','N/orbit = 120','N/orbit = 60','Location','southeast')
