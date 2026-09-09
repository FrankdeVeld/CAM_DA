load('data/batch_results_md.mat');
prc = 100;
rel_err = less_than(rel_err,prc);
dvTot = less_than(dvTot,prc);
max_tca_shift = less_than(max_tca_shift,prc);
comp_time = less_than(comp_time,prc);
start_time = less_than(start_time,prc);
xmin = min(dvTot);

t_thrust = -load('data/Julia_OC_Results_t0.csv');
