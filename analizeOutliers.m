clear
close all
clc

load('data/batch_results_md.mat');
t_thrust = -load('data/Julia_OC_Results_t0.csv');
diff = abs(start_time-t_thrust);
diff_prc = diff./start_time;
% outliers = more_than(diff,95);
% outliers_prc = more_than(diff_prc,95);

f1 = figure;
[a,b] = sort(diff_prc);
prc95y = prctile(diff_prc,95);
prc95x = prctile(start_time,95);
subplot(2,1,1)
plot(start_time(b),a,'.')
hold on
plot(prc95x*[1 1],[0,0.3],'k--')
plot([0 4000],prc95y*[1,1],'k--')
hold off
ylabel('Relative error [%]')
grid on
box on


subplot(2,1,2)
[a,b] = sort(diff);
prc95y = prctile(diff,95);
plot(start_time(b),a,'.')
hold on
plot(prc95x*[1 1],[0,1180],'k--')
plot([0 4000],prc95y*[1,1],'k--')
hold off
xlabel('thrust time [s]')
ylabel('Absolute error [s]')
grid on
box on

saveFigurePDF(f1, 'outliers', 15, 10, 'centimeters')
