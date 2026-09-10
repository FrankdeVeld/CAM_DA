clear
% close all
% clc
addpath(genpath('.\Functions'))
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultAxesFontName','Times');
set(0,'DefaultUicontrolFontName','Times', 'DefaultUicontrolFontSize', 16);
set(0,'DefaultUitableFontName','Times', 'DefaultUitableFontSize', 16);
set(0,'DefaultTextFontName','Times', 'DefaultTextFontSize', 16);
set(0,'DefaultUipanelFontName','Times', 'DefaultUipanelFontSize', 16);

set(0, 'DefaultLineLineWidth', 1);
set(0,'defaultfigurecolor',[1 1 1])
    
%%
load('data/batch_results_md.mat');
t_thrust = -load('data/Julia_OC_Results_t0.csv');
diff = abs(start_time-t_thrust);
diff_prc = diff./start_time;
prc95y1 = prctile(diff_prc,95);
prc95y2 = prctile(diff,95);
prc95x = prctile(start_time,95);
% outliers = more_than(diff,95);
% outliers_prc = more_than(diff_prc,95);

%find outliers not respecting the correlation
ind1 = find(start_time < prc95x & diff_prc > prc95y1);
ind2 = find(start_time < prc95x & diff     > prc95y2);

f1 = figure;
[a,b] = sort(diff_prc);
subplot(2,1,1)
plot(start_time(b),a,'.')
hold on
plot(start_time(ind1),diff_prc(ind1),'.')
plot(prc95x*[1 1],[0,0.3],'k--')
plot([0 4000],prc95y1*[1,1],'k--')
hold off
ylabel('Relative error [%]')
grid on
box on


subplot(2,1,2)
[a,b] = sort(diff);
plot(start_time(b),a,'.')
hold on
plot(start_time(ind2),diff(ind2),'.')
plot(prc95x*[1 1],[0,1180],'k--')
plot([0 4000],prc95y2*[1,1],'k--')
hold off
xlabel('thrust time [s]')
ylabel('Absolute error [s]')
grid on
box on

saveFigurePDF(f1, 'outliers', 15, 10, 'centimeters')
