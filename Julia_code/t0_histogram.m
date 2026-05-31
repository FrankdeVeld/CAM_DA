% 1. Load results 
sub_optimal = readmatrix('C:\Users\frank\Documents\GitHub\CAM_DA\Julia_code\Matlab_Greedy_Results_t0.csv'); 
optimal = readmatrix('C:\Users\frank\Documents\GitHub\CAM_DA\Julia_code\Julia_OC_Results_t0.csv'); 

% 2. Calculate the relative discrepancy
discrepancyrel = (optimal(:) - sub_optimal(:).*-1)./(sub_optimal(:)).*100;

% 3. Create figure and make it much wider than it is tall
%figure;
%set(gcf, 'Position', [100, 100, 800, 300]); 

% 4. Plot Histogram and make it grey
col = [0, 0.4470, 0.7410];

axes( 'FontName', 'Times New Roman', 'FontSize', 16,'NextPlot', 'add');
f1 = figure(1);
hold on
%grid on;
box on;
plot_cdf(discrepancyrel, col);
hold off
%h = histogram(discrepancyrel);
%h.FaceColor = [0.7 0.7 0.7]; % Light grey
%h.EdgeColor = 'black';       % Crisp outlines for the bars

% 5. Format axes (Boxed, Grid, and Times New Roman font)

%set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);

% 6. Add Labels
xlabel('Relative error in t_0 [%]', 'FontName', 'Times New Roman','FontSize',16);
ylabel('Cumulative Probability', 'FontName', 'Times New Roman','FontSize',16);

%ylim([0, max(h.Values) * 1.15]); 

%med_val = median(discrepancyrel);
%xline(med_val, '--k', sprintf('Median: %.3f%%', med_val), ...
%    'FontName', 'Times New Roman', 'FontSize', 11, ...
%    'LabelOrientation', 'horizontal', ...    % Makes text horizontal
%    'LabelHorizontalAlignment', 'right', ... % Moves text to the right side of the line
%    'LabelVerticalAlignment', 'top', ...
%    'LineWidth', 1.5);

saveFigurePDF(f1, 't0_comp_cdf', 15, 6, 'centimeters')