
% 1. DATA LOADING

dr_Julia  = readmatrix('Julia_6_rB.csv');    % Nx2 array
t_Julia = readmatrix('Julia_6_t.csv'); % Nx1 array
dr_Matlab  = readmatrix('Matlab_6_rB.csv');    % Mx2 array
t_Matlab = readmatrix('Matlab_6_t.csv'); % Mx1 array

% dr_Julia  = readmatrix('Julia_6_highecc_rB.csv');    % Nx2 array
% t_Julia = readmatrix('Julia_6_highecc_t.csv'); % Nx1 array
% dr_Matlab  = readmatrix('Matlab_6_highecc_rB.csv');    % Mx2 array
% t_Matlab = readmatrix('Matlab_6_highecc_t.csv'); % Mx1 array

dr_Julia = dr_Julia(1:15:end, :);
t_Julia = t_Julia(1:15:end, :);
% Julia files have too many entries
Fontsize_plots = 10;
f1 = figure('Position', [100, 100, 800, 500]);
hold on;

% Plot the circle with radius 2 centered at (0,0)
theta = linspace(0, 2*pi, 200);
R = 2;
plot(R*cos(theta), R*sin(theta), 'k-', 'LineWidth', 1);

% Plot the little star at the center (0,0)
plot(0, 0, 'k*', 'MarkerSize', 8, 'LineWidth', 1.2);

% Plot the second set of 2D vectors (Scatter plot with BLACK outline)
s1=scatter(dr_Matlab(:,1).*2, dr_Matlab(:,2).*2, 15, t_Matlab.*5928.971459998461, 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.2);

% Plot the first set of 2D vectors (Scatter plot with MAGENTA outline)
% scatter(X, Y, MarkerSize, ColorData)
s2=scatter(dr_Julia(:,1).*2, dr_Julia(:,2).*2, 15, t_Julia.*5928.971459998461, 'filled', ...
    'MarkerEdgeColor', 'm', 'LineWidth', 0.2);


% .*2 to un-normalise, same for the time

% Set the colormap (Winter matches the Blue-to-Green scale in your image)
colormap('winter');

% Add legend in the bottom right ('southeast')
legend([s1, s2], {'Greedy', 'OC Bench'}, 'Location', 'southeast', 'FontSize', Fontsize_plots);

% Add and format the colorbar
cb = colorbar;
clim([0 3500]) 
ylabel(cb, 'Time before TCA [s]', 'Rotation', 270, ...
    'VerticalAlignment', 'bottom', 'FontName', 'Times New Roman', 'FontSize', Fontsize_plots);

% --- Adjust colorbar thickness and position ---
cbPos = cb.Position;        % Get the current position [left, bottom, width, height]
cbPos(3) = cbPos(3) * 0.8;  % Reduce the width (thickness) to 40%
cbPos(1) = cbPos(1) - 0.11; % Shift it to the left (closer to the plot)
cb.Position = cbPos;        % Apply the new position

% Axis labels using LaTeX interpreter for Greek symbols
xlabel('$\xi$ [km]', 'Interpreter', 'latex');
ylabel('$\zeta$ [km]', 'Interpreter', 'latex');

% Enable Grid and Box
grid on;
box on;

% Set 1:1 aspect ratio so the circle doesn't look stretched
axis equal; 

xlim([-2.5, 2.5]);
ylim([-2.5, 2.5]);

% Apply Font Name and Size to the axes
set(gca, 'FontSize', Fontsize_plots, 'FontName', 'Times New Roman');

hold off;

saveFigurePDF(f1, 'bplane_comp_cs6_highecc_traj', 20, 8, 'centimeters')