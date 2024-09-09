FileName = 'N100';

xp     = load(['.\write_read\xp_',FileName,'.dat']);

xs     = load(['.\write_read\xs_',FileName,'.dat']);
xpadj  = load(['.\write_read\xpadj_',FileName,'.dat']);

u      = load(['.\write_read\u_',FileName,'.dat']);
DM     = load(['.\write_read\DM_',FileName,'.dat']);

%% Orbits
figure()
hold on
scatter3(xp(:,1),xp(:,2),xp(:,3))
scatter3(xpadj(:,1),xpadj(:,2),xpadj(:,3))
scatter3(xs(:,1),xs(:,2),xs(:,3))
legend('primary','primary with control','secondary')
axis equal

%% Just the norms of the orbits
normxp = sqrt(xp(2:end,1).^2+ xp(2:end,2).^2+xp(2:end,3).^2);
normxpadj = sqrt(xpadj(1:(end-1),1).^2+ xpadj(1:(end-1),2).^2+xpadj(1:(end-1),3).^2);
t = [1:1:length(normxp)]';
figure()
hold on
plot(t,normxp, 'LineWidth', 2)
plot(t,normxpadj, 'LineWidth', 2)
xlabel('Node number (-)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Orbital distance r (km)', 'FontSize', 12, 'FontWeight', 'bold');

legend({'Nominal', 'Adjusted'}, 'Location', 'northeast', 'FontSize', 10, 'Box', 'off');

title('Orbital distance nominal and controlled orbit')
% Save the figure as a PNG file
saveas(gcf, ['Diff', FileName, '.png']);

% Save the figure as a FIG file
saveas(gcf, ['Diff', FileName, '.fig']);
%% Norm difference over time
figure()
plot(t,abs(normxpadj - normxp), 'LineWidth', 2)
t = [1:1:length(normxp)]';
grid on
xlabel('Node number (-)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Distance to nominal (km)', 'FontSize', 12, 'FontWeight', 'bold');

title('Distance controlled orbit to nominal, per node')
% Save the figure as a PNG file
saveas(gcf, ['DiffPerNode', FileName, '.png']);

% Save the figure as a FIG file
saveas(gcf, ['DiffPerNode', FileName, '.fig']);
%% 
% Assuming 'u' is a 3xN matrix where each row corresponds to R, T, and N directions
% 't' is the time variable
figure;
hold on;

t = [1:1:length(DM)]';

% Plotting each direction of control acceleration with specified line thickness
plot(t, u(:,1), 'r-', 'LineWidth', 2); % R direction in red (solid line)
plot(t, u(:,2), 'g--', 'LineWidth', 2); % T direction in green (dashed line)
plot(t, u(:,3), 'b-.', 'LineWidth', 2); % N direction in blue (dash-dot line)

% Adding grid for clarity
grid on;

% Adding labels to axes
xlabel('Node number (-)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Normalised Control Acceleration', 'FontSize', 12, 'FontWeight', 'bold');

% Adding title to the plot
title('Control Accelerations in R, T, N Directions', 'FontSize', 14, 'FontWeight', 'bold');

% Adding a legend inside the plot without a box
legend({'R Direction', 'T Direction', 'N Direction'}, 'Location', 'northwest', 'FontSize', 10, 'Box', 'off');

% Enhancing the axes for better visibility
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
% Save the figure as a PNG file
saveas(gcf, ['u', FileName, '.png']);

% Save the figure as a FIG file
saveas(gcf, ['u', FileName, '.fig']);
hold off;
