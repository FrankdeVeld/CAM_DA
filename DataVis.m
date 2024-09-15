clc
FileName = 'N100_Julia_Val';

xp     = load(['.\write_read\xp_',FileName,'.dat']);

xs     = load(['.\write_read\xs_',FileName,'.dat']);
xpadj  = load(['.\write_read\xpadj_',FileName,'.dat']);
tca  = load(['.\write_read\tca_',FileName,'.dat']);

DeltaRB_ECI  = load(['.\write_read\DeltaRB_',FileName,'.dat']);

u      = load(['.\write_read\u_',FileName,'.dat']);
DM     = load(['.\write_read\DM_',FileName,'.dat']);

%% Orbits
% figure()
% hold on
% scatter3(xp(:,1),xp(:,2),xp(:,3))
% scatter3(xpadj(:,1),xpadj(:,2),xpadj(:,3))
% scatter3(xs(:,1),xs(:,2),xs(:,3))
% legend('primary','primary with control','secondary')
% axis equal

%% Just the norms of the orbits
%normxp = sqrt(xp(2:end,1).^2+ xp(2:end,2).^2+xp(2:end,3).^2);
DistanceNodes = sqrt((xpadj(1:(end),1)-xp(2:end,1)).^2+ (xpadj(1:(end),2)-xp(2:end,2)).^2+(xpadj(1:(end),3)-xp(2:end,3)).^2);

%% Norm difference over time
figure()
t = [1:1:length(DM)]';
plot(t(1:end),DistanceNodes, 'LineWidth', 2)
grid on
xlabel('Node number (-)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Change in nominal node (km)', 'FontSize', 12, 'FontWeight', 'bold');

title('Distance adjusted node to nominal node')
% Save the figure as a PNG file
saveas(gcf, ['DiffPerNode', FileName, '.png']);

% Save the figure as a FIG file
saveas(gcf, ['DiffPerNode', FileName, '.fig']);
%% Plot control
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

%% DM over time 
t = [1:1:length(DM)]';
for i=1:length(DeltaRB_ECI(:,1))
    DM_From_DeltaRB(i) = norm(DeltaRB_ECI(i,:));
end
figure()
grid on
%hold on
semilogy(t,abs(DM_From_DeltaRB(:)-DM.^(1/2)), 'LineWidth', 2)
%plot(t,DM.^(1/2), 'LineWidth', 2)
%plot(t,DM_From_DeltaRB, 'LineWidth', 2)
xlabel('Node number (-)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Distance Metric difference', 'FontSize', 12, 'FontWeight', 'bold');
title('DM from direct DA and from DeltaRB', 'FontSize', 14, 'FontWeight', 'bold');
legend({'DM is a DA', 'DM from DeltaRB'}, 'Location', 'northeast', 'FontSize', 10, 'Box', 'off');

%% B-plane
figure()
scatter(DeltaRB_ECI(:,1),DeltaRB_ECI(:,3))

for i=1:length(DeltaRB_ECI(:,1))
    RelPos = xp(end,1:3)'-xs(end,1:3)';
    RelVel = xp(end,4:6)'-xs(end,4:6)';
    DeltaRB_BPlane(i,:) = ECI2B(DeltaRB_ECI(i,:)',RelPos,RelVel,xs(end,1:3)');
end

tf = 50:50:5000;
endPoint = 50;

SafeDis = 1;
Angle = 0:0.01:2*pi;

xSafe = cos(Angle).*SafeDis;
ySafe = sin(Angle).*SafeDis;

% Create figure
figure();
hold on;

% Plot safety circle
scatter(xSafe, ySafe, 2, 'k', 'DisplayName', 'Safety');

% Plot relative positions with color gradient based on time
scatter(DeltaRB_BPlane((end-endPoint+1):end,1), DeltaRB_BPlane((end-endPoint+1):end,3), 50, tf(1:endPoint)', 'filled', 'DisplayName', 'Primary');

% Plot secondary object at (0,0)
scatter(0, 0, 100, 'r', 'filled', 'DisplayName', 'Secondary');

% Add labels and title
xlabel('B1 axis [km]');
ylabel('B2 axis [km]');
title('Time till safety, DA first order control, 1e-4 m/s² thrust');

% Add color bar for time
c = colorbar;
c.Label.String = 'Time till safety [s]';

% Add legend
legend();

% Set axis equal for better representation of circles
axis equal;
hold off;

