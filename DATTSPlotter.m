
%% 
clc;
FileName = 'N100_Julia_Val_TTS';

% Load data
xpf     = load(['.\write_read\xpf_',FileName,'.dat']);
xsf     = load(['.\write_read\xsf_',FileName,'.dat']);
tf      = load(['.\write_read\tf_',FileName,'.dat']);

% Invert tf
tf = 4000 - tf;

% Define safe distance and angle for safety circle
SafeDis = 1;
Angle = 0:0.01:2*pi;

xSafe = cos(Angle).*SafeDis;
ySafe = sin(Angle).*SafeDis;

% Calculate relative positions
relPosX = xpf(:,1) - xsf(:,1);
relPosY = xpf(:,3) - xsf(:,3);

% Create figure
figure();
hold on;

% Plot safety circle
scatter(xSafe, ySafe, 2, 'k', 'DisplayName', 'Safety');

% Plot relative positions with color gradient based on time
scatter(relPosX, relPosY, 50, tf, 'filled', 'DisplayName', 'Primary');

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

