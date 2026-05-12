% y0     = [1 0 0.05 0 0 pi/6]';
% alpha  = deg2rad(40);
% beta   = deg2rad(140);
% dv_dir = [sin(alpha)*sin(beta); cos(beta); cos(alpha)*sin(beta)];
% % dv_dir = normalize([1 0 0]','norm');
clear
[primary,secondary] = generateInitShort(1);
N = 100;

% Non-dimensionalization
Lsc   = primary.a;
musc  = 398600.4418;
Vsc   = sqrt(musc/Lsc);
Tsc   = Lsc/Vsc;
Asc   = Vsc/Tsc;
scale = diag([ones(3,1)./Lsc; ones(3,1)./Vsc]);

x_p     = scale*primary.x0;
x_s     = scale*secondary.x0;
ctrlMax = primary.ctrlMax/Asc;
T       = primary.T/Tsc;
HBR     = (primary.HBR + secondary.HBR)/Lsc;

epsilon = 1e-1;
sigma_eps = 0.25;
sigma   = sigma_eps * epsilon;
t_min = -2*pi;

cam = optimal_cam(x_p, x_s, epsilon, sigma, t_min, ...
    'K',        720,    ...
    'Ntime',    100,   ...
    'AutoStop', false,   ...
    'FillTol',  0.02,   ...
    'Verbose',  true);