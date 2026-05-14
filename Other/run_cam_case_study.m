% RUN_CAM_CASE_STUDY.M
% =========================================================================
% Demonstrates optimal_cam() for a synthetic conjunction geometry built
% from the same paper case study (Dell'Elce, de Veld, Pomet AAS 24-458).
%
% We place the secondary at a known position INSIDE the safe circle so that
% the conjunction is well-defined, then call optimal_cam() to find the
% optimal maneuver start time and thrust direction that exits the safe circle.
% =========================================================================
clc; clear; close all;

%% ---- Primary equinoctial elements (same as paper case study) ------------
% y0 = [1; 0; 0.05; 0; 0; pi/6];   % [a, p1, p2, q1, q2, L]
coe = [6780; 0.001; deg2rad(51.6); 0; 0; pi/6];
%% ---- Convert to Cartesian ECI -------------------------------------------
% (equinoctial_to_cartesian is a local function of minimum_warning_time;
%  we replicate the call here for clarity.  You can also inline it.)
% [rp, vp] = eq2cart(y0);   % see local function below
% xp = [rp; vp];
xp = COE2RV(coe);
Lsc   = coe(1);
musc  = 398600.4418;
Vsc   = sqrt(musc/Lsc);
Tsc   = Lsc/Vsc;
Asc   = Vsc/Tsc;
scale = diag([ones(3,1)./Lsc; ones(3,1)./Vsc]);

xp = scale*xp;
rp = xp(1:3);
vp = xp(4:6);
%% ---- Thrust and safe-miss parameters ------------------------------------
epsilon = 4e-2;
sigma   = 10/Lsc;

%% ---- Synthetic secondary state ------------------------------------------
% Place secondary with:
%   - relative velocity in the (alpha=40, beta=140) deg NTH direction
%   - relative position INSIDE the safe circle (|Dr0| = 0.5*sigma)
%   - small out-of-plane offset to make it non-trivial

alpha = 40*pi/180;  beta = 140*pi/180;
T_hat = vp/norm(vp);
H_hat = cross(rp,vp); H_hat = H_hat/norm(H_hat);
N_hat = cross(T_hat,H_hat);
R_nth_to_eci = [N_hat, T_hat, H_hat];

dv_NTH = [sin(alpha)*sin(beta); cos(beta); cos(alpha)*sin(beta)];
dv_NTH = dv_NTH / norm(dv_NTH);
dv_mag = 5e-4;                          % relative speed (normalised)
Dv0    = R_nth_to_eci * dv_NTH * dv_mag;

% Miss vector: 0.5*sigma in the b1 direction (inside the circle)
% b1 is along projection of vp onto B-plane
P_b    = eye(3) - (Dv0/norm(Dv0))*(Dv0/norm(Dv0))';
b1_raw = P_b * vp;
b1_hat_approx = b1_raw / norm(b1_raw);
Dr0    = 0*b1_hat_approx + (P_b * (rp/norm(rp)));
Dr0    = Dr0 / norm(Dr0) * 0.1 * sigma;   % ensure |Dr0| = 0.5*sigma exactly

xs = xp + [Dr0; Dv0];

fprintf('=== Synthetic conjunction ===\n');
fprintf('  |Dr0|         = %.6f (norm. length)\n', norm(Dr0));
fprintf('  |Dr0|/sigma   = %.4f\n', norm(Dr0)/sigma);
fprintf('  |Dv0|         = %.6f\n', norm(Dv0));

%% ---- Run optimal_cam ----------------------------------------------------
cam = optimal_cam(xp, xs, epsilon, sigma, -2*pi, ...
    'K',        360,    ...
    'Ntime',    50,   ...
    'AutoStop', true,   ...
    'FillTol',  0.02,   ...
    'Verbose',  true);

%% ---- Local helper: equinoctial -> cartesian (standalone) ----------------
function [r, v] = eq2cart(x)
a=x(1); p1=x(2); p2=x(3); q1=x(4); q2=x(5); L=x(6);
e=sqrt(p1^2+p2^2); wO=atan2(p1,p2); M=mod(L-wO,2*pi);
E=M; for k=1:50; dE=(M-E+e*sin(E))/(1-e*cos(E)+1e-300); E=E+dE; if abs(dE)<1e-13;break;end;end
nu=2*atan2(sqrt(1+e)*sin(E/2),sqrt(1-e)*cos(E/2)); F=nu+wO;
s2=1+q1^2+q2^2; psl=a*(1-e^2); rm=psl/(1+p1*sin(F)+p2*cos(F));
fh=(1/s2)*[1-q1^2+q2^2;2*q1*q2;-2*q1]; gh=(1/s2)*[2*q1*q2;1+q1^2-q2^2;2*q2];
r=rm*(cos(F)*fh+sin(F)*gh); hm=sqrt(psl);
v=(-(1/hm)*(p1+sin(F)))*fh+((1/hm)*(p2+cos(F)))*gh;
end
