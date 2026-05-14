function results = minimum_warning_time(y0, dv_dir, epsilon, sigma, t_min, varargin)
% MINIMUM_WARNING_TIME  Minimum warning time for low-thrust collision avoidance.
%
%   Implements the approximate solution for fast encounters from:
%   Dell'Elce, de Veld, Pomet — AAS 24-458, 2024.
%
%   Backward-integrates Eq.(21) from 2K+1 equally-spaced final positions on
%   the safe set circle, builds FFT wavefront interpolants (Eq.22), and
%   detects cut/conjugate loci.  Optionally stops automatically when the
%   entire disk of radius sigma is covered (AutoStop).
%
% -------------------------------------------------------------------------
% INPUTS
%   y0      (6x1)  Equinoctial elements of primary at nominal TCA, normalised
%                  (length unit = a, time unit = T_orb/(2*pi)).
%                  Order: [a, e*sin(w+O), e*cos(w+O), tan(i/2)*sin(w),
%                          tan(i/2)*cos(w), L=w+O+M], angles in RADIANS.
%   dv_dir  (3x1)  Unit vector Delta_v(z0)/||Delta_v|| in the
%                  Local-Normal/Tangential/H frame {N,T,H} of the primary:
%                  T = v/|v|,  H = r x v / |r x v|,  N = T x H.
%   epsilon scalar Normalised thrust-to-mass: eps = 4*pi^2*(T/m)*a^2/mu.
%   sigma   scalar Safe miss-distance threshold (normalised by a).
%   t_min   scalar Backward integration horizon (t_min < 0).
%                  Ignored when 'AutoStop' is true and the disk is filled
%                  before |t_min| is reached.
%
% OPTIONAL NAME-VALUE PAIRS
%   'K'        integer  Fourier modes; 2K+1 rays total. Default 360.
%   'Ntime'    integer  Number of integration time steps.  Default 2000.
%   'AutoStop' logical  Stop integration as soon as the entire disk is
%                       covered by the wavefront. Default false.
%   'FillTol'  scalar   Fraction of sigma used as "disk filled" threshold
%                       for AutoStop (default 0.02, i.e. 2% of sigma).
%   'Dr0'      (3x1)    If provided, also solves RED-SHOOT for this initial
%                       miss vector and annotates it on the figure.
%
% OUTPUTS
%   results  struct:
%     .t_vec      (Ntime_eff x 1) time vector actually integrated (0..t_stop)
%     .wf_b1      ((2K+1) x Ntime_eff)  B-plane b1 coord of wavefront / sigma
%     .wf_b2      ((2K+1) x Ntime_eff)  B-plane b2 coord of wavefront / sigma
%     .conj_locus cell(Ntime_eff,1) conjugate locus pts [Mx2] per time step
%     .cut_locus  cell(Ntime_eff,1) cut locus pts [Mx2] per time step
%     .jac_det    ((2K+1) x Ntime_eff) Jacobian determinant
%     .min_dist   (1 x Ntime_eff) min wavefront distance to origin / sigma
%     .t0_fill    scalar  t0 at which disk is fully covered (NaN if not reached)
%     .b1_hat     (3x1)  B-plane basis vector 1
%     .b2_hat     (3x1)  B-plane basis vector 2
%     .t0_query, .theta_query  RED-SHOOT solution (only if Dr0 provided)
%
% -------------------------------------------------------------------------
% NORMALISATION
%   Length : a  (semi-major axis of primary at nominal TCA)
%   Time   : T_orb/(2*pi)  =>  n0 = 1  (normalised mean motion)
%   mu     : 1  (normalised gravitational parameter)
% -------------------------------------------------------------------------

%% ---- Parse optional inputs -----------------------------------------------
p = inputParser;
addParameter(p, 'K',        360,   @(x) isscalar(x) && x>0 && floor(x)==x);
addParameter(p, 'Ntime',    2000,  @(x) isscalar(x) && x>0 && floor(x)==x);
addParameter(p, 'AutoStop', false, @(x) islogical(x) || x==0 || x==1);
addParameter(p, 'FillTol',  0.02,  @(x) isscalar(x) && x>0 && x<1);
addParameter(p, 'Dr0',      [],    @(x) isempty(x) || (isvector(x) && numel(x)==3));
parse(p, varargin{:});

K         = p.Results.K;
Ntime     = p.Results.Ntime;
auto_stop = logical(p.Results.AutoStop);
fill_tol  = p.Results.FillTol;
Dr0       = p.Results.Dr0;
if ~isempty(Dr0); Dr0 = Dr0(:); end

%% ---- Validate & normalise inputs -----------------------------------------
y0     = y0(:);
dv_dir = dv_dir(:);
assert(numel(y0)    == 6, 'y0 must be a 6-element vector (equinoctial elements).');
assert(numel(dv_dir) == 3, 'dv_dir must be a 3-element unit vector.');
assert(t_min < 0,          't_min must be negative (backward integration horizon).');
dv_dir = dv_dir / norm(dv_dir);

fprintf('=== Minimum Warning Time — Dell''Elce, de Veld, Pomet (AAS 24-458) ===\n');
fprintf('  epsilon = %.4e   sigma = %.4e   sigma/eps = %.4f\n', epsilon, sigma, sigma/epsilon);
fprintf('  t_min   = %.4f (norm.)   K = %d   Ntime = %d   AutoStop = %d\n\n', ...
        t_min, K, Ntime, auto_stop);

%% ---- Step 1: Cartesian state at nominal TCA ------------------------------
[r_y0, v_y0] = equinoctial_to_cartesian(y0);

%% ---- Step 2: B-plane basis {b1_hat, b2_hat} ------------------------------
%  B-plane is orthogonal to Delta_v direction in ECI.
%  b1 || projection of v(y0) onto B-plane (paper, below Eq. RED-SHOOT).
%  b2 = dv x b1  (right-hand, in B-plane).

dv_ECI = lntf_to_eci(dv_dir, r_y0, v_y0);
P_b    = eye(3) - dv_ECI * dv_ECI';      % projector onto B-plane

b1_hat = P_b * v_y0;
assert(norm(b1_hat) > 1e-10, ...
    'v(y0) is parallel to Delta_v: cannot define b1. Choose a different dv_dir or orbital configuration.');
b1_hat = b1_hat / norm(b1_hat);
b2_hat = cross(dv_ECI, b1_hat);
b2_hat = b2_hat / norm(b2_hat);

fprintf('B-plane basis:\n  b1 = [%+.6f %+.6f %+.6f]\n  b2 = [%+.6f %+.6f %+.6f]\n\n', ...
        b1_hat', b2_hat');

%% ---- Step 3: Jacobian dr/dy at y0 (3x6, numerical) ----------------------
drdy = drdy_equinoctial(y0);

%% ---- Step 4: Time grid and U(t) matrix (Eq. 20) -------------------------
%  U(t) = P_b * (dr/dy)(y0) * D(-t, y0)    [3x3]
%  D(tau, y0) = Phi(tau, phi^{-tau}(y0)) * G(phi^{-tau}(y0))  [6x3]
%
%  For the backward integration (t runs from 0 to t_min < 0):
%    tau = -t >= 0  (time interval from x = phi^{-tau}(y0) to y0)

t_vec = linspace(0, t_min, Ntime);   % [1 x Ntime], t_vec(1)=0, t_vec(end)=t_min
dt    = t_vec(2) - t_vec(1);          % dt < 0

fprintf('Pre-computing U(t) ... ');
U_all = zeros(3, 3, Ntime);           % U_all(:,:,k) = U(t_vec(k))
for k = 1:Ntime
    tau = -t_vec(k);                   % >= 0
    U_all(:,:,k) = P_b * drdy * D_matrix(tau, y0);
end
fprintf('done.\n');

%% ---- Step 5: Backward integration of 2K+1 extremal trajectories (Eq.21) -

% Eq.(21): for each ray k with FIXED final direction n_k = b1*cos(theta_k)+b2*sin(theta_k),
% d(Dr_k/sigma)/dt = (epsilon/sigma) * U(t) * [U(t)' * n_k] / ||U(t)' * n_k||

% KEY: the RHS uses n_k (the FIXED boundary direction at t=0), NOT the current
% position Dr_norm. The control u*(t) = U(t)'*n_k / ||...|| is completely
% determined by the costate n_k which is constant along each ray (Lemma 2).

% Starting condition: Dr_k(t=0)/sigma = n_k.
% dt < 0 => Euler step moves backward in time.
N_rays = 2*K + 1;
theta_k = (2*pi * (-K:K)') / N_rays; % [N_rays x 1]

% n_hat_all [3 x N_rays]: FIXED final directions. Constant for each ray.
n_hat_all = b1_hat * cos(theta_k') + b2_hat * sin(theta_k'); % [3 x N_rays]

fprintf('Backward integration (%d rays) ... ', N_rays);

% Pre-allocate storage
wf_b1_all = zeros(N_rays, Ntime);
wf_b2_all = zeros(N_rays, Ntime);
min_dist_all = zeros(1, Ntime);

% Running state [3 x N_rays]: Delta_r / sigma, starts on the boundary.
Dr_norm = n_hat_all;

wf_b1_all(:,1) = b1_hat' * Dr_norm;
wf_b2_all(:,1) = b2_hat' * Dr_norm;
min_dist_all(1) = min(sqrt(wf_b1_all(:,1).^2 + wf_b2_all(:,1).^2));

t_stop_idx = Ntime;
t0_fill = NaN;

for k = 2:Ntime
    U = U_all(:,:,k); % U at current time step
    
    % RHS: uses fixed n_hat_all (costate), NOT the current Dr_norm
    Utn = U' * n_hat_all; % [3 x N_rays]
    Utn_norms = sqrt(sum(Utn.^2, 1)); % [1 x N_rays]
    
    valid = Utn_norms > 1e-14;
    rhs = zeros(3, N_rays);
    
    % CORRECTED LINE: Added the (epsilon / sigma) scaling factor
    rhs(:,valid) = (epsilon / sigma) * U * (Utn(:,valid) ./ Utn_norms(valid));
    
    Dr_norm = Dr_norm + dt * rhs; % Euler: dt<0 => backward

    wf_b1_all(:,k) = b1_hat' * Dr_norm;
    wf_b2_all(:,k) = b2_hat' * Dr_norm;
    min_dist_all(k) = min(sqrt(wf_b1_all(:,k).^2 + wf_b2_all(:,k).^2));

    if auto_stop && min_dist_all(k) < fill_tol
        t0_fill = t_vec(k);
        t_stop_idx = k;
        fprintf('\n AutoStop triggered at t0 = %.4f (disk filled, min_dist/sigma = %.4f)', ...
            t0_fill, min_dist_all(k));
        break;
    end
end
fprintf(' done.\n\n');

% Truncate to actual integrated range
t_vec     = t_vec(1:t_stop_idx);
wf_b1     = wf_b1_all(:, 1:t_stop_idx);
wf_b2     = wf_b2_all(:, 1:t_stop_idx);
min_dist  = min_dist_all(1:t_stop_idx);
Ntime_eff = t_stop_idx;

if isnan(t0_fill)
    % AutoStop not triggered; check final min_dist anyway
    idx_fill = find(min_dist < fill_tol, 1, 'first');
    if ~isempty(idx_fill)
        t0_fill = t_vec(idx_fill);
        fprintf('Disk filled at t0 = %.4f  (|t0| = %.4f orb. periods)\n\n', ...
                t0_fill, abs(t0_fill)/(2*pi));
    else
        fprintf('WARNING: Disk not fully covered within integration window.\n');
        fprintf('  Consider decreasing t_min or enabling AutoStop.\n\n');
    end
end

%% ---- Step 6: FFT derivatives -> Jacobian determinant -> loci -------------
%  Conjugate locus: det[d(c1,c2)/d(t0,theta)] = 0  (Eq. 23, paper)
%    = (dc1/dt0)(dc2/dtheta) - (dc1/dtheta)(dc2/dt0) = 0
%  Cut locus: self-intersections of the wavefront polygon.

fprintf('Computing cut and conjugate loci ... ');

% d/dtheta via spectral differentiation in Fourier space
fft_modes = ifftshift(-K:K)';                  % [N_rays x 1]

dwf_b1_dth = zeros(N_rays, Ntime_eff);
dwf_b2_dth = zeros(N_rays, Ntime_eff);
for k = 1:Ntime_eff
    dwf_b1_dth(:,k) = real(ifft(fft(wf_b1(:,k)) .* (1i * fft_modes)));
    dwf_b2_dth(:,k) = real(ifft(fft(wf_b2(:,k)) .* (1i * fft_modes)));
end

% d/dt0 via central finite differences
dwf_b1_dt = central_diff(wf_b1, dt, 2);   % diff along dim 2 (time)
dwf_b2_dt = central_diff(wf_b2, dt, 2);

% Jacobian determinant [N_rays x Ntime_eff]
jac_det = dwf_b1_dt .* dwf_b2_dth - dwf_b1_dth .* dwf_b2_dt;

% Conjugate locus: sign changes of jac_det along theta (for each time)
conj_locus = cell(Ntime_eff, 1);
cut_locus  = cell(Ntime_eff, 1);

for k = 2:Ntime_eff
    jd  = jac_det(:, k);
    c1k = wf_b1(:, k);
    c2k = wf_b2(:, k);

    % -- Conjugate locus --
    sc = sign(jd); sc(sc==0) = 1;
    idx_zc = find(diff(sc) ~= 0);
    cpts = zeros(numel(idx_zc), 2);
    for m = 1:numel(idx_zc)
        ii  = idx_zc(m);
        alp = jd(ii) / (jd(ii) - jd(ii+1));
        cpts(m,1) = c1k(ii) + alp*(c1k(ii+1) - c1k(ii));
        cpts(m,2) = c2k(ii) + alp*(c2k(ii+1) - c2k(ii));
    end
    conj_locus{k} = cpts;

    % -- Cut locus (self-intersections of wavefront polygon) --
    cut_locus{k} = find_self_intersections(c1k, c2k);
end

fprintf('done.\n\n');

%% ---- Step 7 (optional): RED-SHOOT for query point -----------------------
results.t0_query    = [];
results.theta_query = [];
if ~isempty(Dr0)
    [t0_q, theta_q] = solve_red_shoot(Dr0, sigma, b1_hat, b2_hat, ...
                                       wf_b1, wf_b2, t_vec, theta_k);
    results.t0_query    = t0_q;
    results.theta_query = theta_q;
    fprintf('RED-SHOOT for Dr0 = [%.4f %.4f %.4f]:\n', Dr0);
    fprintf('  t0 = %.4f   theta = %.4f rad (%.1f deg)\n\n', t0_q, theta_q, theta_q*180/pi);
end

%% ---- Pack outputs --------------------------------------------------------
results.t_vec      = t_vec(:);
results.wf_b1      = wf_b1;
results.wf_b2      = wf_b2;
results.conj_locus = conj_locus;
results.cut_locus  = cut_locus;
results.jac_det    = jac_det;
results.min_dist   = min_dist;
results.t0_fill    = t0_fill;
results.b1_hat     = b1_hat;
results.b2_hat     = b2_hat;

%% ---- Step 8: Generate figure (analog to paper Fig. 2) -------------------
generate_figure(results, sigma, epsilon, theta_k, auto_stop);

end  % ===== END main function ===============================================


% =============================================================================
%                             LOCAL FUNCTIONS
% =============================================================================

% -----------------------------------------------------------------------------
function [r, v] = equinoctial_to_cartesian(x)
% Equinoctial elements -> Cartesian ECI (mu=1, normalised).
% x = [a, p1=e*sin(w+O), p2=e*cos(w+O), q1=tan(i/2)*sin(w),
%       q2=tan(i/2)*cos(w), L=w+O+M]  (L = mean longitude, angles in rad).

a  = x(1);  p1 = x(2);  p2 = x(3);
q1 = x(4);  q2 = x(5);  L  = x(6);

e        = sqrt(p1^2 + p2^2);
w_plus_O = atan2(p1, p2);
M_anom   = mod(L - w_plus_O, 2*pi);
E_anom   = kepler_newton(M_anom, e);
nu       = 2*atan2(sqrt(1+e)*sin(E_anom/2), sqrt(1-e)*cos(E_anom/2));
F        = nu + w_plus_O;       % true longitude

s2    = 1 + q1^2 + q2^2;
p_sl  = a * (1 - e^2);         % semi-latus rectum
r_mag = p_sl / (1 + p1*sin(F) + p2*cos(F));

f_hat = (1/s2) * [ 1 - q1^2 + q2^2;  2*q1*q2;  -2*q1 ];
g_hat = (1/s2) * [ 2*q1*q2;  1 + q1^2 - q2^2;   2*q2 ];

r = r_mag*(cos(F)*f_hat + sin(F)*g_hat);

h_mag  = sqrt(p_sl);           % mu=1
X_dot  = -(1/h_mag)*(p1 + sin(F));
Y_dot  =  (1/h_mag)*(p2 + cos(F));
v      = X_dot*f_hat + Y_dot*g_hat;
end

% -----------------------------------------------------------------------------
function E = kepler_newton(M, e)
% Newton solver for Kepler's equation  M = E - e*sin(E).
E = M;
for k = 1:50
    dE = (M - E + e*sin(E)) / (1 - e*cos(E) + 1e-300);
    E  = E + dE;
    if abs(dE) < 1e-13; break; end
end
end

% -----------------------------------------------------------------------------
function drdy = drdy_equinoctial(y0)
% Numerical Jacobian  dr/dy  [3x6]  via central finite differences.
h = 1e-6;
drdy = zeros(3, 6);
for j = 1:6
    yp = y0; yp(j) = yp(j) + h;
    ym = y0; ym(j) = ym(j) - h;
    [rp,~] = equinoctial_to_cartesian(yp);
    [rm,~] = equinoctial_to_cartesian(ym);
    drdy(:,j) = (rp - rm) / (2*h);
end
end

% -----------------------------------------------------------------------------
function G = gve_equinoctial(x_true)
% Gauss Variational Equations for equinoctial elements [6x3].
% Thrust components in RTN (radial, transverse, normal) frame.
% Reference: Battin (1999) Ch.10 pp.492-493 / Schaub & Junkins (2018).
%
% x_true = [a, p1, p2, q1, q2, F]  where F is the TRUE longitude.

a  = x_true(1);  p1 = x_true(2);  p2 = x_true(3);
q1 = x_true(4);  q2 = x_true(5);  F  = x_true(6);

s2 = 1 + q1^2 + q2^2;
e  = sqrt(p1^2 + p2^2);
p  = a*(1 - e^2);              % semi-latus rectum
h  = sqrt(p);                  % mu=1  =>  h = sqrt(mu*p) = sqrt(p)
w  = 1 + p1*sin(F) + p2*cos(F);   % r = p/w

% --- Radial (R) ---
G_R    = zeros(6,1);
G_R(1) =  2*a^2*(p2*sin(F) - p1*cos(F)) / h;
G_R(2) =  sin(F) / h;
G_R(3) = -cos(F) / h;
G_R(4) =  0;
G_R(5) =  0;
G_R(6) = -(q1*sin(F) - q2*cos(F)) / (h*w);

% --- Transverse (T) ---
G_T    = zeros(6,1);
G_T(1) =  2*a^2*w / h;
G_T(2) =  ((w+1)*cos(F) + p1) / h;
G_T(3) = -((w+1)*sin(F) - p2) / h;
G_T(4) =  0;
G_T(5) =  0;
G_T(6) =  (q1*cos(F) + q2*sin(F)) / (h*w);

% --- Normal (N) ---
G_N    = zeros(6,1);
G_N(1) =  0;
G_N(2) =  0;
G_N(3) =  0;
G_N(4) =  (s2*sin(F)) / (2*h);
G_N(5) =  (s2*cos(F)) / (2*h);  % note: this should be negative for some conventions
G_N(6) = -(s2*(q1*sin(F) - q2*cos(F))) / (2*h*w);  % note sign

G = [G_R, G_T, G_N];
end

% -----------------------------------------------------------------------------
function x = keplerian_flow(y0, t)
% Keplerian flow phi^t(y0): only mean longitude L advances.
x    = y0;
x(6) = y0(6) + t;    % n0 = 1 (normalised)
end

% -----------------------------------------------------------------------------
function Phi = keplerian_STM(t, x)
% Analytic Keplerian STM for equinoctial elements [6x6] (Example 1, paper).
% Only off-diagonal entry: Phi(6,1) = dn/da * t = -3/2 * a^{-5/2} * t.
Phi      = eye(6);
a        = x(1);
Phi(6,1) = -1.5 * t * a^(-2.5);
end

% -----------------------------------------------------------------------------
function [F_true] = mean_long_to_true(x)
% Extract true longitude F from equinoctial state x with mean longitude L.
p1 = x(2); p2 = x(3); L = x(6);
e        = sqrt(p1^2 + p2^2);
w_plus_O = atan2(p1, p2);
M_anom   = mod(L - w_plus_O, 2*pi);
E_anom   = kepler_newton(M_anom, e);
nu       = 2*atan2(sqrt(1+e)*sin(E_anom/2), sqrt(1-e)*cos(E_anom/2));
F_true   = nu + w_plus_O;
end

% -----------------------------------------------------------------------------
function D = D_matrix(tau, y0)
% D(tau, y0) = Phi(tau, phi^{-tau}(y0)) * G(phi^{-tau}(y0))  [6x3]
% (Definition after Eq.(5) in paper.)
x       = keplerian_flow(y0, -tau);     % phi^{-tau}(y0)
Phi     = keplerian_STM(tau, x);        % [6x6]
F_true  = mean_long_to_true(x);
x_true  = x;  x_true(6) = F_true;
G       = gve_equinoctial(x_true);      % [6x3]
D       = Phi * G;
end

% -----------------------------------------------------------------------------
function dv_ECI = lntf_to_eci(dv_lntf, r_ECI, v_ECI)
% Convert vector from Local-Normal/Tangential/H frame to ECI.
% Paper (footnote p.15): T=v/|v|, H=rxv/|rxv|, N=TxH.
T_hat = v_ECI / norm(v_ECI);
H_hat = cross(r_ECI, v_ECI);  H_hat = H_hat / norm(H_hat);
N_hat = cross(T_hat, H_hat);
R_to_ECI = [N_hat, T_hat, H_hat];   % columns: N, T, H
dv_ECI   = R_to_ECI * dv_lntf;
dv_ECI   = dv_ECI / norm(dv_ECI);
end

% -----------------------------------------------------------------------------
function df = central_diff(f, d, dim)
% Central finite differences along dimension dim.
df = zeros(size(f));
idx_start = repmat({':'}, 1, ndims(f));
idx_end   = idx_start;
idx_mid_p = idx_start;
idx_mid_m = idx_start;
N = size(f, dim);

s = idx_start; s{dim} = 1;       df_s = idx_start; df_s{dim} = 1;
e_i = idx_start; e_i{dim} = N;   df_e = idx_start; df_e{dim} = N;
% Forward diff at start
sp = idx_start; sp{dim} = 2;
df(df_s{:}) = (f(sp{:}) - f(s{:})) / d;
% Backward diff at end
em = idx_start; em{dim} = N-1;
df(df_e{:}) = (f(e_i{:}) - f(em{:})) / d;
% Central diff in interior
for k = 2:N-1
    ip = idx_start; ip{dim} = k+1;
    im = idx_start; im{dim} = k-1;
    ik = idx_start; ik{dim} = k;
    df(ik{:}) = (f(ip{:}) - f(im{:})) / (2*d);
end
end

% -----------------------------------------------------------------------------
function pts = find_self_intersections(c1, c2)
% Self-intersections of closed polygon (c1,c2) using segment-segment test.
% Returns [Mx2] array of intersection B-plane coordinates.
N   = numel(c1);
pts = [];
for i = 1:N
    i2 = mod(i, N) + 1;
    P1 = [c1(i); c2(i)];   P2 = [c1(i2); c2(i2)];
    for j = i+2:N
        if j == N && i == 1; continue; end
        j2 = mod(j, N) + 1;
        Q1 = [c1(j); c2(j)];  Q2 = [c1(j2); c2(j2)];
        pt = seg_intersect(P1, P2, Q1, Q2);
        if ~isempty(pt); pts = [pts; pt']; end %#ok<AGROW>
    end
end
end

% -----------------------------------------------------------------------------
function pt = seg_intersect(P1, P2, Q1, Q2)
% Intersection of segments P1-P2 and Q1-Q2. Returns [] if none.
pt    = [];
d1    = P2 - P1;  d2 = Q2 - Q1;
denom = d1(1)*d2(2) - d1(2)*d2(1);
if abs(denom) < 1e-14; return; end
t = ((Q1(1)-P1(1))*d2(2) - (Q1(2)-P1(2))*d2(1)) / denom;
s = ((Q1(1)-P1(1))*d1(2) - (Q1(2)-P1(2))*d1(1)) / denom;
if t >= 0 && t <= 1 && s >= 0 && s <= 1
    pt = P1 + t*d1;
end
end

% -----------------------------------------------------------------------------
function [t0_q, theta_q] = solve_red_shoot(Dr0, sigma, b1_hat, b2_hat, ...
                                             wf_b1, wf_b2, t_vec, theta_k)
% Solve RED-SHOOT: find (t0,theta) such that wavefront passes through Dr0.
% Scans all (t_vec, theta_k) and returns the pair with minimum distance.
c1_tgt = (b1_hat' * Dr0) / sigma;
c2_tgt = (b2_hat' * Dr0) / sigma;
Ntime  = numel(t_vec);
best   = Inf;
t0_q   = t_vec(1);  theta_q = theta_k(1);
for k = 1:Ntime
    d2 = (wf_b1(:,k) - c1_tgt).^2 + (wf_b2(:,k) - c2_tgt).^2;
    [md, idx] = min(d2);
    if md < best
        best    = md;
        t0_q    = t_vec(k);
        theta_q = theta_k(idx);
    end
end
end

% =============================================================================
%                          FIGURE (analog to paper Fig.2)
% =============================================================================
function generate_figure(results, sigma, epsilon, theta_k, auto_stop)

t_vec      = results.t_vec;
wf_b1      = results.wf_b1;
wf_b2      = results.wf_b2;
conj_locus = results.conj_locus;
cut_locus  = results.cut_locus;
min_dist   = results.min_dist;
t0_fill    = results.t0_fill;
Ntime_eff  = numel(t_vec);
N_rays     = size(wf_b1, 1);

% Choose 4 representative time snapshots evenly distributed,
% with the last snapshot at the disk-fill time (if known) or at t_min.
if ~isnan(t0_fill)
    [~, idx_fill] = min(abs(t_vec - t0_fill));
else
    idx_fill = Ntime_eff;
end
t_snaps = unique(round(linspace(2, idx_fill, 4)));
if numel(t_snaps) < 4
    t_snaps = round(linspace(2, idx_fill, 4));
end
panel_lbl = {'(a)','(b)','(c)','(d)'};
th_circ = linspace(0, 2*pi, 500);
cmap    = parula(N_rays);

fig = figure('Name','Minimum Warning Time — Wavefront Evolution', ...
             'Position', [60, 60, 960, 900]);

for ip = 1:4
    ax = subplot(2, 2, ip);
    k  = t_snaps(ip);

    % ---- Safe-set circle (unit circle in normalised B-plane coords) ----
    fill(ax, cos(th_circ), sin(th_circ), [0.92 0.92 0.92], ...
         'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 1.5, 'FaceAlpha', 0.55); hold on;

    % ---- Extremal trajectories (thinned for readability) ----
    step = max(1, floor(N_rays / 80));
    for ir = 1:step:N_rays
        plot(ax, wf_b1(ir, 1:k), wf_b2(ir, 1:k), '-', ...
             'Color', [cmap(ir,:), 0.55], 'LineWidth', 0.7);
    end

    % ---- Wavefront at this snapshot ----
    wf_c1_closed = [wf_b1(:,k); wf_b1(1,k)];
    wf_c2_closed = [wf_b2(:,k); wf_b2(1,k)];
    plot(ax, wf_c1_closed, wf_c2_closed, 'r-', 'LineWidth', 2.0, ...
         'DisplayName','Wavefront');

    % ---- Conjugate locus ----
    if ~isempty(conj_locus{k})
        scatter(ax, conj_locus{k}(:,1), conj_locus{k}(:,2), 36, ...
                'g', 'filled', 'DisplayName','Conjugate locus');
    end

    % ---- Cut locus ----
    if ~isempty(cut_locus{k}) && size(cut_locus{k},1) > 0
        scatter(ax, cut_locus{k}(:,1), cut_locus{k}(:,2), 36, ...
                'b', 'filled', 'DisplayName','Cut locus');
    end

    % ---- Query point (RED-SHOOT) if available ----
    if ~isempty(results.t0_query) && t_vec(k) <= results.t0_query + eps
        scatter(ax, (results.b1_hat' * results.t0_query) / sigma, ...
                    (results.b2_hat' * results.t0_query) / sigma, ...
                60, 'm', 'filled', 'DisplayName','Dr0 (query)');
    end

    axis(ax, 'equal'); grid(ax, 'on');
    xlim(ax, [-1.35, 1.35]); ylim(ax, [-1.35, 1.35]);
    xlabel(ax, '\Deltar \cdot \bfb_1 / \sigma', 'FontSize', 10);
    ylabel(ax, '\Deltar \cdot \bfb_2 / \sigma', 'FontSize', 10);

    t_label = t_vec(k);
    if ~isnan(t0_fill) && k == t_snaps(end)
        title(ax, sprintf('%s t_0 = %.3f  [DISK FILLED]', panel_lbl{ip}, t_label), ...
              'FontSize', 10, 'Color', [0.7 0 0]);
    else
        title(ax, sprintf('%s t_0 = %.3f (norm.)', panel_lbl{ip}, t_label), 'FontSize', 10);
    end

    % Legend only on first panel
    if ip == 1
        legend(ax, 'Location', 'northeast', 'FontSize', 7);
    end
    hold off;
end
end
