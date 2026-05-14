function cam = optimal_cam(xp, xs, epsilon, sigma, t_min, varargin)
% OPTIMAL_CAM  Optimal low-thrust CAM for a given conjunction geometry.
%
%   Given the Cartesian states of the primary and secondary at nominal TCA,
%   finds the optimal CONSTANT costate lambda* and the resulting TIME-VARYING
%   thrust direction profile u*(t) from maneuver start t0* to TCA (t=0).
%
%   From the paper (Eq. 13 + Lemma 2 / Section "Approximate solution"):
%     - The linearised system is DRIFTLESS, so the costate lambda is CONSTANT.
%     - But u*(t) = B(t,z0)' lambda / ||B(t,z0)' lambda|| IS time-varying
%       because B(t,z0) = U(t) depends on the STM and GVE evaluated along
%       the nominal trajectory.
%     - Specifically, d(Dr)/dt = U(t) u*(t) * epsilon  (Eq. 20),
%       where U(t) = P_b * (dr/dy)(y0) * D(-t, y0).
%     - The costate is lambda = P_b * n_hat(theta*)  (from transversality).
%     - Therefore: u*(t) = U(t)' * (P_b * n_hat(theta*))
%                          / ||U(t)' * (P_b * n_hat(theta*))||
%       which rotates continuously in RTN as t goes from t0* to 0.
%
%   The function:
%     1. Converts primary ECI state -> equinoctial elements y0.
%     2. Extracts relative velocity direction dv_dir in NTH frame.
%     3. Calls minimum_warning_time() to build the B-plane wavefront.
%     4. Solves RED-SHOOT (Newton refinement) to find (t0*, theta*).
%     5. Evaluates the FULL TIME HISTORY of u*(t) for t in [t0*, 0].
%     6. Produces a summary figure with:
%          Panel 1: B-plane (wavefront + Dr0 + optimal ray)
%          Panel 2: RTN components of u*(t) vs time
%          Panel 3: Thrust direction angles (azimuth + elevation) vs time
%
% -------------------------------------------------------------------------
% INPUTS
%   rp_ECI  (3x1)  Primary position at nominal TCA  [normalised]
%   vp_ECI  (3x1)  Primary velocity at nominal TCA  [normalised]
%   rs_ECI  (3x1)  Secondary position at nominal TCA [normalised]
%   vs_ECI  (3x1)  Secondary velocity at nominal TCA [normalised]
%   epsilon scalar Normalised thrust-to-mass ratio
%   sigma   scalar Safe miss-distance threshold (normalised)
%   t_min   scalar Backward integration horizon (< 0). Used as hard limit;
%                  use 'AutoStop',true for automatic termination.
%
% OPTIONAL NAME-VALUE PAIRS
%   'K'        integer  Fourier modes (default 360)
%   'Ntime'    integer  Integration time steps (default 2000)
%   'AutoStop' logical  Stop when disk is filled (default true)
%   'FillTol'  scalar   Disk-fill threshold (default 0.02)
%   'Nprofile' integer  Number of time points for the control profile
%                       output, t in [t0*, 0]. Default 500.
%   'Verbose'  logical  Print detailed output (default true)
%
% -------------------------------------------------------------------------
% OUTPUTS  cam  struct:
%
%   Conjunction geometry
%   .Dr0_ECI      (3x1)    Initial miss vector (ECI)
%   .Dr0_norm     scalar   |Dr0|/sigma
%   .dv_ECI       (3x1)    Unit relative velocity (ECI)
%   .b1_hat       (3x1)    B-plane basis vector 1 (ECI)
%   .b2_hat       (3x1)    B-plane basis vector 2 (ECI)
%   .Dr0_b1       scalar   Dr0.b1 / sigma
%   .Dr0_b2       scalar   Dr0.b2 / sigma
%
%   Optimal maneuver summary
%   .t0_star      scalar   Optimal maneuver start time (normalised, <0)
%   .t0_star_orb  scalar   |t0*| in orbital periods
%   .theta_star   scalar   Optimal boundary angle [rad]
%   .lambda_star  (3x1)    Constant costate = P_b * n_hat(theta*) (ECI)
%
%   Time-varying control profile  (Nprofile points, t in [t0*,0])
%   .t_profile    (Nprofile x 1)  Time vector [t0*, ..., 0]
%   .u_RTN        (3 x Nprofile)  Optimal thrust direction in RTN of primary
%                                 Rows: [R; T; N]
%   .u_ECI        (3 x Nprofile)  Same in ECI
%   .azimuth_deg  (1 x Nprofile)  Angle from T toward R in R-T plane [deg]
%   .elevation_deg(1 x Nprofile)  Out-of-plane angle from R-T plane [deg]
%   .u_RTN_start  (3x1)  u*(t0*)  — thrust direction at maneuver start
%   .u_RTN_end    (3x1)  u*(0)    — thrust direction at TCA
%
%   Wavefront intermediate results
%   .wf           struct  Full output of minimum_warning_time()
%   .y0           (6x1)   Primary equinoctial elements at TCA
%
% -------------------------------------------------------------------------

%% ---- Parse inputs --------------------------------------------------------
p = inputParser;
addParameter(p,'K',        360,  @(x) isscalar(x)&&x>0);
addParameter(p,'Ntime',    2000, @(x) isscalar(x)&&x>0);
addParameter(p,'AutoStop', true, @(x) islogical(x)||x==0||x==1);
addParameter(p,'FillTol',  0.02, @(x) isscalar(x)&&x>0&&x<1);
addParameter(p,'Nprofile', 500,  @(x) isscalar(x)&&x>1);
addParameter(p,'Verbose',  true, @(x) islogical(x)||x==0||x==1);
parse(p, varargin{:});
opt = p.Results;
verbose  = logical(opt.Verbose);
Nprofile = opt.Nprofile;

%% ---- Column vectors ------------------------------------------------------
rp = xp(1:3); vp = xp(4:6);
rs = xs(1:3); vs = xs(4:6);

%% ---- Conjunction geometry ------------------------------------------------
Dr0     = rs - rp;
Dv0     = vs - vp;
Dr0_mag = norm(Dr0);
Dv0_mag = norm(Dv0);
assert(Dv0_mag > 1e-14, 'Relative velocity is zero.');
if Dr0_mag >= sigma
    warning('optimal_cam:notConjunction', ...
        '|Dr0| = %.4f >= sigma = %.4f. Not strictly a conjunction.', ...
        Dr0_mag, sigma);
end
dv_unit_ECI = Dv0 / Dv0_mag;

%% ---- Primary equinoctial elements ----------------------------------------
y0 = cartesian_to_equinoctial(rp, vp);

%% ---- dv_dir in NTH frame -------------------------------------------------
T_hat = vp/norm(vp);
H_hat = cross(rp,vp); H_hat = H_hat/norm(H_hat);
N_hat = cross(T_hat,H_hat);
R_ECI_to_NTH = [N_hat, T_hat, H_hat]';
dv_NTH = R_ECI_to_NTH * dv_unit_ECI;

if verbose
    fprintf('\n=== OPTIMAL CAM ===\n');
    fprintf('  |Dr0|/sigma = %.4f  (%s)\n', Dr0_mag/sigma, ...
        ternary(Dr0_mag<sigma,'CONJUNCTION','outside safe circle'));
    fprintf('  epsilon = %.4e,  sigma = %.4e\n\n', epsilon, sigma);
end

%% ---- Build wavefront -------------------------------------------------------
if verbose; fprintf('--- Building wavefront ---\n'); end
wf = minimum_warning_time(y0, dv_NTH, epsilon, sigma, t_min, ...
    'K',        opt.K,        ...
    'Ntime',    opt.Ntime,    ...
    'AutoStop', opt.AutoStop, ...
    'FillTol',  opt.FillTol,  ...
    'Dr0',      Dr0);

b1_hat = wf.b1_hat;
b2_hat = wf.b2_hat;
t_vec  = wf.t_vec;          % [Ntime_eff x 1], from 0 to t_stop (<=t_min)
wf_b1  = wf.wf_b1;
wf_b2  = wf.wf_b2;
K      = opt.K;
theta_k = (2*pi*(-(K):(K))')/(2*K+1);

%% ---- RED-SHOOT: find (t0*, theta*) ----------------------------------------
c1_tgt = (b1_hat' * Dr0) / sigma;
c2_tgt = (b2_hat' * Dr0) / sigma;
N_rays  = numel(theta_k);
Ntime_eff = numel(t_vec);

% Coarse grid search
dist2 = (wf_b1 - c1_tgt).^2 + (wf_b2 - c2_tgt).^2;
[~, lin_idx] = min(dist2(:));
[ir_best, kt_best] = ind2sub(size(dist2), lin_idx);

% Newton refinement
[t0_star, theta_star] = refine_red_shoot( ...
    c1_tgt, c2_tgt, wf_b1, wf_b2, t_vec, theta_k, kt_best, ir_best);

if verbose
    fprintf('\nRED-SHOOT result:\n');
    fprintf('  t0*     = %.6f (norm.)  = %.4f orbital periods\n', ...
            t0_star, abs(t0_star)/(2*pi));
    fprintf('  theta*  = %.6f rad  (%.2f deg)\n\n', theta_star, theta_star*180/pi);
end

%% ---- Constant costate lambda* (ECI, in B-plane) --------------------------
%  From transversality (Eq. 15 + 17):  lambda = -P_b * n_hat(theta*)
%  The sign convention is chosen so that u* exits the unsafe disk, i.e.
%  n . d(Dr)/dt > 0.  (Paper, after Eq. 20: n U(t) n > 0 for all t < 0)
P_b      = eye(3) - (dv_unit_ECI * dv_unit_ECI');
n_star   = b1_hat*cos(theta_star) + b2_hat*sin(theta_star);  % [3x1] ECI
lambda   = P_b * n_star;
lambda   = lambda / norm(lambda);

%% ---- dr/dy at y0 (3x6) ---------------------------------------------------
drdy = drdy_equinoctial(y0);

%% ---- Full time-varying control profile u*(t), t in [t0*, 0] --------------
%  From Eq. (20) and PMP Eq. (13):
%    u*(t) = U(t)' * lambda / || U(t)' * lambda ||
%  where U(t) = P_b * drdy * D(-t, y0)   [3x3]
%
%  Note: U(t)' * lambda = D(-t,y0)' * drdy' * P_b * lambda  [3x1]  in RTN.
%  This is the "costate in GVE input space" and rotates as the primary moves
%  along its orbit from t0* to TCA.

t_profile = linspace(t0_star, 0, Nprofile)';   % [Nprofile x 1]

u_ECI_all    = zeros(3, Nprofile);
u_RTN_all    = zeros(3, Nprofile);
azimuth_all  = zeros(1, Nprofile);
elevation_all= zeros(1, Nprofile);

if verbose; fprintf('Computing control profile over %d time steps ... ', Nprofile); end
for k = 1:Nprofile
    t = t_profile(k);
    tau = -t;                       % tau >= 0

    % U(t) = P_b * drdy * D(tau, y0)  [3x3]
    D   = D_matrix_local(tau, y0);  % [6x3]
    U_t = P_b * drdy * D;           % [3x3]

    % Costate in GVE-input (RTN) space
    psi = U_t' * lambda;            % [3x1]
    psi_norm = norm(psi);
    if psi_norm < 1e-14
        % degenerate (conjugate point neighbourhood): keep previous direction
        if k > 1
            u_RTN_all(:,k)  = u_RTN_all(:,k-1);
            u_ECI_all(:,k)  = u_ECI_all(:,k-1);
        end
        azimuth_all(k)   = azimuth_all(max(1,k-1));
        elevation_all(k) = elevation_all(max(1,k-1));
        continue;
    end
    u_rtn = psi / psi_norm;         % unit thrust direction in RTN

    % RTN basis of primary at time t (propagate along Keplerian orbit)
    y_t   = keplerian_flow_local(y0, t);
    [r_t, v_t] = equinoctial_to_cartesian_local(y_t);
    R_hat = r_t / norm(r_t);
    H_t   = cross(r_t, v_t); H_t = H_t / norm(H_t);
    T_hat_t = cross(H_t, R_hat);
    R_RTN_to_ECI = [R_hat, T_hat_t, H_t];   % columns: R, T, N

    u_eci = R_RTN_to_ECI * u_rtn;

    u_RTN_all(:,k)   = u_rtn;
    u_ECI_all(:,k)   = u_eci;
    azimuth_all(k)   = atan2d(u_rtn(1), u_rtn(2));   % from T toward R
    elevation_all(k) = asind(u_rtn(3));               % out-of-plane
end
if verbose; fprintf('done.\n\n'); end

%% ---- Print summary --------------------------------------------------------
if verbose
    fprintf('=== OPTIMAL MANEUVER SUMMARY ===\n');
    fprintf('  t0*           = %.6f (norm.)  = %.4f orb. periods\n', ...
            t0_star, abs(t0_star)/(2*pi));
    fprintf('  theta*        = %.4f rad  (%.2f deg)\n', theta_star, theta_star*180/pi);
    fprintf('  u*(t0*) RTN   = [R=%+.4f, T=%+.4f, N=%+.4f]\n', u_RTN_all(:,1));
    fprintf('  u*(0)   RTN   = [R=%+.4f, T=%+.4f, N=%+.4f]\n', u_RTN_all(:,end));
    fprintf('  Azimuth range : %.2f to %.2f deg\n', ...
            min(azimuth_all), max(azimuth_all));
    fprintf('  Elevation range: %.2f to %.2f deg\n\n', ...
            min(elevation_all), max(elevation_all));
end

%% ---- Pack outputs ---------------------------------------------------------
cam.Dr0_ECI       = Dr0;
cam.Dr0_norm      = Dr0_mag/sigma;
cam.dv_ECI        = dv_unit_ECI;
cam.b1_hat        = b1_hat;
cam.b2_hat        = b2_hat;
cam.Dr0_b1        = c1_tgt;
cam.Dr0_b2        = c2_tgt;
cam.t0_star       = t0_star;
cam.t0_star_orb   = abs(t0_star)/(2*pi);
cam.theta_star    = theta_star;
cam.lambda_star   = lambda;
cam.t_profile     = t_profile;
cam.u_RTN         = u_RTN_all;
cam.u_ECI         = u_ECI_all;
cam.azimuth_deg   = azimuth_all;
cam.elevation_deg = elevation_all;
cam.u_RTN_start   = u_RTN_all(:,1);
cam.u_RTN_end     = u_RTN_all(:,end);
cam.wf            = wf;
cam.y0            = y0;

%% ---- Figure ---------------------------------------------------------------
generate_cam_figure(cam, sigma, epsilon);

end   % ===== END main function ==============================================


% =============================================================================
%                           LOCAL FUNCTIONS
% =============================================================================

function y = cartesian_to_equinoctial(r, v)
rv=norm(r); vv=norm(v); hv=cross(r,v); h=norm(hv); H_hat=hv/h;
e_vec=cross(v,hv)-r/rv; e=norm(e_vec); a=1/(2/rv-vv^2);
q2= H_hat(1)/(1+H_hat(3)+1e-300);
q1=-H_hat(2)/(1+H_hat(3)+1e-300);
f_hat=[1-q1^2+q2^2;2*q1*q2;-2*q1]/(1+q1^2+q2^2);
g_hat=[2*q1*q2;1+q1^2-q2^2;2*q2]/(1+q1^2+q2^2);
p1=dot(e_vec,g_hat); p2=dot(e_vec,f_hat);
cosF=dot(r/rv,f_hat); sinF=dot(r/rv,g_hat); F=atan2(sinF,cosF);
w_plus_O=atan2(p1,p2); nu=mod(F-w_plus_O,2*pi);
E_anom=2*atan2(sqrt(1-e)*sin(nu/2),sqrt(1+e)*cos(nu/2));
M_anom=E_anom-e*sin(E_anom); L=mod(M_anom+w_plus_O,2*pi);
y=[a;p1;p2;q1;q2;L];
end

function drdy = drdy_equinoctial(y0)
h=1e-6; drdy=zeros(3,6);
for j=1:6
    yp=y0; yp(j)=yp(j)+h; ym=y0; ym(j)=ym(j)-h;
    [rp,~]=equinoctial_to_cartesian_local(yp);
    [rm,~]=equinoctial_to_cartesian_local(ym);
    drdy(:,j)=(rp-rm)/(2*h);
end
end

function [r,v]=equinoctial_to_cartesian_local(x)
a=x(1);p1=x(2);p2=x(3);q1=x(4);q2=x(5);L=x(6);
e=sqrt(p1^2+p2^2); wO=atan2(p1,p2); M=mod(L-wO,2*pi);
E=M; for k=1:50; dE=(M-E+e*sin(E))/(1-e*cos(E)+1e-300); E=E+dE; if abs(dE)<1e-13;break;end;end
nu=2*atan2(sqrt(1+e)*sin(E/2),sqrt(1-e)*cos(E/2)); F=nu+wO;
s2=1+q1^2+q2^2; psl=a*(1-e^2); rm=psl/(1+p1*sin(F)+p2*cos(F));
fh=(1/s2)*[1-q1^2+q2^2;2*q1*q2;-2*q1]; gh=(1/s2)*[2*q1*q2;1+q1^2-q2^2;2*q2];
r=rm*(cos(F)*fh+sin(F)*gh); hm=sqrt(psl);
v=(-(1/hm)*(p1+sin(F)))*fh+((1/hm)*(p2+cos(F)))*gh;
end

function E=kepler_newton_local(M,e)
E=M; for k=1:50; dE=(M-E+e*sin(E))/(1-e*cos(E)+1e-300); E=E+dE; if abs(dE)<1e-13;break;end;end
end

function x=keplerian_flow_local(y0,t)
x=y0; x(6)=y0(6)+t;
end

function Phi=keplerian_STM_local(t,x)
Phi=eye(6); a=x(1); Phi(6,1)=-1.5*t*a^(-2.5);
end

function G=gve_equinoctial_local(x_true)
a=x_true(1);p1=x_true(2);p2=x_true(3);
q1=x_true(4);q2=x_true(5);F=x_true(6);
s2=1+q1^2+q2^2; e=sqrt(p1^2+p2^2); p=a*(1-e^2); h=sqrt(p); w=1+p1*sin(F)+p2*cos(F);
G_R=zeros(6,1);G_T=zeros(6,1);G_N=zeros(6,1);
G_R(1)=2*a^2*(p2*sin(F)-p1*cos(F))/h; G_R(2)=sin(F)/h; G_R(3)=-cos(F)/h;
G_R(6)=-(q1*sin(F)-q2*cos(F))/(h*w);
G_T(1)=2*a^2*w/h; G_T(2)=((w+1)*cos(F)+p1)/h; G_T(3)=-((w+1)*sin(F)-p2)/h;
G_T(6)=(q1*cos(F)+q2*sin(F))/(h*w);
G_N(4)=(s2*sin(F))/(2*h); G_N(5)=(s2*cos(F))/(2*h);
G_N(6)=-(s2*(q1*sin(F)-q2*cos(F)))/(2*h*w);
G=[G_R,G_T,G_N];
end

function F_true=mean_long_to_true_local(x)
p1=x(2);p2=x(3);L=x(6); e=sqrt(p1^2+p2^2); wO=atan2(p1,p2);
M=mod(L-wO,2*pi); E=kepler_newton_local(M,e);
nu=2*atan2(sqrt(1+e)*sin(E/2),sqrt(1-e)*cos(E/2)); F_true=nu+wO;
end

function D=D_matrix_local(tau,y0)
x=keplerian_flow_local(y0,-tau); Phi=keplerian_STM_local(tau,x);
Ft=mean_long_to_true_local(x); xt=x; xt(6)=Ft;
G=gve_equinoctial_local(xt); D=Phi*G;
end

function s=ternary(cond,a,b); if cond; s=a; else; s=b; end; end

% -----------------------------------------------------------------------------
function [t0_ref,theta_ref]=refine_red_shoot( ...
        c1_tgt,c2_tgt,wf_b1,wf_b2,t_vec,theta_k,kt0,ir0)
N_rays=numel(theta_k); Ntime_eff=numel(t_vec);
dt=t_vec(2)-t_vec(1); dtheta=theta_k(2)-theta_k(1);
clamp=@(x,lo,hi) max(lo,min(hi,x));
t0=t_vec(kt0); theta=theta_k(ir0);
for iter=1:20
    kt=clamp(interp_idx(t0,t_vec),2,Ntime_eff-1);
    ir=clamp(interp_idx(theta,theta_k),2,N_rays-1);
    c1_00=wf_b1(ir,kt);   c1_10=wf_b1(ir+1,kt);
    c1_01=wf_b1(ir,kt+1); c1_11=wf_b1(ir+1,kt+1);
    c2_00=wf_b2(ir,kt);   c2_10=wf_b2(ir+1,kt);
    c2_01=wf_b2(ir,kt+1); c2_11=wf_b2(ir+1,kt+1);
    ft=clamp((t0-t_vec(kt))/dt,0,1);
    fr=clamp((theta-theta_k(ir))/dtheta,0,1);
    c1_v=(1-fr)*(1-ft)*c1_00+fr*(1-ft)*c1_10+(1-fr)*ft*c1_01+fr*ft*c1_11;
    c2_v=(1-fr)*(1-ft)*c2_00+fr*(1-ft)*c2_10+(1-fr)*ft*c2_01+fr*ft*c2_11;
    res1=c1_v-c1_tgt; res2=c2_v-c2_tgt;
    if sqrt(res1^2+res2^2)<1e-8; break; end
    dc1_dt=((1-fr)*(c1_01-c1_00)+fr*(c1_11-c1_10))/dt;
    dc1_dth=((1-ft)*(c1_10-c1_00)+ft*(c1_11-c1_01))/dtheta;
    dc2_dt=((1-fr)*(c2_01-c2_00)+fr*(c2_11-c2_10))/dt;
    dc2_dth=((1-ft)*(c2_10-c2_00)+ft*(c2_11-c2_01))/dtheta;
    J=[dc1_dt,dc1_dth;dc2_dt,dc2_dth];
    detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);
    if abs(detJ)<1e-14; break; end
    step=J\[-res1;-res2];
    for s=1:2; if abs(step(s))>abs([dt;dtheta]*3); step(s)=sign(step(s))*abs([dt;dtheta(s)])*3; end; end
    t0=clamp(t0+step(1),t_vec(end),t_vec(1));
    theta=mod(theta+step(2)+pi,2*pi)-pi;
end
t0_ref=t0; theta_ref=theta;
end

function idx=interp_idx(val,vec)
[~,idx]=min(abs(vec-val)); idx=max(1,min(numel(vec)-1,idx));
if numel(vec)>1&&vec(2)<vec(1); if val>vec(idx); idx=max(1,idx-1); end
else; if numel(vec)>idx&&val>vec(idx+1); idx=min(numel(vec)-1,idx+1); end; end
end

% =============================================================================
%                              FIGURE
% =============================================================================
function generate_cam_figure(cam, sigma, epsilon)

wf         = cam.wf;
t_vec      = wf.t_vec;
wf_b1      = wf.wf_b1;
wf_b2      = wf.wf_b2;
t0_star    = cam.t0_star;
theta_star = cam.theta_star;
t_profile  = cam.t_profile;
u_RTN      = cam.u_RTN;
az         = cam.azimuth_deg;
el         = cam.elevation_deg;

[~, idx_t0] = min(abs(t_vec - t0_star));
N_rays  = size(wf_b1,1);
cmap    = parula(N_rays);
th_circ = linspace(0,2*pi,500);
step    = max(1,floor(N_rays/80));

% Normalised time in orbital periods for the profile axis
t_orb = t_profile / (2*pi);

figure('Name','Optimal CAM','Position',[60,60,1200,850]);

%% Panel 1: B-plane -----------------------------------------------------------
ax1 = subplot(2,2,[1,3]);   % left column, full height
fill(ax1,cos(th_circ),sin(th_circ),[0.92 0.92 0.92],...
    'EdgeColor',[0.5 0.5 0.5],'LineWidth',1.5,'FaceAlpha',0.55); hold(ax1,'on');

Ntime_eff = numel(t_vec);
for ir=1:step:N_rays
    plot(ax1, wf_b1(ir,1:min(idx_t0,Ntime_eff)), ...
              wf_b2(ir,1:min(idx_t0,Ntime_eff)), '-', ...
         'Color',[cmap(ir,:),0.40],'LineWidth',0.6);
end
wfc1=[wf_b1(:,min(idx_t0,Ntime_eff));wf_b1(1,min(idx_t0,Ntime_eff))];
wfc2=[wf_b2(:,min(idx_t0,Ntime_eff));wf_b2(1,min(idx_t0,Ntime_eff))];
plot(ax1,wfc1,wfc2,'r-','LineWidth',2,'DisplayName', ...
     sprintf('Wavefront at t_0^* = %.3f',t0_star));
scatter(ax1,cos(theta_star),sin(theta_star),100,'r','filled',...
    'MarkerEdgeColor','k','DisplayName','\hat{n}(\theta^*)');
scatter(ax1,cam.Dr0_b1,cam.Dr0_b2,120,'m','filled',...
    'MarkerEdgeColor','k','DisplayName','\Deltar_0/\sigma');
plot(ax1,[cam.Dr0_b1,cos(theta_star)],[cam.Dr0_b2,sin(theta_star)],...
    'm--','LineWidth',1.5,'DisplayName','Optimal extremal ray');
scatter(ax1,0,0,50,'k','filled','DisplayName','Origin');
axis(ax1,'equal'); grid(ax1,'on');
wmax=max(max(abs([wfc1;wfc2]))*1.15,1.4);
xlim(ax1,[-wmax,wmax]); ylim(ax1,[-wmax,wmax]);
xlabel(ax1,'\Deltar \cdot \bfb_1 / \sigma','FontSize',10);
ylabel(ax1,'\Deltar \cdot \bfb_2 / \sigma','FontSize',10);
title(ax1,sprintf('B-plane — t_0^* = %.4f norm. (%.3f orb. per.)', ...
    t0_star,cam.t0_star_orb),'FontSize',10);

%% Panel 2: RTN components of u*(t) ------------------------------------------
ax2 = subplot(2,2,2);
plot(ax2, t_orb, u_RTN(1,:), 'b-',  'LineWidth',1.5, 'DisplayName','u_R'); hold(ax2,'on');
plot(ax2, t_orb, u_RTN(2,:), 'r-',  'LineWidth',1.5, 'DisplayName','u_T');
plot(ax2, t_orb, u_RTN(3,:), 'g-',  'LineWidth',1.5, 'DisplayName','u_N');
yline(ax2, 0,'k:','LineWidth',0.8);
xline(ax2, 0,'k--','LineWidth',1.0,'DisplayName','TCA');
grid(ax2,'on');
xlabel(ax2,'Time [orbital periods] (0 = TCA)','FontSize',9);
ylabel(ax2,'u^* component (unit)','FontSize',9);
title(ax2,'Optimal thrust direction u^*(t) — RTN components','FontSize',10);
legend(ax2,'Location','best','FontSize',8); hold(ax2,'off');
xlim(ax2,[t_orb(1), 0]);

end
