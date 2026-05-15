######### LIN CON, PROP Δr, LINEAR PROP:  Δr trajectory (B Plane) #########
include("MinTimeDynamics.jl")
include("Conversions.jl")

### Propagating Δr
function FP_L_Δr(Mu, Epsilon, Sigma, t0, tf, Theta, Δvtca, yp0, ys0, NTS,  NuES0, ODB, EclipseBool, TuningP, RE)

    t0                        = t0
    tspan                     = (tf, t0);                                # Forward integration
    B1Ref, B3U                = B1B3(Δvtca, yp0, Mu)
    T                         = eltype(Δvtca)
    zero_T                    = zero(T)
    cosT_T                    = cos(Theta); sinT_T = sin(Theta)
    NBPlane_T                 = SA[cosT_T, sinT_T, zero_T]
    NECI                      = B2EU(NBPlane_T,B1Ref, B3U)
    OTE                       = ∂r∂yNoT(Mu, yp0)

    # --- Create a NamedTuple for parameters ---
    sim_params = (
        Mu = Mu,
        yp0 = SVector{6, Float64}(yp0), # Convert to SVector if using StaticArrays
        ys0 = SVector{6, Float64}(ys0),
        Δvtca = SVector{3, Float64}(Δvtca),
        NTheta = NECI,
        OTE = OTE,
        NuES0 = NuES0,
        ODB = ODB, # Keep as Int
        EclipseBool = EclipseBool, # Keep as Int
        LuxBool = 0,
        TuningP = TuningP,
        Re_const = RE,
        forbidden_intervals = [(0.0, 0.01)]
        # Add other parameters needed inside IntδΔr! (like forbidden_intervals if passed this way)
    )

#    C                         = [Mu; yp0; Δvtca; NECI; OTE[1,1:6]; OTE[2,1:6]; OTE[3,1:6]; NuES0; ODB; TuningP];                     # Constants

    δΔr0                      = zeros(T, 3)
    STA                       = 1e-10                            # Absolute solver tolerance
    STR                       = 1e-10                            # Relative solver tolerance

    C_Δr_NL_Prob              = ODEProblem((dδΔr, δΔr, C, t) -> IntδΔr!(dδΔr, δΔr, C, t), δΔr0, tspan, sim_params);
    C_Δr_NL_Sol               = solve(C_Δr_NL_Prob, Tsit5(), abstol = STA, reltol = STR, saveat = range(tf, t0, length=100));

    tInt                      = C_Δr_NL_Sol.t;           # Unpacking solution (t0)
    δΔr_E                     = C_Δr_NL_Sol[1:3,:];      # Unpacking solution (deltaz)  
    
    Δr_E                      = zeros(T, NTS, 3)# #zeros(NTS,3);
    Δr_B                      = zeros(T, NTS, 3)# #zeros(NTS,3);
    N_B                       = zeros(T, NTS, 3)# #zeros(NTS,3)
    ΔR_ESS_F                  = zeros(T, NTS, 3)# #zeros(NTS,3); 
    u_RSW                     = zeros(T, NTS, 3)# #zeros(NTS,3)
    InnerProduct              = zeros(T, NTS)#zeros(NTS);

    scaling_factor::T = Epsilon / Sigma

    for i in 1:length(δΔr_E[1,:])
        # Extract δΔr(t) as SVector{T} for calculations
        # Use @view for efficiency, then convert to SVector{T}
        δΔr_eci = SVector{3, T}(@view δΔr_E[:, i])

        Δr_S_E = scaling_factor * δΔr_eci # SVector{3, T}
        # E2BU takes SVector{3, T}, returns SVector{3, T}
        Δr_S_B = E2BU(Δr_S_E, B1Ref, B3U)
        N_B[i,:]              = SA[cosT_T, sinT_T, zero_T]

        # Assign result (SVector{3, T}) back to standard Array{T, 3} slice
        ΔR_ESS_F[i, :] = N_B[i,:]  + Δr_S_B
        u_RSW[i,:], InnerProduct[i]            = GetCon(sim_params,tInt[i])
    end

    return tInt, ΔR_ESS_F, u_RSW
end


### Propagating x
function EOM_C_NL(X, p,t) 
    Mu      = p.Mu
    epsilon = p.epsilon
    tstop   = p.tstop
    
    x       = SVector{6}(X[1:6]);
    
    if t < tstop
        # Pass the whole params object to GetCon if it needs more info
        #uOptRSW, _ = GetCon(p, t)
        uOptRSW = p.control_itp(t)
    else 
        uOptRSW = SA[0.0, 0.0, 0.0]
    end

    r = SVector{3}(x[1:3])
    v = SVector{3}(x[4:6])
    uDyn = if norm(uOptRSW) != 0
        uOptECI_int = R2EU(uOptRSW, r, v)
        # Construct an SVector directly
        @SVector [0.0, 0.0, 0.0, uOptECI_int[1], uOptECI_int[2], uOptECI_int[3]]
    else
        @SVector zeros(6) # A clear way to get a zero SVector
    end
    # Calculate the derivative and return it
    dx = @SVector[x[4], x[5], x[6], -Mu * x[1] / norm(x[1:3])^3, -Mu * x[2] / norm(x[1:3])^3, -Mu * x[3] / norm(x[1:3])^3]

    return dx + uDyn * epsilon
end

function EOM_NC_NL(X, C,t) # Nonlinear dynamics
    Mu            = C[1]; 
    x             = X[1:6];
    # Calculate and return the derivative as an SVector
    return @SVector [x[4], x[5], x[6], -Mu * x[1] / norm(x[1:3])^3, -Mu * x[2] / norm(x[1:3])^3, -Mu * x[3] / norm(x[1:3])^3]
end

function SolUnpack(Sol,bool)
    if bool
        t                = Sol.t;
        Pos              = Sol[1:3,:] + Sol[7:9,:];
        #Xf               = Sol[1:6,end]
    else
        t                = Sol.t;
        Pos              = Sol[1:3,:];
        #Xf               = Sol[1:6,end]
    end    
    return t,Pos#,Xf
end

function FP_NL_x(Mu, Epsilon, t0, tf, tstop, Theta, Δvtca, yp0, yp0Int, ys0, NTS, NuES0, ODB,TuningP)
    # IntTime: positive (abs(t0))
    # Xp0, Xs0: states of primary, secondary at t0 
    # If methodology correct: propagation for time t0 will ensure delta R(t_f) = sigma

    tspan_fw                  = (t0, tf);                                # Forward integration
    NECI                      = Theta2NThetaE(Δvtca, yp0, Mu, Theta)
    OTE                       = ∂r∂yNoT(Mu, yp0)

    C                         = [Mu; yp0; Δvtca; NECI; OTE[1,1:6]; OTE[2,1:6]; OTE[3,1:6]; NuES0; ODB; TuningP; Epsilon; tstop];                     # Constants



    # 1. Load the pre-calculated data
    filename               = "test"
    filenamet              = filename * "t.jld2"
    filenameU              = filename * "u.jld2"

    @load filenamet       FW_t
    @load filenameU       FW_u_RSW

    t_grid = reverse(FW_t)        # Must be increasing for interpolation
    u_grid = FW_u_RSW[end:-1:1, :] # Reverses the order of the rows

    # 2. Prepare the data for Interpolations.jl
    u_data_svectors = [SVector{3}(row) for row in eachrow(u_grid)]

    # 3. Create the interpolation object
    itp = LinearInterpolation(t_grid, u_data_svectors)

    # 4. Get the minimum and maximum time of the grid
    t_min = first(t_grid)
    t_max = last(t_grid)

    # 5. Create our custom "Flat" extrapolation wrapper using clamp
    # If t goes outside [t_min, t_max], clamp safely forces it to the nearest edge
    control_function = t -> itp(clamp(t, t_min, t_max))

        ode_params = (
        Mu = C[1],
        yp0 = SVector{6}(C[2:7]),
        Δvtca = SVector{3}(C[8:10]),
        NTheta = SVector{3}(C[11:13]),
        OTE = SMatrix{3,6}(C[14:31]), 
        NuES0 = C[32],
        ODB = 0,
        EclipseBool = 0,
        LuxBool = 0,
        TuningP= 0.0,
        Re_const = 6.378e6,
        forbidden_intervals = [(0.0, 0.01)],
        epsilon = C[35],
        tstop = C[36],
        control_itp = control_function
    )

    Xp0                       = SVector{6}(equin2cart(Flow(yp0Int, t0, Mu),Mu))
    X0                        = Xp0

    STA                       = 1e-10                            # Absolute solver tolerance
    STR                       = 1e-10                            # Relative solver tolerance

    C_Xp_NL_Prob              = ODEProblem{false}(EOM_C_NL, X0, tspan_fw, ode_params);
    save_points               = t0 .+ (tf - t0) .* range(0, 1, length=NTS).^3 #Ensure it mostly saves at the end
    C_Xp_NL_Sol               = solve(C_Xp_NL_Prob, Tsit5(), abstol = STA, reltol = STR, saveat = save_points);

    NC_Xp_NL_Prob             = ODEProblem{false}(EOM_NC_NL, X0, tspan_fw, ode_params)
    # save_points can be reused
    NC_Xp_NL_Sol              = solve(NC_Xp_NL_Prob, Tsit5(), abstol = STA, reltol = STR, saveat = save_points)

    Xs0                       = equin2cart(Flow(ys0, t0, Mu),Mu)

    NC_Xs_NL_Prob             = ODEProblem{false}(EOM_NC_NL, Xs0, tspan_fw, C)
    # save_points can be reused
    NC_Xs_NL_Sol              = solve(NC_Xs_NL_Prob, Tsit5(), abstol = STA, reltol = STR, saveat = save_points)

    C_Xp_NL_t, C_Xp_NL_Pos    = SolUnpack(C_Xp_NL_Sol,false);
    NC_Xp_NL_t, NC_Xp_NL_Pos  = SolUnpack(NC_Xp_NL_Sol,false);
    NC_Xs_NL_t, NC_Xs_NL_Pos  = SolUnpack(NC_Xs_NL_Sol,false);

    return C_Xp_NL_t, C_Xp_NL_Pos, NC_Xp_NL_t, NC_Xp_NL_Pos, NC_Xs_NL_t, NC_Xs_NL_Pos
end