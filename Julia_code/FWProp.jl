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
    C_Δr_NL_Sol               = solve(C_Δr_NL_Prob, Tsit5(), abstol = STA, reltol = STR, saveat = range(tf, t0, length=NTS));

    tInt                      = C_Δr_NL_Sol.t;           # Unpacking solution (t0)
    δΔr_E                     = C_Δr_NL_Sol[1:3,:];      # Unpacking solution (deltaz)  
    
    Δr_E                      = zeros(T, NTS, 3)
    Δr_B                      = zeros(T, NTS, 3)
    N_B                       = zeros(T, NTS, 3)
    ΔR_ESS_F                  = zeros(T, NTS, 3)
    u_RSW                     = zeros(T, NTS, 3)
    InnerProduct              = zeros(T, NTS)

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
    
    x       = SVector{7}(X[1:7]);
    
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
        @SVector [0.0, 0.0, 0.0, uOptECI_int[1], uOptECI_int[2], uOptECI_int[3],norm(uOptECI_int) * epsilon]
    else
        @SVector zeros(7) # A clear way to get a zero SVector
    end
    # Calculate the derivative and return it
    dx = @SVector[x[4], x[5], x[6], -Mu * x[1] / norm(x[1:3])^3, -Mu * x[2] / norm(x[1:3])^3, -Mu * x[3] / norm(x[1:3])^3, norm(uOptRSW) * epsilon]

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
        State            = Sol[1:6,:];
        #Xf               = Sol[1:6,end]
    end    
    return t,State#,Xf
end

function FP_NL_x(Mu, Epsilon, t0, tf, tstop, Theta, Δvtca, yp0, yp0Int, ys0, NuES0, ODB,TuningP, JuliaBool, TimeDL)

    # IntTime: positive (abs(t0))
    # Xp0, Xs0: states of primary, secondary at t0 
    # If methodology correct: propagation for time t0 will ensure delta R(t_f) = sigma

    tspan_fw                  = (t0, tf);                                # Forward integration
    NECI                      = Theta2NThetaE(Δvtca, yp0, Mu, Theta)
    OTE                       = ∂r∂yNoT(Mu, yp0)

    C                         = [Mu; yp0; Δvtca; NECI; OTE[1,1:6]; OTE[2,1:6]; OTE[3,1:6]; NuES0; ODB; TuningP; Epsilon; tstop];                     # Constants


    if JuliaBool
        # 1. Load the pre-calculated data
        filename               = "DAComp1"
        filenamet              = filename * "t.jld2"
        filenameU              = filename * "u.jld2"

        @load filenamet       FW_t
        @load filenameU       FW_u_RSW

        t_grid = reverse(FW_t)        # Must be increasing for interpolation
        u_grid = FW_u_RSW[end:-1:1, :] # Reverses the order of the rows
    else 
        # 1. Load the pre-calculated data
        filename               = "Matlab1"
        filenamet              = filename * "_t.csv"
        filenameU              = filename * "_u.csv"

        # Read the time data and convert it to a 1D vector
        FW_t = vec(readdlm(filenamet, ',', Float64))
        
        # Read the U data which will automatically parse as an N x 3 Matrix
        FW_u_RSW = readdlm(filenameU, ',', Float64) 

        t_grid = FW_t
        u_grid = FW_u_RSW
    end
    
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
        control_itp = control_function,
    )

    Xp0                       = SVector{6}(equin2cart(Flow(yp0Int, t0, Mu),Mu))
    X0                        = Xp0 

    X0C                       = SVector{7}(X0[1],X0[2],X0[3],X0[4],X0[5],X0[6],0)

    STA                       = 1e-10                            # Absolute solver tolerance
    STR                       = 1e-10                            # Relative solver tolerance

    dt_save = abs(tf-t0)/1e6

    tgrid = t0:dt_save:tf
    
    C_Xp_NL_Prob              = ODEProblem{false}(EOM_C_NL, X0C, tspan_fw, ode_params)
    C_Xp_NL_Sol               = solve(
        C_Xp_NL_Prob,
        Tsit5(),
        abstol = STA,
        reltol = STR,
        saveat = tgrid
    
        )   
    NC_Xp_NL_Prob             = ODEProblem{false}(EOM_NC_NL, X0, tspan_fw, ode_params)
    NC_Xp_NL_Sol              = solve(
        NC_Xp_NL_Prob,
        Tsit5(),
        abstol = STA,
        reltol = STR,
        saveat = tgrid
    )
    Xs0                       = equin2cart(Flow(ys0, t0, Mu),Mu)


    NC_Xs_NL_Prob             = ODEProblem{false}(EOM_NC_NL, Xs0, tspan_fw, C)
    NC_Xs_NL_Sol              = solve(
        NC_Xs_NL_Prob,
        Tsit5(),
        abstol = STA,
        reltol = STR,
        saveat = tgrid
    )
    # Create a continuous function for the distance at ANY time `t`
    distance_NC_xp_NC_Xs(t)              = norm(NC_Xp_NL_Sol(t)[1:3] - NC_Xs_NL_Sol(t)[1:3])
    distance_C_xp_NC_Xs(t)               = norm( C_Xp_NL_Sol(t)[1:3] - NC_Xs_NL_Sol(t)[1:3])

    # Use an optimization package (like Optim.jl) to find the absolute true minimum
    trange = abs(tf-t0)
    res_NC                               = optimize(distance_NC_xp_NC_Xs, tf - trange/100, tf)
    res_C                                = optimize(distance_C_xp_NC_Xs, tf - trange/100, tf)
    true_tCA_NC                          = res_NC.minimizer
    true_mindis_NC                       = res_NC.minimum
    true_tCA_C                           = res_C.minimizer
    true_mindis_C                        = res_C.minimum
    total_Deltav                         = C_Xp_NL_Sol(true_tCA_C)[7]

    # C_Xp_NL_t, C_Xp_NL_State    = SolUnpack(C_Xp_NL_Sol,false);
    # NC_Xp_NL_t, NC_Xp_NL_State  = SolUnpack(NC_Xp_NL_Sol,false);
    # NC_Xs_NL_t, NC_Xs_NL_State  = SolUnpack(NC_Xs_NL_Sol,false);

    # # Finding the nominal separation distance and tca 
    # RelDisNoControl                      = NC_Xp_NL_State[1:3,:] - NC_Xs_NL_State[1:3,:]
    # RelVelNoControl                      = NC_Xp_NL_State[4:6,:] - NC_Xs_NL_State[4:6,:]
    # tCA_NC, mindis_NC, mininner_NC       = PlotRelDis(NC_Xp_NL_t, RelDisNoControl,RelVelNoControl, TimeDL, true, JuliaBool)
    
    
    # # Finding the adjusted separation distance and tca 
    # RelDisWithControl                    = C_Xp_NL_State[1:3,:] - NC_Xs_NL_State[1:3,:]
    # RelVelWithControl                    = C_Xp_NL_State[4:6,:] - NC_Xs_NL_State[4:6,:]
    # tCA_C, mindis_C, mininner_C,         = PlotRelDis(C_Xp_NL_t, RelDisWithControl, RelVelWithControl, TimeDL ,true, JuliaBool)

    return true_tCA_NC, true_mindis_NC, true_tCA_C, true_mindis_C,total_Deltav     #C_Xp_NL_t, C_Xp_NL_State, NC_Xp_NL_t, NC_Xp_NL_State, NC_Xs_NL_t, NC_Xs_NL_State
end