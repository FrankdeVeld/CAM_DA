using LinearAlgebra, ForwardDiff, DifferentialEquations, StaticArrays, Plots, Interpolations, FFTW, Roots, SciMLBase, JLD2

default(fontfamily="Computer Modern")
default(grid=false)

include("Geometries.jl")
include("Conversions.jl")
include("LinearBackwardPropagation.jl")
include("MinTimeDynamics.jl")
include("DifferentialContinuation.jl") 
include("FWProp.jl")
include("BasicPlots.jl")

# Important parameters to change:
# MaxT: maximum integration time [s or revs], ensure it is above t0, but ideally not far. Line: 60
# Epsilon: maximum thrust acceleration magnitude [m/s²]. Line: 83
# Sigma: safe distance [m]. Line: 84
# Geometry at tca and Δvtca [equinoctial elements and m/s]. Line: 45
# Scaled initial condition on B-plane [b-plane coords]. Line: 129
# Note; evaluations assume perfect collision in first-order, above parameter is for setting actual trajectory

function run_simulation()
    cd("C:/Users/frank/Documents/GitHub/CAM_DA/Julia_code")   # Not sure if needed but change to own directory when running
    DimBool             = true    # Dimensionless units or not 
    ProblemBool         = true    # Integrating the problem 
    FWBool              = true    # Forward propagation for plots and validation. Optional (takes a while and is not written very nicely)
    ODB                 = 0;      # 1 True, 0 False. ODB means One-Dimensional Bool, e.g. one-dimensional control, prescribed three-axis attitude control. Leave false
    EclipseBool         = 0;      # 1 True, 0 False. Eclipse model on or off. Leave false
    # For this to be a real bool, make constants an Any[] object. For later

    Debug_mode          = true    # Explicit debug info print
    
    # Choose geometry
    # Scenario 1: free parameters, requires 6 orbital elements (Kepler), Delta V magnitude and direction
    # Scenario 2: toy model; primary in XY plane, secondary in XZ plane, Delta V practically in Z direction
    # Scenario 3: ISS type orbits for primary, secondary
    # Scenario 4: Lamberto validation data
    # Scenario 5: another toy model with different Delta V 
    # Scenario 6: another toy model. Looks similar to scn 2, not sure why this was made 
    # Scenario 7: idea is Iridium NEXT type orbit, with right now requires all variables
    # DimBool: Whether units are dimensionless or not
    Mu                     = 3.986004418 * 10^14;
    Scn                    = 1 
    
    #yp0, ys0, Δvtca        = GeometryPicker(Scn,Mu,aC = 6378.1e3 + 780e3, eC = 0.0001, iC= deg2rad(86.4), nuC=deg2rad(70.16), OmC = deg2rad(0.0), omC = 0.0, deltaVcMag = 1.0554e4, deltaVAngA=deg2rad(0.0), deltaVAngB=deg2rad(-45.0))
    yp0, ys0, Δvtca        = GeometryPicker(Scn,Mu,aC = 6378.1e3 + 705e3, eC = 0.0001, iC= deg2rad(98.2), nuC=deg2rad(107.3), OmC = deg2rad(70.5), omC = 0.0, deltaVcMag = 1.06094e4, deltaVAngA=deg2rad(0.0), deltaVAngB=deg2rad(-45.0))
    println("Sanity check; primary state at tCA in equinoctial elements: ", yp0) 
    println("Sanity check; secondary state at tCA in equinoctial elements: ", ys0) 
    println("Sanity check; primary state at tCA in ECI, Cartesian: ", equin2cart(yp0,Mu)) 
    println("Sanity check; secondary state at tCA in ECI, Cartesian: ", equin2cart(ys0,Mu)) 

    println("Sanity check; Δvtca in ECI, Cartesian: ", Δvtca) 
 
    Period                 = 2*pi*sqrt(yp0[1]^3/Mu)

    N                      = 100                              # Number of nodes temp
    Numθ                   = 2*N + 1                          # Number of nodes on circles (odd for FFT)
    NumT                   = 1000;                            # Number of time grid points to save during integration

    MaxT                   = -3*Period                        # Maximum integration time (seconds)
    MaxTAbs                = abs(MaxT)                        # Absolute value of MaxT (only for file saving purposes; minuses not appreciated there)

    # If DimBool true, remove dimensions of time, distance 
    if DimBool
        TimeDL             = 2*pi*sqrt(yp0[1]^3/Mu)
        DistanceDL         = yp0[1]
        MaxT               = MaxT/TimeDL;                    # Integration time in dimensionless units
        MaxTAbs            = MaxTAbs/TimeDL
        Mu                 = 4*pi^2                          # Dimensionless gravitational parameter
        # --- Create NEW SVectors using StaticArrays.setindex ---
        yp0_new_a = yp0[1] / DistanceDL
        yp0 = StaticArrays.setindex(yp0, yp0_new_a, 1)       # Returns a NEW SVector

        ys0_new_a = ys0[1] / DistanceDL
        ys0 = StaticArrays.setindex(ys0, ys0_new_a, 1)       # Returns a NEW SVector
        Δvtca              = Δvtca./(DistanceDL/TimeDL)      #/1000
    end

    Constants              = [Mu; yp0; Δvtca];               # Constants for integration
    day_of_year = 62
    NuES0                  = mod(2*pi * (day_of_year - 81) / 365.25, 2*pi)      # Initial true anomaly of Earth around Sun (0 = spring). Only relevant for yaw-steering

    Epsilon                = 1e-4                             # Thrust magnitude in m/s²
    Sigma                  = 10000                            # Safe distance for CAM in m
    RE                     = 6.378e6 

    if DimBool
        Epsilon            = Epsilon/(DistanceDL/TimeDL^2);  # Epsilon in dimensionless units
        Sigma              = Sigma/DistanceDL                # Sigma in dimensionless units 
        RE                 = RE/DistanceDL
    end

    TuningP = 0.0         # For control smoothing, not fully implemented currently
    
    # --- Create a NamedTuple for parameters ---
    sim_params = (
        Mu = Mu,
        yp0 = SVector{6, Float64}(yp0), # Convert to SVector if using StaticArrays
        ys0 = SVector{6, Float64}(ys0),
        Δvtca = SVector{3, Float64}(Δvtca),
        Epsilon = Epsilon,
        Sigma = Sigma,
        NuES0 = NuES0,
        RE = RE,
        ODB = ODB, # Keep as Int
        EclipseBool = EclipseBool, # Keep as Int
        LuxBool = 0,
        TuningP = TuningP,
        Re_const = RE,
        forbidden_intervals = [(0.0, 0.01)] # Forbidden intervals in positive t range mean they are not used
        # Add other parameters needed inside IntδΔr! (like forbidden_intervals if passed this way)
    )

    if ProblemBool
        tProp, δΔrtCA, uOptRSW, B1Ref, B3U, fft_coeffs_b1_t, fft_coeffs_b1_dt_t, fft_coeffs_b2_t, fft_coeffs_b2_dt_t = IntegrateCircleδΔr(Numθ, NumT, MaxT, sim_params) 
        common_time_grid_for_dc = vec(tProp[1,:])   

        println("Sanity check: B1 unit vector: ", B1Ref)
        println("Sanity check: B2 unit vector: ", cross(B3U, B1Ref) )
        println("Sanity check: B3 unit vector: ", B3U)
        ΔR_ESS_F                           = δΔrScaling(Numθ, NumT, δΔrtCA, Epsilon, Sigma, B1Ref, B3U)
    else
        println("ProblemBool is false. Skipping.")
    end

    println("Starting Differential Continuation (FFT-based)...")

    # InitCon_physical is the target in the B-plane, scaled by Sigma (so it's dimensionless, target on unit circle if no perturbation)
    InitCon_B_scaled = SVector{3,Float64}(0.0, -0.1, 0.0) # Using the one defined above
    
    # The time grid for FFT coefficient series is tProp[1,:]
    # Ensure Numθ is odd for get_wavenumbers, or adapt get_wavenumbers
    if !isodd(Numθ)
        @warn "Numθ ($Numθ) is even. DC_FFT works best with odd Numθ. Consider adjusting."
        # If you must use even Numθ, the get_wavenumbers and FFT indexing needs careful review.
        # For now, it might error out in solve_dc_fft.
    end

    # Initialize variables to store the best result found
    best_t0_from_dc = -Inf # Initialize with a value that any valid t0 (negative) will be less than
    best_theta_f_from_dc = NaN
    best_retcode_from_dc = :NotRun # Or some other initial state

    # Solve for both branches
    for branch_idx_dc in 1:2
        current_t0_dc, current_theta_f_dc, current_retcode_dc = solve_dc_fft(
            InitCon_B_scaled,
            common_time_grid_for_dc,
            fft_coeffs_b1_t,
            fft_coeffs_b1_dt_t,
            fft_coeffs_b2_t,
            fft_coeffs_b2_dt_t,
            Numθ, # This is num_theta_dc
            initial_theta_branch=branch_idx_dc,
            orientation=:horizontal, # or :vertical depending on your problem
            #orientation=:vertical,
        )
        if Debug_mode
            println("DC Branch $branch_idx_dc Result: t0=$(current_t0_dc), theta_f=$(rad2deg(current_theta_f_dc)) deg, retcode=$(current_retcode_dc)")
        end
        # Check if this branch yielded a valid, better solution
        if current_retcode_dc == SciMLBase.ReturnCode.Terminated && !isnan(current_t0_dc)
            # We are looking for the most negative t0 (latest time to start)
            if current_t0_dc > best_t0_from_dc # Note: t0 is negative, so '>' means "less negative" or "later start"
                best_t0_from_dc = current_t0_dc
                best_theta_f_from_dc = current_theta_f_dc
                best_retcode_from_dc = current_retcode_dc
                if Debug_mode
                    println("Branch $branch_idx_dc provided a new best solution.")
                end
            end
        end
    end
    
    # Assign the final selected values
    t0_from_dc = best_t0_from_dc
    theta_f_from_dc = best_theta_f_from_dc 
    retcode_from_dc = best_retcode_from_dc

    # Check if any valid solution was found
    if isinf(t0_from_dc) # Means no branch terminated successfully with a valid t0
        println("Warning: Differential Continuation did not find a valid terminating solution for any branch.")
        t0_from_dc = NaN # Set to NaN if no solution found
        # theta_f_from_dc will remain NaN
        retcode_from_dc = :Failure_NoValidBranchSolution # Custom retcode
    end
    
    println("Selected DC Solution: t0=$(t0_from_dc), theta_f=$(rad2deg(best_theta_f_from_dc)) deg, retcode=$(retcode_from_dc)")

    ### Other code which is not needed has been removed
    if FWBool
        ### Propagating in the linearised system contained on the B-plane
        
        t0                                = t0_from_dc
        Theta                             = theta_f_from_dc
        tf                                = 0
        NTS                               = 100 # Number of time steps for forward propagation
        FW_t, FW_Δr_Traj, FW_u_RSW        = FP_L_Δr(Mu, Epsilon, Sigma, t0, tf, Theta, Δvtca, yp0, ys0, NTS,  NuES0, ODB, EclipseBool, 0.0, RE) 
        
        filename               = "test"
        filenamet              = filename * "t.jld2"
        filenameTraj           = filename * "Dr.jld2"
        filenameU              = filename * "u.jld2"

        @save filenamet       FW_t
        @save filenameTraj    FW_Δr_Traj
        @save filenameU       FW_u_RSW

        BPlaneTrajectoryPlot(FW_Δr_Traj)
        DeltaRMagPlot(FW_Δr_Traj, FW_t)
        uRSWPlot(FW_t, FW_u_RSW)

        ### Propagating xp and xs from t0 and tf to check true tca change and true separation distance
        ## Obtaining xp(t0) and xs(t0), with xp(t0) from the initial condition on the B-plane 
        rs0E                 = equin2cart(ys0,Mu)[1:3]
        r0E                  = B2EU(InitCon_B_scaled, B1Ref, B3U)*Sigma + rs0E
        y0V                  = equin2cart(yp0,Mu)[4:6]
        xtcaE                = SVector{6}([r0E;y0V])
        yp0Int               = cart2equin(xtcaE,Mu)

        ## Propagating xp (no u) and xs from t0 to tf=0 + a little bit 
        ## Propagating xp (with u) and xs from t0 to tf=0 + a little bit 
        NTS = Int(1e6) #Fine grid
        Buffertime = 0.0001
        t_NL_C_Xp, Traj_NL_C_Xp, t_NL_NC_Xp, Traj_NL_NC_Xp, t_NL_Xs, Traj_NL_Xs   = FP_NL_x(Mu, Epsilon, t0_from_dc, tf+Buffertime, tf+Buffertime, theta_f_from_dc, Δvtca, yp0, yp0Int, ys0, NTS, NuES0, ODB, TuningP)
        
        ## Finding the nominal separation distance and tca 
        RelDisNoControl                      = Traj_NL_NC_Xp - Traj_NL_Xs
        tCA_NC, mindis_NC                    = PlotRelDis(t_NL_NC_Xp, RelDisNoControl, false)
        println("Nominal tCA at: ", tCA_NC*TimeDL, " s")
        println("Nominal miss distance: ", mindis_NC*DistanceDL, " m")

        ## Finding the adjusted separation distance and tca 
        RelDisWithControl                    = Traj_NL_C_Xp - Traj_NL_Xs
        tCA_C, mindis_C                      = PlotRelDis(t_NL_C_Xp, RelDisWithControl, false)
        println("Adjusted tCA at: ", tCA_C*TimeDL, " s")
        println("Adjusted miss distance: ", mindis_C*DistanceDL, " m")

        ## Note: tf=0 works in first-order approximation, but the method tracks states at the 'real' adjusted tca
        ## Even if tca evolution itself is not extracted from the linearised solution (could be extracted from full non-linear propagation)
        ## The code for the non-linear propagation has been removed from this code version for clarity
    end
end 

run_simulation()
println("Script finished.")