using LinearAlgebra, ForwardDiff, DifferentialEquations, StaticArrays, Plots, Interpolations, FFTW, Roots, SciMLBase, JLD2, Optim, DelimitedFiles

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
# MaxT: maximum integration time [s or revs], ensure it is above t0
# Note; evaluations assume perfect collision in first-order

Length_Array = 2170
Skip = 0
t0_Array = zeros(Length_Array - Skip)
Dv_Array = zeros(Length_Array - Skip)
cd("C:/Users/frank/Documents/GitHub/CAM_DA/Julia_code")   # Not sure if needed but change to own directory when running
xpxs = readdlm("xpxsdata.csv", ',', Float64)

#for i=1:Length_Array
for i=6
    Epsilon                = 0.3/800                          # Thrust magnitude in m/s²
    Sigma                  = 2000.0                            # Safe distance for CAM in m
    Factor                 = 1.225;
    for j=1:3
        xpxs[i,j+3] = xpxs[i,j+3] * Factor;  # Higher ecc
    end
    xptcaTrue              = SVector{6}(xpxs[i,1:6].*1000)                 # m, m/s 
    xstcaECI               = SVector{6}(xpxs[i,7:12].*1000)                 # m, m/s
    xptcaECI               = SVector{6}(xstcaECI[1],xstcaECI[2],xstcaECI[3],xptcaTrue[4],xptcaTrue[5],xptcaTrue[6])
    # for j=1:3
    #     xptcaTrue[j+3] = xptcaTrue[j+3] * 1.2248;
    #     xptcaECI[j+3] = xptcaECI[j+3] * 1.2248;
    # end
        
    t0_Array[i], Dv_Array[i]       = run_simulation(xptcaECI,xstcaECI,xptcaTrue, Epsilon, Sigma)
    @show i
end
#writedlm("AllResultst0.csv", t0_Array, ',')
#writedlm("AllResultsDv.csv", Dv_Array, ',')

println("Script finished.")

function run_simulation(xptcaECI,xstcaECI,xptcaTrue,Epsilon, Sigma)
    
    DimBool             = true   # Dimensionless units or not 
    ProblemBool         = true    # Integrating the problem 
    FWBool              = true    # Forward propagation for plots and validation. Optional (takes a while and is not written very nicely)
    ODB                 = 0;      # 1 True, 0 False. ODB means One-Dimensional Bool, e.g. one-dimensional control, prescribed three-axis attitude control. Leave false
    EclipseBool         = 0;      # 1 True, 0 False. Eclipse model on or off. Leave false
    # For this to be a real bool, make constants an Any[] object. For later

    Debug_mode          = true   # Explicit debug info print
    
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
    #Scn                    = 1 
    
    #yp0, ys0, Δvtca        = GeometryPicker(Scn,Mu,aC = 6378.1e3 + 780e3, eC = 0.0001, iC= deg2rad(86.4), nuC=deg2rad(70.16), OmC = deg2rad(0.0), omC = 0.0, deltaVcMag = 1.0554e4, deltaVAngA=deg2rad(0.0), deltaVAngB=deg2rad(-45.0))
    #yp0, ys0, Δvtca        = GeometryPicker(Scn,Mu,aC = 6378.1e3 + 705e3, eC = 0.0001, iC= deg2rad(98.2), nuC=deg2rad(107.3), OmC = deg2rad(70.5), omC = 0.0, deltaVcMag = 1.06094e4, deltaVAngA=deg2rad(0.0), deltaVAngB=deg2rad(-45.0))
    
    ### DA COMP 1
    
    # Assume perfect collision for evaluations
    #xptcaECI = SVector{6}([2.333465506263321e+03, -1.103671212478364e+06, 7.105914958099038e+06, -7.442862828717730e3, -6.137347436526596e-01,  0.003951361392933e3])
    #xstcaECI = SVector{6}([2.333465506263321e+03, -1.103671212478364e+06, 7.105914958099038e+06,  7.353740487126315e3, -1.142814049765362e+02, -0.198247225911377e3])
    ΔrtcaECI = SVector{3}([xptcaECI[1] - xstcaECI[1],  xptcaECI[2] - xstcaECI[2], xptcaECI[3] - xstcaECI[3]])
    ΔvtcaECI = SVector{3}([xptcaECI[4] - xstcaECI[4],  xptcaECI[5] - xstcaECI[5], xptcaECI[6] - xstcaECI[6]])
    
    yp0 = cart2equin(xptcaECI, Mu)
    ys0 = cart2equin(xstcaECI, Mu)
    if Debug_mode
        println("Sanity check; (idealised) primary state at tCA in equinoctial elements: ", yp0) 
        println("Sanity check; secondary state at tCA in equinoctial elements: ", ys0) 
        println("Sanity check; (idealised) primary state at tCA in ECI, Cartesian: ", equin2cart(yp0,Mu)) 
        println("Sanity check; secondary state at tCA in ECI, Cartesian: ", equin2cart(ys0,Mu)) 

        println("Sanity check; Δvtca in ECI, Cartesian: ", ΔvtcaECI) 
    end
    Period                 = 2*pi*sqrt(yp0[1]^3/Mu)

    N                      = 50                              # Number of nodes temp
    Numθ                   = 2*N + 1                          # Number of nodes on circles (odd for FFT)
    NumT                   = 500;                            # Number of time grid points to save during integration

    MaxT                   = -1*Period                        # Maximum integration time (seconds)
    MaxTAbs                = abs(MaxT)                        # Absolute value of MaxT (only for file saving purposes; minuses not appreciated there)
    @show DimBool
    # If DimBool true, remove dimensions of time, distance 
    if DimBool
        TimeDL             = 2*pi*sqrt(yp0[1]^3/Mu)
        println("Orbital period: ", TimeDL)
        DistanceDL         = yp0[1]
        MaxT               = MaxT/TimeDL;                    # Integration time in dimensionless units
        MaxTAbs            = MaxTAbs/TimeDL
        Mu                 = 4*pi^2                          # Dimensionless gravitational parameter
        # --- Create NEW SVectors using StaticArrays.setindex ---
        yp0_new_a = yp0[1] / DistanceDL
        yp0 = StaticArrays.setindex(yp0, yp0_new_a, 1)       # Returns a NEW SVector

        ys0_new_a = ys0[1] / DistanceDL
        ys0 = StaticArrays.setindex(ys0, ys0_new_a, 1)       # Returns a NEW SVector
        ΔvtcaECI              = ΔvtcaECI./(DistanceDL/TimeDL)      #/1000
        if Debug_mode
            @show TimeDL 
            @show DistanceDL
        end
    end

    Constants              = [Mu; yp0; ΔvtcaECI];               # Constants for integration
    day_of_year = 62
    NuES0                  = mod(2*pi * (day_of_year - 81) / 365.25, 2*pi)      # Initial true anomaly of Earth around Sun (0 = spring). Only relevant for yaw-steering

    ### DA COMP 1
    #Epsilon                = 0.3/800                          # Thrust magnitude in m/s²
    #Sigma                  = 2000.0                            # Safe distance for CAM in m
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
        Δvtca = SVector{3, Float64}(ΔvtcaECI),
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
        if Debug_mode
            println("Sanity check: B1 unit vector: ", B1Ref)
            println("Sanity check: B2 unit vector: ", cross(B3U, B1Ref) )
            println("Sanity check: B3 unit vector: ", B3U)
        end
        ΔR_ESS_F                           = δΔrScaling(Numθ, NumT, δΔrtCA, Epsilon, Sigma, B1Ref, B3U)
    else
        println("ProblemBool is false. Skipping.")
    end

    if Debug_mode
        println("Starting Differential Continuation (FFT-based)...")
    end

    # InitCon_physical is the target in the B-plane, scaled by Sigma (so it's dimensionless, target on unit circle if no perturbation)
    ### DA COMP 1
    # The true orbit from which we will extract the position on the B-plane
    #xptcaTrue = SVector{6}([2.330521851751368e+03, -1.103704510502015e+06, 7.105887642997178e+06, -7.442862828717730e3, -6.137347436526596e-01,  0.003951361392933e3])
    if DimBool
        InitCon_B_scaled = E2BU( SVector{3}(xptcaTrue[1] - xstcaECI[1],xptcaTrue[2]- xstcaECI[2],xptcaTrue[3]- xstcaECI[3]), B1Ref, B3U)./DistanceDL./Sigma
        if Debug_mode
            println("Sanity check; true primary state at tCA on B-plane: ", InitCon_B_scaled, ", initial miss distance: ", norm(InitCon_B_scaled)*Sigma*DistanceDL , " m") 
        end
    else 
        InitCon_B_scaled = E2BU( SVector{3}(xptcaTrue[1] - xstcaECI[1],xptcaTrue[2]- xstcaECI[2],xptcaTrue[3]- xstcaECI[3]), B1Ref, B3U)./Sigma
        if Debug_mode
            println("Sanity check; true primary state at tCA on B-plane: ", InitCon_B_scaled, ", initial miss distance: ", norm(InitCon_B_scaled)*Sigma , " m") 
        end
    end
    #InitCon_B_scaled = SVector{3,Float64}(0.0, 0.0, 0.0)

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
    if Debug_mode
        println("Selected DC Solution: t0=$(t0_from_dc), theta_f=$(rad2deg(best_theta_f_from_dc)) deg, retcode=$(retcode_from_dc)")
    end
    if Debug_mode
        if DimBool
            println("Manoeuvre starts " , abs(t0_from_dc * TimeDL) , " s before nominal tca")
        else
            println("Manoeuvre starts " , abs(t0_from_dc) , " s before nominal tca")
        end
    end
    ### Other code which is not needed has been removed
    if FWBool
        ### Propagating in the linearised system contained on the B-plane
        
        t0                                = t0_from_dc
        Theta                             = theta_f_from_dc
        tf                                = 0
        NTS_Δr                            = 500 # Number of time steps for forward propagation
        FW_t, FW_Δr_Traj, FW_u_RSW        = FP_L_Δr(Mu, Epsilon, Sigma, t0, tf, Theta, ΔvtcaECI, yp0, ys0, NTS_Δr,  NuES0, ODB, EclipseBool, 0.0, RE) 
        
        filename               = "DAComp1"
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
        #NTS_X = Int(1e6) #Fine grid
        NTS_x = 100
        if DimBool
            Buffertime = 2.5e-5
        else 
            Buffertime = 1 # s
        end
        TimeDL = 1 ## REMOVE
        #t_NL_C_Xp, Traj_NL_C_Xp, t_NL_NC_Xp, Traj_NL_NC_Xp, t_NL_Xs, Traj_NL_Xs   = FP_NL_x(Mu, Epsilon, t0_from_dc, tf+Buffertime, tf+Buffertime, theta_f_from_dc, ΔvtcaECI, yp0, yp0Int, ys0, NTS_X, NuES0, ODB, TuningP)
        true_tCA_NC, true_mindis_NC, true_tCA_C, true_mindis_C, total_Deltav = FP_NL_x(Mu, Epsilon, t0_from_dc, tf+Buffertime, tf, theta_f_from_dc, ΔvtcaECI, yp0, yp0Int, ys0, NuES0, ODB, TuningP, true, TimeDL)

        if DimBool
            println("Nominal tCA at: ", true_tCA_NC*TimeDL, " s")
            println("Nominal miss distance: ", true_mindis_NC*DistanceDL, " m")

            println("Adjusted tCA at: ", true_tCA_C*TimeDL, " s")
            println("Adjusted miss distance: ", true_mindis_C*DistanceDL, " m")

            println("Total DeltaV: ", total_Deltav * DistanceDL/TimeDL, " m/s")
        else 
            println("Nominal tCA at: ", true_tCA_NC, " s")
            println("Nominal miss distance: ", true_mindis_NC, " m")

            println("Adjusted tCA at: ", true_tCA_C, " s")
            println("Adjusted miss distance: ", true_mindis_C, " m")
            println("Total DeltaV: ", total_Deltav, " m/s")
        end

        # true_tCA_NC_ML, true_mindis_NC_ML, true_tCA_C_ML, true_mindis_C_ML, total_Deltav_ML = FP_NL_x(Mu, Epsilon, t0_from_dc, tf+Buffertime, tf, theta_f_from_dc, ΔvtcaECI, yp0, yp0Int, ys0, NuES0, ODB, TuningP, false,  TimeDL)
        # if DimBool
        #     println("Nominal tCA (Matlab) at: ", true_tCA_NC_ML*TimeDL, " s")
        #     println("Nominal miss distance (Matlab): ", true_mindis_NC_ML*DistanceDL, " m")

        #     println("Adjusted tCA (Matlab) at: ", true_tCA_C_ML*TimeDL, " s")
        #     println("Adjusted miss distance (Matlab): ", true_mindis_C_ML*DistanceDL, " m")

        #     println("Total DeltaV (Matlab): ", total_Deltav_ML * DistanceDL/TimeDL, " m/s")
        # else 
        #     println("Nominal tCA (Matlab) at: ", true_tCA_NC_ML, " s")
        #     println("Nominal miss distance (Matlab): ", true_mindis_NC_ML, " m")

        #     println("Adjusted tCA (Matlab) at: ", true_tCA_C_ML, " s")
        #     println("Adjusted miss distance (Matlab): ", true_mindis_C_ML, " m")

        #     println("Total DeltaV (Matlab): ", total_Deltav_ML, " m/s")

        # end

        ## Note: tf=0 works in first-order approximation, but the method tracks states at the 'real' adjusted tca
        ## Even if tca evolution itself is not extracted from the linearised solution (could be extracted from full non-linear propagation)
        ## The code for the non-linear propagation has been removed from this code version for clarity
    end
    return t0_from_dc, total_Deltav
end 

