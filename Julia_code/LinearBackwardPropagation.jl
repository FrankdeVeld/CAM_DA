# ============================================================================
# File: LinearBackwardPropagation.jl (Adjusted Section)
# ============================================================================
using DifferentialEquations, StaticArrays, FFTW

# Assume Conversions.jl (or elements.jl) is included above
# Assume MinTimeDynamics.jl provides generic versions of IntδΔr!, GetCon, ∂r∂yNoT etc.
include("MinTimeDynamics.jl") # Provides IntδΔr!, GetCon, B1B3 etc.
include("Conversions.jl") # Provides E2BU, B2EU, SA, etc.
# const TOL_ZERO_NORM = 1e-12 # Defined above
"""
MODIFIED: IntegrateCircleδΔr
Integrates the linearized relative motion dynamics backward in time.
Now also computes and returns FFT coefficients of the scaled B-plane wavefront.

Returns:
- `tIntTot::Matrix{Float64}`: (NumTheta x NumPoints) time points.
- `δΔrIntTot::Array{T, 3}`: (NumTheta x 3 x NumPoints) integrated δΔr (ECI, physical / epsilon).
- `uOptTot::Array{T, 3}`: (NumTheta x 3 x NumPoints) optimal control (RSW).
- `B1Ref::SVector{3, T}`
- `B3U::SVector{3, T}`
- `fft_coeffs_b1_scaled_t::Matrix{Complex{Float64}}`: (NumTheta x NumPoints) p_k,1(t)
- `fft_coeffs_b1_scaled_dt_t::Matrix{Complex{Float64}}`: (NumTheta x NumPoints) d(p_k,1)/dt (numerically derived)
- `fft_coeffs_b2_scaled_t::Matrix{Complex{Float64}}`: (NumTheta x NumPoints) p_k,2(t)
- `fft_coeffs_b2_scaled_dt_t::Matrix{Complex{Float64}}`: (NumTheta x NumPoints) d(p_k,2)/dt (numerically derived)
"""
function IntegrateCircleδΔr(NumTheta::Int, NumPoints::Int, t0_physical::Float64, # t0_physical is negative
                           sim_params::NamedTuple)

    # --- Input Unpacking and Type Inference ---
    if t0_physical >= 0
        @warn "IntegrateCircleδΔr: Initial time t0_physical should be negative. Got t0_physical = $t0_physical."
    end
    Mu    = sim_params.Mu
    yp0   = sim_params.yp0
    Δvtca = sim_params.Δvtca
    ODB   = sim_params.ODB
    EclipseBool = sim_params.EclipseBool
    LuxBool = sim_params.LuxBool
    Epsilon = sim_params.Epsilon # Physical thrust
    Sigma = sim_params.Sigma     # Physical safe distance

    T = eltype(Δvtca)
    zero_T = zero(T); one_T = one(T); two_T = T(2.0); pi_T = T(pi)

    if norm(Δvtca) < TOL_ZERO_NORM
        error("IntegrateCircleδΔr Error: norm(Δvtca) is near zero.")
    end

    B1Ref, B3U = B1B3(Δvtca, yp0, Mu)

    local OTE::SMatrix{3, 6, T, 18}
    try
        OTE = ∂r∂yNoT(Mu, yp0)
    catch e
        error("IntegrateCircleδΔr requires a generic function '∂r∂yNoT(Mu, yp0)' Error: $e")
    end

    tspan_ode = (0.0, t0_physical) # Physical time for ODE solver (integrates backward)
    δΔr0_vec = zeros(T, 3)
    STA = 1e-10; STR = 1e-10
    solver = Tsit5()
    
   # Enforce a LinRange for save_points_physical
    # These are physical time points, from 0 down to t0_physical (most negative)
    save_points_physical = LinRange(0.0, t0_physical, NumPoints)
    if length(save_points_physical) > NumPoints 
        save_points_physical = LinRange(0.0, t0_physical, NumPoints)
    end


    tIntTot   = zeros(Float64, NumTheta, NumPoints)
    δΔrIntTot = zeros(T, NumTheta, 3, NumPoints) # This is delta_r_ECI_physical / epsilon
    uOptTot   = zeros(T, NumTheta, 3, NumPoints)
    
    # Store Ar_B_scaled components temporarily to compute FFTs per time step
    Ar_B1_scaled_profile_vs_theta = zeros(Float64, NumTheta)
    Ar_B2_scaled_profile_vs_theta = zeros(Float64, NumTheta)

    fft_coeffs_b1_scaled_t = zeros(Complex{Float64}, NumTheta, NumPoints)
    fft_coeffs_b2_scaled_t = zeros(Complex{Float64}, NumTheta, NumPoints)
    
    theta_ang_grid_rad = range(0.0, stop=2pi, length=NumTheta + 1)[1:NumTheta]

    println("Starting backward integration for $NumTheta angles...")
    for n_angle_idx in 0:(NumTheta-1) # Corresponds to row in output arrays
        theta_float = n_angle_idx * 2.0 * pi / NumTheta
        theta_T::T = T(theta_float)

        cosT_T = cos(theta_T); sinT_T = sin(theta_T)
        NBPlane_T = SA[cosT_T, sinT_T, zero_T]
        NECI_T = B2EU(NBPlane_T, B1Ref, B3U)

        ode_params_i = (
            Mu = Mu, yp0 = yp0, Δvtca = Δvtca, NTheta = NECI_T, OTE = OTE,
            NuES0 = sim_params.NuES0, ODB = ODB, EclipseBool = EclipseBool, LuxBool = sim_params.LuxBool, TuningP = sim_params.TuningP,
            Re_const = sim_params.Re_const,
            forbidden_intervals = sim_params.forbidden_intervals
            
        )

        prob = ODEProblem((dX, X, p, t) -> IntδΔr!(dX, X, p, t), δΔr0_vec, tspan_ode, ode_params_i);
        sol = solve(prob, solver, dtmax = abs(t0_physical)/100, abstol = STA, reltol = STR, saveat = save_points_physical, progress=false)

        if sol.retcode ∉ (SciMLBase.ReturnCode.Success,SciMLBase.ReturnCode.Terminated)
             @warn "ODE solver failed for angle index $n_angle_idx. Filling with NaN."
             nan_val_T = T(NaN)
             tIntTot[n_angle_idx+1, :] .= NaN
             δΔrIntTot[n_angle_idx+1, :, :] .= nan_val_T
             uOptTot[n_angle_idx+1, :, :] .= nan_val_T
             # FFT coeffs will be NaN later if tIntTot is NaN
             continue
        end
        
        # sol.t is physical time, sol.u is delta_r_ECI_physical / epsilon
        num_sol_points = length(sol.t)
        
        # Store results (handle length mismatch with pre-allocated NumPoints if saveat was adjusted)
        if num_sol_points != NumPoints
            @warn "Solver for angle $n_angle_idx returned $num_sol_points points, expected $NumPoints. Check saveat logic."
            # For now, assume they match or take min length
            len_to_copy = min(num_sol_points, NumPoints)
            tIntTot[n_angle_idx+1, 1:len_to_copy] = sol.t[1:len_to_copy]
            δΔrIntTot[n_angle_idx+1, :, 1:len_to_copy] = sol[1:3, 1:len_to_copy]
        else
            tIntTot[n_angle_idx+1, :] = sol.t
            δΔrIntTot[n_angle_idx+1, :, :] = sol[1:3,:]
        end

        # Calculate and store control (T)
        uOpt_temp_for_angle = zeros(T, 3, NumPoints)
        nan_val_T = T(NaN)
        fill!(uOpt_temp_for_angle, nan_val_T)

        try
            for j_time_idx = 1:num_sol_points # Iterate over time points for this solution
                if j_time_idx <= NumPoints # Ensure we don't write out of bounds for uOpt_temp_for_angle
                    u_rsw, _ = GetCon(ode_params_i, sol.t[j_time_idx])
                    uOpt_temp_for_angle[:, j_time_idx] .= u_rsw
                end
            end
            len_to_copy_u = min(num_sol_points, NumPoints)
            uOptTot[n_angle_idx+1, :, 1:len_to_copy_u] = uOpt_temp_for_angle[:, 1:len_to_copy_u]
            if num_sol_points < NumPoints; uOptTot[n_angle_idx+1, :, (num_sol_points+1):end] .= nan_val_T; end
        catch e
             @warn "Failed to calculate control for angle index $n_angle_idx. Error: $e"
             uOptTot[n_angle_idx+1, :, :] .= nan_val_T
        end
    end # End loop over angles

    # --- Post-process ALL integrated trajectories to get FFT coefficients vs time ---
    # At this point, tIntTot[1,:] should contain the common time grid (physical, 0 to t0_physical)
    # And δΔrIntTot contains all ECI deviations (phys/eps)
    
    common_t_grid_physical = vec(tIntTot[1,:]) # Assuming all are same
    if any(isnan, common_t_grid_physical)
        @warn "NaNs found in common_t_grid_physical. FFT calculation might fail or be inaccurate."
        # Potentially find first NaN and truncate, or handle error
        first_nan_idx = findfirst(isnan, common_t_grid_physical)
        if !isnothing(first_nan_idx) && first_nan_idx > 1
            common_t_grid_physical = common_t_grid_physical[1:first_nan_idx-1]
            # Also truncate δΔrIntTot and output FFT matrices accordingly
            # This means NumPoints effectively changes. For simplicity, error out or ensure no NaNs.
        else # All NaN or first is NaN
             error("Cannot proceed with FFT due to NaNs in time grid.")
        end
    end
    num_actual_points_for_fft = length(common_t_grid_physical)
    # Resize output FFT matrices if num_actual_points_for_fft < NumPoints
    if num_actual_points_for_fft < NumPoints
        fft_coeffs_b1_scaled_t = zeros(Complex{Float64}, NumTheta, num_actual_points_for_fft)
        fft_coeffs_b2_scaled_t = zeros(Complex{Float64}, NumTheta, num_actual_points_for_fft)
    end


    scaling_pert_b_plane = Epsilon / Sigma

    for i_time in 1:num_actual_points_for_fft
        # For each time point, gather the Ar_B_scaled profile vs theta
        for j_angle in 1:NumTheta
            theta_val = theta_ang_grid_rad[j_angle] # This is theta_k from paper, or just theta
            N_B_x = cos(theta_val)
            N_B_y = sin(theta_val)

            # δΔr_ECI_phys_over_eps for this (theta_val, common_t_grid_physical[i_time])
            # The j_angle here corresponds to the initial N(theta) direction for that trajectory
            # So δΔrIntTot[j_angle, :, i_time] is the ECI deviation for the trajectory
            # that was started with initial condition corresponding to theta_ang_grid_rad[j_angle].
            # This is Ar(t, theta_k) in the paper's notation if theta_k is the initial B-plane angle.
            
            pert_ECI_div_eps_vec = SVector{3,T}(
                δΔrIntTot[j_angle, 1, i_time],
                δΔrIntTot[j_angle, 2, i_time],
                δΔrIntTot[j_angle, 3, i_time]
            )
            if any(isnan, pert_ECI_div_eps_vec)
                 Ar_B1_scaled_profile_vs_theta[j_angle] = NaN
                 Ar_B2_scaled_profile_vs_theta[j_angle] = NaN
                 continue
            end
            
            pert_B_plane_div_eps_vec = E2BU(pert_ECI_div_eps_vec, B1Ref, B3U) # Type SVector{3,T}
            
            # Ar_B_scaled = N_B(theta) + (Epsilon/Sigma) * Ar_B_pert(t,theta)
            # Ar_B_pert(t,theta) = E2BU( δΔr_ECI_phys_over_eps(t,theta) , B1Ref, B3U)
            Ar_B1_scaled_profile_vs_theta[j_angle] = N_B_x + scaling_pert_b_plane * pert_B_plane_div_eps_vec[1]
            Ar_B2_scaled_profile_vs_theta[j_angle] = N_B_y + scaling_pert_b_plane * pert_B_plane_div_eps_vec[2]
        end

        if any(isnan, Ar_B1_scaled_profile_vs_theta) || any(isnan, Ar_B2_scaled_profile_vs_theta)
            @warn "NaNs in B-plane profile at time index $i_time. FFT coeffs will be NaN."
            fft_coeffs_b1_scaled_t[:, i_time] .= NaN
            fft_coeffs_b2_scaled_t[:, i_time] .= NaN
        else
            # Compute FFT for this time slice

            fft_b1_raw = fft(Ar_B1_scaled_profile_vs_theta) ./ NumTheta
            fft_b2_raw = fft(Ar_B2_scaled_profile_vs_theta) ./ NumTheta
            
            fft_coeffs_b1_scaled_t[:, i_time] = fftshift(fft_b1_raw)
            fft_coeffs_b2_scaled_t[:, i_time] = fftshift(fft_b2_raw)
        end
    end

    # Numerically differentiate fft_coeffs_b2_scaled_t to get d(p_k,2)/dt
    fft_coeffs_b2_scaled_dt_t = zeros(Complex{Float64}, NumTheta, num_actual_points_for_fft)
    t_interp_for_deriv = reverse(common_t_grid_physical) 
    if num_actual_points_for_fft > 1
        # Check if t_interp_for_deriv can be made into a range
        local t_interp_range_for_deriv :: AbstractRange{Float64}
        is_uniform_deriv_grid = true
        if num_actual_points_for_fft > 1
            first_step_deriv = t_interp_for_deriv[2] - t_interp_for_deriv[1] # Should be positive
            for i_idx in 2:(num_actual_points_for_fft - 1)
                if !isapprox(t_interp_for_deriv[i_idx+1] - t_interp_for_deriv[i_idx], first_step_deriv, atol=1e-9 * abs(first_step_deriv) + 1e-12)
                    is_uniform_deriv_grid = false
                    break
                end
            end
        end

        if is_uniform_deriv_grid
            t_interp_range_for_deriv = range(t_interp_for_deriv[1], stop=t_interp_for_deriv[end], length=num_actual_points_for_fft)
        else
            @warn "Time grid for pk2 derivative calculation is not perfectly uniform. Spline differentiation might be less accurate or have issues if it expects a range."
            t_interp_range_for_deriv = t_interp_for_deriv # Use the vector if not uniform; might error or be less accurate
            error("Time grid for pk2 derivative calculation must be uniform to form an AbstractRange for reliable spline differentiation.")

        end


        for k_coeff_idx in 1:NumTheta # For each k-th coefficient
            p_k2_real_vs_t_reversed = real.(vec(fft_coeffs_b2_scaled_t[k_coeff_idx, end:-1:1]))
            p_k2_imag_vs_t_reversed = imag.(vec(fft_coeffs_b2_scaled_t[k_coeff_idx, end:-1:1]))

            if any(isnan, p_k2_real_vs_t_reversed) || any(isnan, p_k2_imag_vs_t_reversed)
                fft_coeffs_b2_scaled_dt_t[k_coeff_idx, :] .= NaN
                continue
            end

            # Create splines
            spl_real = cubic_spline_interpolation(t_interp_range_for_deriv, p_k2_real_vs_t_reversed, extrapolation_bc=Throw())
            spl_imag = cubic_spline_interpolation(t_interp_range_for_deriv, p_k2_imag_vs_t_reversed, extrapolation_bc=Throw())

            # Differentiate the splines and evaluate at the original (reversed) time points
            # Interpolations.gradient gives a tuple (even for 1D). We need the first element.
            dpk2_real_dt_interp = [Interpolations.gradient(spl_real, t_val)[1] for t_val in t_interp_for_deriv]
            dpk2_imag_dt_interp = [Interpolations.gradient(spl_imag, t_val)[1] for t_val in t_interp_for_deriv]
            
            fft_coeffs_b2_scaled_dt_t[k_coeff_idx, :] = reverse(dpk2_real_dt_interp .+ im .* dpk2_imag_dt_interp)
        end
    else 
        fft_coeffs_b2_scaled_dt_t .= 0.0 
    end

    
    # Numerically differentiate fft_coeffs_b1_scaled_t to get d(p_k,2)/dt
    fft_coeffs_b1_scaled_dt_t = zeros(Complex{Float64}, NumTheta, num_actual_points_for_fft)
    t_interp_for_deriv = reverse(common_t_grid_physical) 
    if num_actual_points_for_fft > 1
        # Check if t_interp_for_deriv can be made into a range
        #local t_interp_range_for_deriv :: AbstractRange{Float64}
        is_uniform_deriv_grid = true
        if num_actual_points_for_fft > 1
            first_step_deriv = t_interp_for_deriv[2] - t_interp_for_deriv[1] # Should be positive
            for i_idx in 2:(num_actual_points_for_fft - 1)
                if !isapprox(t_interp_for_deriv[i_idx+1] - t_interp_for_deriv[i_idx], first_step_deriv, atol=1e-9 * abs(first_step_deriv) + 1e-12)
                    is_uniform_deriv_grid = false
                    break
                end
            end
        end

        if is_uniform_deriv_grid
            t_interp_range_for_deriv = range(t_interp_for_deriv[1], stop=t_interp_for_deriv[end], length=num_actual_points_for_fft)
        else
            @warn "Time grid for pk2 derivative calculation is not perfectly uniform. Spline differentiation might be less accurate or have issues if it expects a range."
            t_interp_range_for_deriv = t_interp_for_deriv # Use the vector if not uniform; might error or be less accurate
            error("Time grid for pk2 derivative calculation must be uniform to form an AbstractRange for reliable spline differentiation.")

        end


        for k_coeff_idx in 1:NumTheta # For each k-th coefficient
            p_k1_real_vs_t_reversed = real.(vec(fft_coeffs_b1_scaled_t[k_coeff_idx, end:-1:1]))
            p_k1_imag_vs_t_reversed = imag.(vec(fft_coeffs_b1_scaled_t[k_coeff_idx, end:-1:1]))

            if any(isnan, p_k1_real_vs_t_reversed) || any(isnan, p_k1_imag_vs_t_reversed)
                fft_coeffs_b1_scaled_dt_t[k_coeff_idx, :] .= NaN
                continue
            end

            # Create splines
            spl_real = cubic_spline_interpolation(t_interp_range_for_deriv, p_k1_real_vs_t_reversed, extrapolation_bc=Throw())
            spl_imag = cubic_spline_interpolation(t_interp_range_for_deriv, p_k1_imag_vs_t_reversed, extrapolation_bc=Throw())
            # Interpolations.gradient gives a tuple (even for 1D). We need the first element.
            dpk1_real_dt_interp = [Interpolations.gradient(spl_real, t_val)[1] for t_val in t_interp_for_deriv]
            dpk1_imag_dt_interp = [Interpolations.gradient(spl_imag, t_val)[1] for t_val in t_interp_for_deriv]
            
            fft_coeffs_b1_scaled_dt_t[k_coeff_idx, :] = reverse(dpk1_real_dt_interp .+ im .* dpk1_imag_dt_interp)
        end
    else 
        fft_coeffs_b1_scaled_dt_t .= 0.0 
    end
    
    # Trim tIntTot and δΔrIntTot if num_actual_points_for_fft < NumPoints
    if num_actual_points_for_fft < NumPoints
        tIntTot_final = tIntTot[:, 1:num_actual_points_for_fft]
        δΔrIntTot_final = δΔrIntTot[:, :, 1:num_actual_points_for_fft]
        uOptTot_final = uOptTot[:, :, 1:num_actual_points_for_fft]
    else
        tIntTot_final = tIntTot
        δΔrIntTot_final = δΔrIntTot
        uOptTot_final = uOptTot
    end


    println("Backward integration and FFT processing finished.")
    return tIntTot_final, δΔrIntTot_final, uOptTot_final, B1Ref, B3U, fft_coeffs_b1_scaled_t, fft_coeffs_b1_scaled_dt_t, fft_coeffs_b2_scaled_t, fft_coeffs_b2_scaled_dt_t
end


"""
    δΔrScaling(NumTheta::Int, NumPoints::Int, δΔrIntTot::Array{T, 3},
               Epsilon::T, Sigma::T,
               B1Ref::SVector{3, T}, B3U::SVector{3, T}) where {T<:Real} -> Array{T, 3}

Scales integrated δΔr trajectories (type T) to the normalized B-plane.
"""
function δΔrScaling(NumTheta::Int, NumPoints::Int, δΔrIntTot::Array{T, 3},
                   Epsilon::T, Sigma::T,
                   B1Ref::SVector{3, T}, B3U::SVector{3, T}) where {T<:Real}

    zero_T = zero(T); two_T = T(2.0); pi_T = T(pi)

    # Use eps(T) for zero check if Sigma could be Dual
    if abs(Sigma) < eps(T); error("δΔrScaling Error: Sigma is near zero."); end
    if size(δΔrIntTot) != (NumTheta, 3, NumPoints); error("δΔrScaling Error: Input dimensions mismatch."); end

    # Calculate angles using Float64, then convert to T for trig/vector construction
    theta_vals_float = [n * 2.0 * pi / NumTheta for n in 0:(NumTheta-1)]
    N_B_unit_vectors = [SA[T(cos(th)), T(sin(th)), zero_T] for th in theta_vals_float] # List of SVector{3, T}

    # Allocate output arrays with type T
    ΔR_ESS_F     = zeros(T, NumTheta, NumPoints, 3)
    ΔR_ESS_F_app = zeros(T, NumTheta + 1, NumPoints, 3)

    println("Scaling integrated deviations to B-plane reachable set...")
    scaling_factor::T = Epsilon / Sigma

    for n in 0:(NumTheta - 1)
        N_B::SVector{3, T} = N_B_unit_vectors[n+1]

        for i in 1:NumPoints
            # Extract δΔr(t) as SVector{T} for calculations
            # Use @view for efficiency, then convert to SVector{T}
            δΔr_eci = SVector{3, T}(@view δΔrIntTot[n+1, :, i])

            Δr_S_E = scaling_factor * δΔr_eci # SVector{3, T}
            # E2BU takes SVector{3, T}, returns SVector{3, T}
            Δr_S_B = E2BU(Δr_S_E, B1Ref, B3U)

            # Assign result (SVector{3, T}) back to standard Array{T, 3} slice
            ΔR_ESS_F[n+1, i, :] = N_B + Δr_S_B
        end
    end

    ΔR_ESS_F_app[1:NumTheta, :, :] = ΔR_ESS_F
    # Use view for efficient assignment of the wrapped-around point
    ΔR_ESS_F_app[NumTheta + 1, :, :] = @view ΔR_ESS_F[1, :, :]

    println("Scaling finished.")
    return ΔR_ESS_F_app; # Returns Array{T, 3}
end