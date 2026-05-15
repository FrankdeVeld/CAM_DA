using DifferentialEquations
using Interpolations # For interpolating FFT coeffs vs time
using FFTW # For reconstructing Ar_B_1 from coeffs
using StaticArrays
using Roots 
using Plots
using SciMLBase

# Wavenumbers for N points, shifted for sum from k=-K to K
# K = (N-1)/2 for N odd.
# fftshift maps indices [1,...,N] to [-floor(N/2)+1, ..., ceil(N/2)] (approx)
# For N odd, K=(N-1)/2. Indices for fftshift output are effectively -K, ..., 0, ..., K
# if the original array was indexed 0 to N-1 for physical space.
# If original array X was length N, fftshift(fft(X))[j] corresponds to wavenumber m_j.
# The standard wavenumbers for an N-point FFT, after fftshift, are typically:
# N even: [-(N/2)+1, ..., -1, 0, 1, ..., N/2] (N/2 is Nyquist, often special)
# N odd:  [-(N-1)/2, ..., -1, 0, 1, ..., (N-1)/2]
# These are the 'k' values in sum p_k e^(ik theta)
function get_wavenumbers(NumTheta::Int)
    if isodd(NumTheta)
        K = (NumTheta - 1) ÷ 2
        return SVector{NumTheta, Int}((-K):K)
    else # NumTheta even
        error("NumTheta must be odd for the current wavenumber setup in DC_FFT.")
    end
end


"""
Reconstruct Ar_B_1_scaled(t_physical, theta_c) from its FFT coefficients p_k,1(t_physical)
"""
function reconstruct_from_fft_coeffs(
    theta_c_rad,#::Float64,
    coeffs_p_k_at_t,#::AbstractVector{Complex{Float64}}, # fftshifted coeffs for a specific time
    wavenumbers_k,#::SVector # fftshifted wavenumbers
    )
    
    val = 0.0 + 0.0im
    for (idx, k) in enumerate(wavenumbers_k)
        val += coeffs_p_k_at_t[idx] * exp(im * k * theta_c_rad)
    end
    # Ar_B should be real. Small imag part due to numerics.
    if abs(imag(val)) > 1e-9
        @warn "Significant imaginary part in reconstructed value: $(imag(val))"
    end
    return real(val)
end

"""
ODE function for theta_c using FFT coefficients, for a HORIZONTAL constraint (Ar_B2 = const).
`t_ode` is physical time (decreasing from 0 to t0_physical).
"""
function dc_ode_fft_horizontal!(dXc, Xc, p, t_ode_physical)
    theta_c_rad = Xc[1]
    
    # Unpack parameters for horizontal constraint (Ar_B2 = const)
    itp_pk2_coeffs_vs_t = p.itp_pk2_coeffs_vs_t
    itp_dpk2_dt_coeffs_vs_t = p.itp_dpk2_dt_coeffs_vs_t
    wavenumbers_k = p.wavenumbers_k
    
    current_pk2_coeffs = itp_pk2_coeffs_vs_t(t_ode_physical)
    current_dpk2_dt_coeffs = itp_dpk2_dt_coeffs_vs_t(t_ode_physical)

    # Numerator of Eq. 68: sum_k {d(p_k,2)/dt * e^(ik theta_c)}
    numerator_sum = sum(current_dpk2_dt_coeffs[idx] * exp(im * k * theta_c_rad) for (idx, k) in enumerate(wavenumbers_k))
    
    # Denominator of Eq. 68: sum_k {i*k*p_k,2 * e^(ik theta_c)}
    denominator_sum = sum(im * k * current_pk2_coeffs[idx] * exp(im * k * theta_c_rad) for (idx, k) in enumerate(wavenumbers_k))

    if abs(denominator_sum) < 1e-9
        dXc[1] = 0.0
    else
        val = -numerator_sum / denominator_sum
        dXc[1] = real(val)
    end
end

"""
ODE function for theta_c using FFT coefficients, for a VERTICAL constraint (Ar_B1 = const).
`t_ode` is physical time (decreasing from 0 to t0_physical).
"""
function dc_ode_fft_vertical!(dXc, Xc, p, t_ode_physical)
    theta_c_rad = Xc[1]
    
    # Unpack parameters for vertical constraint (Ar_B1 = const)
    itp_pk1_coeffs_vs_t = p.itp_pk1_coeffs_vs_t
    itp_dpk1_dt_coeffs_vs_t = p.itp_dpk1_dt_coeffs_vs_t
    wavenumbers_k = p.wavenumbers_k

    current_pk1_coeffs = itp_pk1_coeffs_vs_t(t_ode_physical)
    current_dpk1_dt_coeffs = itp_dpk1_dt_coeffs_vs_t(t_ode_physical)

    # Numerator: sum_k {d(p_k,1)/dt * e^(ik theta_c)}
    numerator_sum = sum(current_dpk1_dt_coeffs[idx] * exp(im * k * theta_c_rad) for (idx, k) in enumerate(wavenumbers_k))
    
    # Denominator: sum_k {i*k*p_k,1 * e^(ik theta_c)}
    denominator_sum = sum(im * k * current_pk1_coeffs[idx] * exp(im * k * theta_c_rad) for (idx, k) in enumerate(wavenumbers_k))
    
    if abs(denominator_sum) < 1e-9
        dXc[1] = 0.0
    else
        val = -numerator_sum / denominator_sum
        dXc[1] = real(val)
    end
end

"""
Top-level wrapper for Differential Continuation.
Tries the standard horizontal approach first. If both initial branches fail,
it attempts a fail-over using a vertical approach.

Returns a vector of successful solutions. A successful solution is a tuple:
(t0_physical, theta_f_rad, retcode, solution_time_points, solution_states)
"""
function run_dc_fft_with_failover(
    InitCon_B_scaled::SVector{3,Float64},
    t_physical_grid::AbstractVector{Float64},
    fft_coeffs_b1_t_series::Matrix{Complex{Float64}},
    fft_coeffs_b1_dt_t_series::Matrix{Complex{Float64}},
    fft_coeffs_b2_t_series::Matrix{Complex{Float64}},
    fft_coeffs_b2_dt_t_series::Matrix{Complex{Float64}},
    NumTheta::Int
)
    successful_solutions = []
    args = (
        InitCon_B_scaled, t_physical_grid,
        fft_coeffs_b1_t_series, fft_coeffs_b1_dt_t_series,
        fft_coeffs_b2_t_series, fft_coeffs_b2_dt_t_series,
        NumTheta
    )

    # --- Step 1: Try Horizontal approach ---
    println("Attempting horizontal continuation...")
    for branch in 1:2
        println("  - Trying branch $branch...")
        try
            kwargs = (initial_theta_branch=branch, orientation=:horizontal)
            t0, theta_f, retcode, sol_t, sol_u = solve_dc_fft(args...; kwargs...)

            if retcode == SciMLBase.ReturnCode.Terminated && !isnan(t0)
                println("    Success on horizontal branch $branch. t0 = $t0")
                push!(successful_solutions, (t0, theta_f, retcode, sol_t, sol_u))
            else
                println("    Horizontal branch $branch did not yield a valid solution (retcode: $retcode).")
            end
        catch e
            println("    Error on horizontal branch $branch: $e")
        end
    end

    # --- Step 2: If horizontal failed, try vertical (fail-safe) ---
    if isempty(successful_solutions)
        println("\nHorizontal continuation failed. Attempting vertical fail-over...")
        for branch in 1:2
            println("  - Trying branch $branch...")
            try
                kwargs = (initial_theta_branch=branch, orientation=:vertical)
                t0, theta_f, retcode, sol_t, sol_u = solve_dc_fft(args...; kwargs...)

                if retcode == SciMLBase.ReturnCode.Terminated && !isnan(t0)
                    println("    Success on vertical branch $branch. t0 = $t0")
                    push!(successful_solutions, (t0, theta_f, retcode, sol_t, sol_u))
                else
                    println("    Vertical branch $branch did not yield a valid solution (retcode: $retcode).")
                end
            catch e
                println("    Error on vertical branch $branch: $e")
            end
        end
    end
    
    if isempty(successful_solutions)
        println("\nAll continuation attempts failed.")
    else
        println("\nFound $(length(successful_solutions)) valid solution(s).")
    end

    return successful_solutions
end


"""
    solve_dc_fft(
        InitCon_B_scaled::SVector{3,Float64}, # Target [b1, b2, 0]
        t_physical_grid::AbstractVector{Float64},
        fft_coeffs_b1_t_series::Matrix{Complex{Float64}},
        fft_coeffs_b1_dt_t_series::Matrix{Complex{Float64}}, # NEW: for vertical continuation
        fft_coeffs_b2_t_series::Matrix{Complex{Float64}},
        fft_coeffs_b2_dt_t_series::Matrix{Complex{Float64}},
        NumTheta::Int;
        initial_theta_branch::Int = 1,
        orientation::Symbol = :horizontal # NEW: :horizontal or :vertical
    ) -> Tuple{Float64, Float64, Symbol, Vector{Float64}, Vector{Vector{Float64}}}
"""
function solve_dc_fft(
    InitCon_B_scaled::SVector{3,Float64},
    t_physical_grid::AbstractVector{Float64},
    fft_coeffs_b1_t_series::Matrix{Complex{Float64}},
    fft_coeffs_b1_dt_t_series::Matrix{Complex{Float64}},
    fft_coeffs_b2_t_series::Matrix{Complex{Float64}},
    fft_coeffs_b2_dt_t_series::Matrix{Complex{Float64}},
    NumTheta::Int;
    initial_theta_branch::Int = 1,
    orientation::Symbol = :horizontal
)
    if !isodd(NumTheta)
        error("solve_dc_fft currently requires NumTheta to be odd.")
    end
    wavenumbers_k = get_wavenumbers(NumTheta)
    Debug_mode = false

    if length(t_physical_grid) < 2
        @warn "t_physical_grid has fewer than 2 points. Cannot solve."
        return NaN, NaN, :Failure_TimeGridEmpty, Float64[], Vector{Float64}[]
    end

    t_interp_range = range(t_physical_grid[end], stop=t_physical_grid[1], length=length(t_physical_grid))

    itp_pk1_real_list = Vector{Any}(undef, NumTheta)
    itp_pk1_imag_list = Vector{Any}(undef, NumTheta)
    itp_dpk1_dt_real_list = Vector{Any}(undef, NumTheta)
    itp_dpk1_dt_imag_list = Vector{Any}(undef, NumTheta)
    itp_pk2_real_list = Vector{Any}(undef, NumTheta)
    itp_pk2_imag_list = Vector{Any}(undef, NumTheta)
    itp_dpk2_dt_real_list = Vector{Any}(undef, NumTheta)
    itp_dpk2_dt_imag_list = Vector{Any}(undef, NumTheta)

    for k_idx in 1:NumTheta
        pk1_series = vec(fft_coeffs_b1_t_series[k_idx, end:-1:1])
        dpk1_dt_series = vec(fft_coeffs_b1_dt_t_series[k_idx, end:-1:1])
        pk2_series = vec(fft_coeffs_b2_t_series[k_idx, end:-1:1])
        dpk2_dt_series = vec(fft_coeffs_b2_dt_t_series[k_idx, end:-1:1])
        
        itp_pk1_real_list[k_idx] = cubic_spline_interpolation(t_interp_range, real.(pk1_series), extrapolation_bc=Throw())
        itp_pk1_imag_list[k_idx] = cubic_spline_interpolation(t_interp_range, imag.(pk1_series), extrapolation_bc=Throw())
        itp_dpk1_dt_real_list[k_idx] = cubic_spline_interpolation(t_interp_range, real.(dpk1_dt_series), extrapolation_bc=Throw())
        itp_dpk1_dt_imag_list[k_idx] = cubic_spline_interpolation(t_interp_range, imag.(dpk1_dt_series), extrapolation_bc=Throw())

        itp_pk2_real_list[k_idx] = cubic_spline_interpolation(t_interp_range, real.(pk2_series), extrapolation_bc=Throw())
        itp_pk2_imag_list[k_idx] = cubic_spline_interpolation(t_interp_range, imag.(pk2_series), extrapolation_bc=Throw())
        itp_dpk2_dt_real_list[k_idx] = cubic_spline_interpolation(t_interp_range, real.(dpk2_dt_series), extrapolation_bc=Throw())
        itp_dpk2_dt_imag_list[k_idx] = cubic_spline_interpolation(t_interp_range, imag.(dpk2_dt_series), extrapolation_bc=Throw())
    end

    get_all_pk1_at_t(t) = SVector{NumTheta, Complex{eltype(t)}}(itp_pk1_real_list[i](t) + im*itp_pk1_imag_list[i](t) for i=1:NumTheta)
    get_all_dpk1_dt_at_t(t) = SVector{NumTheta, Complex{eltype(t)}}(itp_dpk1_dt_real_list[i](t) + im*itp_dpk1_dt_imag_list[i](t) for i=1:NumTheta)
    get_all_pk2_at_t(t) = SVector{NumTheta, Complex{eltype(t)}}(itp_pk2_real_list[i](t) + im*itp_pk2_imag_list[i](t) for i=1:NumTheta)
    get_all_dpk2_dt_at_t(t) = SVector{NumTheta, Complex{eltype(t)}}(itp_dpk2_dt_real_list[i](t) + im*itp_dpk2_dt_imag_list[i](t) for i=1:NumTheta)

    # --- Generic setup based on orientation ---
    local ode_func, params_ode, initial_condition_target_func, condition_cb

    if orientation == :horizontal
        ode_func = dc_ode_fft_horizontal!
        params_ode = (itp_pk2_coeffs_vs_t = get_all_pk2_at_t, itp_dpk2_dt_coeffs_vs_t = get_all_dpk2_dt_at_t, wavenumbers_k = wavenumbers_k)
        initial_coeffs_at_t0 = get_all_pk2_at_t(0.0)
        initial_target_val = InitCon_B_scaled[2]
        callback_coeffs_func = get_all_pk1_at_t
        callback_target_val = InitCon_B_scaled[1]
    elseif orientation == :vertical
        ode_func = dc_ode_fft_vertical!
        params_ode = (itp_pk1_coeffs_vs_t = get_all_pk1_at_t, itp_dpk1_dt_coeffs_vs_t = get_all_dpk1_dt_at_t, wavenumbers_k = wavenumbers_k)
        initial_coeffs_at_t0 = get_all_pk1_at_t(0.0)
        initial_target_val = InitCon_B_scaled[1]
        callback_coeffs_func = get_all_pk2_at_t
        callback_target_val = InitCon_B_scaled[2]
    else
        error("Invalid orientation: $orientation. Must be :horizontal or :vertical.")
    end

    initial_condition_target_func(theta) = reconstruct_from_fft_coeffs(theta, initial_coeffs_at_t0, wavenumbers_k) - initial_target_val
    condition_cb = (Xc, t, integrator) -> reconstruct_from_fft_coeffs(Xc[1], callback_coeffs_func(t), wavenumbers_k) - callback_target_val
    
    # --- Find initial theta_c ---
    initial_theta_roots = find_zeros(initial_condition_target_func, 0.0, 2*pi)
    if isempty(initial_theta_roots)
        @warn "DC_FFT ($orientation): No initial theta_c found satisfying the constraint. Target might be out of range at t=0."
        return NaN, NaN, :Failure_InitialTheta_NoRoots, Float64[], Vector{Float64}[]
    end
    sort!(initial_theta_roots)

    theta_c_initial_rad = 0.0
    if initial_theta_branch == 1
        theta_c_initial_rad = initial_theta_roots[1]
    elseif initial_theta_branch == 2
        if length(initial_theta_roots) >= 2
            theta_c_initial_rad = initial_theta_roots[2]
        else
            @warn "DC_FFT ($orientation): Only one initial theta_c root found, cannot select branch 2. Using the first root."
            theta_c_initial_rad = initial_theta_roots[1]
        end
    else
        error("Invalid initial_theta_branch value: $initial_theta_branch")
    end
    theta_c_initial_rad = mod(theta_c_initial_rad, 2*pi)
    
    # --- ODE Problem and Solver ---
    t_max_neg = t_physical_grid[end]
    tspan_ode_physical = (0.0, t_max_neg)
    
    affect_cb! = (integrator) -> terminate!(integrator)
    cb = ContinuousCallback(condition_cb, affect_cb!, abstol=1e-7)
    
    prob = ODEProblem(ode_func, [theta_c_initial_rad], tspan_ode_physical, params_ode)
    sol = solve(prob, Tsit5(), callback=cb, abstol=1e-9, reltol=1e-9, dtmin=1e-6, maxiters=10000, save_everystep=true, dense=false, dtmax=abs(t_max_neg)/200.0)

    # --- Process solution ---
    if sol.retcode == SciMLBase.ReturnCode.Terminated
        t0_physical_event = sol.t[end]
        theta_f_event_rad = sol.u[end][1]
        return t0_physical_event, mod(theta_f_event_rad, 2*pi), sol.retcode, sol.t, sol.u
    elseif sol.retcode == SciMLBase.ReturnCode.Success
        @warn "DC_FFT ($orientation) reached MaxTime (t=$(sol.t[end])) without meeting condition for branch $initial_theta_branch."
        return t_max_neg, mod(sol.u[end][1], 2*pi), sol.retcode, sol.t, sol.u
    else 
        @warn "DC_FFT ($orientation) ODE solve failed with unexpected retcode for branch $initial_theta_branch. Retcode: $(sol.retcode)"
        return NaN, NaN, sol.retcode, sol.t, sol.u
    end
end