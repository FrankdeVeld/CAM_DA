# File: YawSteering.jl
# Purpose: Implements yaw steering logic based on Sun direction, including
#          perfect and rate-limited options.

include("Conversions.jl") # Provides coordinate transforms, safe_normalize etc.
using LinearAlgebra
using Printf # For warnings
using ForwardDiff 
using StaticArrays

const TOL_ZERO = 1e-12 # Tolerance for near-zero checks

"""
    calculate_perfect_yaw(alpha::Float64, beta::Float64) -> Float64

Calculates the "perfect" or unconstrained yaw steering angle (zeta_e or Zeta)
based on Sun azimuth (alpha) and elevation (beta).

Matches the formula: zeta_e = atan2(sin(beta), cos(beta)*sin(alpha))

# Arguments
- `alpha::Float64`: Sun azimuth angle relative to noon (radians), typically in [-π, π].
- `beta::Float64`: Sun elevation angle relative to orbit plane (radians), typically in [-π/2, π/2].

# Returns
- `zeta_e::Float64`: Perfect yaw steering angle (radians), in range [-π, π].
"""
function calculate_perfect_yaw(alpha, beta)
    # atan2(y, x) is generally preferred for handling quadrants correctly.
    # y = sin(beta)
    # x = cos(beta)*sin(alpha)
    zeta_e = atan(sin(beta), cos(beta)*sin(alpha))
    return zeta_e # Output range is [-pi, pi]
end


"""
    rESECI(t::Float64, NuES0::T) where T<:Real -> SVector{3, <Promoted Type>}
Calculates the unit vector from Earth to Sun in ECI coordinates.
Returns SVector. NuES0 can be generic.
"""
function rESECI(t::TimeT, NuES0::N0T) where {TimeT<:Real, N0T<:Real}
    n_Earth = 2.0 * pi / (365.25 * 24.0 * 3600.0) # Float64
    EpsilonS = deg2rad(23.44) # Float64

    # lambda_Sun will be promote_type(TimeT, N0T, Float64)
    lambda_Sun = NuES0 + n_Earth * t # t is now TimeT

    rES_x = cos(lambda_Sun)
    rES_y = sin(lambda_Sun) * cos(EpsilonS) # EpsilonS is Float64, result promotes
    rES_z = sin(lambda_Sun) * sin(EpsilonS)

    rES_vec = SA[rES_x, rES_y, rES_z]
    norm_rES = norm(rES_vec)
    local result_vec

    if norm_rES < TOL_ZERO
        @warn "rESECI: Calculated Sun vector has near-zero norm. Returning [1,0,0]."
        zero_el = zero(eltype(rES_vec))
        one_el  = one(eltype(rES_vec))
        result_vec = SA[one_el, zero_el, zero_el]
    else
        result_vec = rES_vec / norm_rES
    end
    return result_vec
end
function AlphaCalc(x_eq::AbstractVector{EqT}, t::TimeT, NuES0::N0T, Mu::MuParamT) where {EqT<:Real, TimeT<:Real, N0T<:Real, MuParamT<:Real}
    if length(x_eq) != 6; error("AlphaCalc: Input x_eq must have 6 elements."); end

    # rESECI needs to handle t::TimeT. If its signature is rESECI(t::Float64, ...),
    # then TimeT MUST be Float64 or you need to convert t, or make rESECI take t::TimeT.
    # Assuming rESECI's t argument is also made generic or TimeT is Float64 here.
    # For now, let's assume rESECI was also updated for t::TimeT or TimeT happens to be Float64 in this call stack.
    # If rESECI strictly needs Float64 for t, and TimeT can be Dual:
    # rES_unit = rESECI(Float64(value(t)), NuES0) # If t is Dual, strips dual part for rESECI
    # This ^ would break AD if rESECI's result depends on t's derivative.
    # Best is to make rESECI also accept t::TimeT
    rES_unit = rESECI(t, NuES0) # Assumes rESECI can handle t::TimeT

    x_cart = equin2cart(x_eq, Mu) # Returns SVector{6, promote(EqT,MuParamT)}
    Rp = x_cart[SA[1, 2, 3]]
    Vp = x_cart[SA[4, 5, 6]]
    h = cross(Rp, Vp)

    if norm(h) < TOL_ZERO
        error("AlphaCalc: Cannot calculate orbit normal, angular momentum is near zero.")
    end

    Rp_OPF = E2OU(Rp, h, rES_unit)


    Alpha = atan(Rp_OPF[2], Rp_OPF[1])
    return Alpha
end
"""
    BetaCalc(x_eq::Vector{Float64}, t::Float64, NuES0::Float64, Mu::Float64) -> Float64

Calculates the Sun elevation angle 'beta' relative to the orbital plane.
(Implementation largely unchanged, added safety check)
"""
function BetaCalc(x_eq::AbstractVector{EqT}, t::TimeT, NuES0::N0T, Mu::MuParamT) where {EqT<:Real, TimeT<:Real, N0T<:Real, MuParamT<:Real}
    if length(x_eq) != 6; error("BetaCalc: Input x_eq must have 6 elements."); end

    rES_unit = rESECI(t, NuES0) # Assumes rESECI can handle t::TimeT
    x_cart = equin2cart(x_eq, Mu)
    Rp = x_cart[SA[1, 2, 3]]
    Vp = x_cart[SA[4, 5, 6]]
    h = cross(Rp, Vp)

    if norm(h) < TOL_ZERO
         error("BetaCalc: Cannot calculate orbit normal, angular momentum is near zero.")
    end
    h_unit = safe_normalize(h, "h in BetaCalc")
    sin_beta = dot(h_unit, rES_unit)
    sin_beta_clamped = clamp(sin_beta, -one(sin_beta), one(sin_beta))
    Beta = asin(sin_beta_clamped)
    # @show h_unit 
    # @show rES_unit
    # @show Beta*180/pi
    # sleep(1)
    return Beta
end

"""
    apply_yaw_rate_limit(zeta_perfect::Float64, alpha::Float64, beta::Float64, beta0::Float64, zeta_rate_yaw::Float64) -> Float64

Applies the full yaw steering logic, including rate limiting and the
comparison between perfect and rate-limited profiles when |beta| < beta0.

Implements the logic from the "Yaw steering motion", "Rate limited yaw angle",
and "To sum-up" sections of the provided PDF.

# Arguments
- `zeta_perfect::Float64`: Perfect yaw angle (psi1) (radians).
- `alpha::Float64`: Sun azimuth angle relative to noon (radians), in [-π, π].
- `beta::Float64`: Sun elevation angle (radians).
- `beta0::Float64`: Threshold Sun elevation angle for rate limiting (radians, > 0).
- `zeta_rate_yaw::Float64`: Kyaw = yaw_rate_max / omega (rad/rad, > 0).

# Returns
- `zeta_final::Float64`: The final applied yaw angle (psi) (radians).
"""
function apply_yaw_rate_limit(zeta_perfect::ZPT, alpha::AT, beta::BT, beta0::Float64, zeta_rate_yaw::Float64) where {ZPT<:Real, AT<:Real, BT<:Real}
    # --- Input Validation ---
    if zeta_rate_yaw <= TOL_ZERO # TOL_ZERO is Float64
        error("apply_yaw_rate_limit: zeta_rate_yaw (Kyaw) must be positive.")
    end

    # Promoted pi/2
    pi_half = promote_type(ZPT, AT, BT)(pi/2)
    pi_val = promote_type(ZPT, AT, BT)(pi)


    zeta2_positive_beta::promote_type(AT, Float64) = zero(promote_type(AT,Float64))
    rate_term = zeta_rate_yaw # Float64

    if -pi_val <= alpha < -pi_half # Comparisons work with Duals
        zeta2_positive_beta = pi_half + rate_term * (alpha + pi_val)
    elseif -pi_half <= alpha < pi_half
        zeta2_positive_beta = pi_half - rate_term * alpha
    elseif pi_half <= alpha <= pi_val
        zeta2_positive_beta = pi_half - rate_term * (pi_val - alpha)
    else
        alpha_wrapped = mod(alpha + pi_val, 2*pi_val) - pi_val
         @warn "apply_yaw_rate_limit: alpha ($alpha -> $alpha_wrapped) outside expected. Using wrapped."
        if -pi_val <= alpha_wrapped < -pi_half
            zeta2_positive_beta = pi_half + rate_term * (alpha_wrapped + pi_val)
        # ... (rest of wrapped logic, ensure use of pi_half, pi_val) ...
        else # Fallback
            zeta2_positive_beta = zeta_perfect # Fallback, types must match or promote
        end
    end

    zeta2 = (beta < zero(beta)) ? -zeta2_positive_beta : zeta2_positive_beta

    zeta_final::promote_type(ZPT, AT, BT, Float64) = zero(promote_type(ZPT,AT,BT,Float64))

    if abs(beta) >= beta0 # beta0 is Float64, abs(beta) will be BT or Dual
        zeta_final = zeta_perfect
    else
        psi1 = zeta_perfect
        psi2 = zeta2

        if beta >= zero(beta)
            if abs(psi1 - pi_half) <= abs(psi2 - pi_half)
                zeta_final = psi1
            else
                zeta_final = psi2
            end
        else
             if abs(psi1 + pi_half) <= abs(psi2 + pi_half)
                 zeta_final = psi1
             else
                 zeta_final = psi2
             end
        end
    end
    return zeta_final
end

"""
    YawThrust(x_eq::Vector{Float64}, t::Float64, NuES0::Float64, Mu::Float64;
              use_rate_limit::Bool = false,
              yaw_rate_max::Float64 = Inf, # rad/s
              orbital_period::Float64 = Inf # seconds
             ) -> Tuple{Vector{Float64}, Float64}

Computes the unit thrust direction vector in RSW coordinates based on yaw steering logic.

Can compute either the "perfect" yaw steering or a rate-limited version.

# Arguments
- `x_eq::Vector{Float64}`: Equinoctial state vector of the spacecraft.
- `t::Float64`: Time (seconds).
- `NuES0::Float64`: Initial Earth true anomaly (radians).
- `Mu::Float64`: Gravitational parameter.

# Keyword Arguments
- `use_rate_limit::Bool`: If true, applies rate limiting. Defaults to false.
- `yaw_rate_max::Float64`: Maximum allowed yaw rate (rad/s). Required if `use_rate_limit` is true.
- `orbital_period::Float64`: Orbital period (seconds). Required if `use_rate_limit` is true.

# Returns
- `thrust_dir_rsw::Vector{Float64}`: Unit thrust vector in RSW coordinates `[u_r, u_s, u_w]`.
- `yaw_angle_used::Float64`: The actual yaw angle (Zeta/zeta) used (either perfect or limited) (radians).

# Throws
- Error if rate limiting is requested without valid `yaw_rate_max` or `orbital_period`.
- See `AlphaCalc`, `BetaCalc`, `apply_yaw_rate_limit`.
"""
# In YawSteering.jl
function YawThrust(x_eq::AbstractVector{EqT}, t::TimeT, NuES0::N0T, Mu::MuParamT;
                  use_rate_limit::Bool = false,
                  yaw_rate_max = Inf,
                  orbital_period = Inf
                 ) where {EqT<:Real, TimeT<:Real, N0T<:Real, MuParamT<:Real}

    local alpha, beta, zeta_perfect, zeta_final, thrust_dir_rsw

    try
        alpha = AlphaCalc(x_eq, t, NuES0, Mu)
        beta = BetaCalc(x_eq, t, NuES0, Mu)
        zeta_perfect = calculate_perfect_yaw(alpha, beta)
        # @show alpha*180/pi
        # @show beta*180/pi
        # @show t
        # sleep(1)
        # @show zeta_perfect*180/pi
        zeta_final = zeta_perfect

        if use_rate_limit
            omega = 2.0 * pi / orbital_period
            if abs(omega) < TOL_ZERO; error("Omega near zero"); end
            zeta_rate_yaw = yaw_rate_max / omega
            local beta0_val::Float64 # beta0 is always Float64
            if abs(zeta_rate_yaw) < TOL_ZERO
                 beta0_val = pi/2
            else
                 beta0_val = abs(atan(1.0 / zeta_rate_yaw))
            end

            zeta_final = apply_yaw_rate_limit(zeta_perfect, alpha, beta, beta0_val, zeta_rate_yaw)

        end

        cosZ = cos(zeta_final)
        sinZ = sin(zeta_final)

        #thrust_dir_rsw = SA[zero(cosZ), -cosZ, -sinZ] 
        thrust_dir_rsw = SA[zero(cosZ), cosZ, -sinZ] #first
        #@show SA[zero(cosZ), cosZ, -sinZ]
        #@show thrust_dir_rsw
        #@show LVLHCorrection(SA[zero(cosZ), cosZ, sinZ])
        # sleep(1)
        
    catch e
        if e isa DimensionMismatch
            println("!!! DimensionMismatch caught INSIDE YawThrust !!!")
            println("Error: ", e)
            # Potentially print variables around the point of SVector construction
            # This might require more granular try-catch blocks if the error isn't at the final SA[]
        end
        println("Error during YawThrust execution: ", e)
        rethrow(e) # Re-throw to see the original stack trace from G
    end

    return thrust_dir_rsw, zeta_final
end