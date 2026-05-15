# File: YawSteering.jl (or MinTimeDynamics.jl if preferred)
# Purpose: Defines dynamics, constraints, sensitivities, and Hamiltonian
#          for optimal control problems, often involving yaw steering.
# Refactored for StaticArrays and generic types (T<:Real) for AD compatibility.

using StaticArrays
using LinearAlgebra
using ForwardDiff

# Assumes Conversions.jl is included and contains GENERIC versions of:
# equin2cart, safe_normalize, Flow, f0, Phi0, B1B3, B2EU, E2BU, etc.
include("Conversions.jl") # Make sure this is the generic version
include("YawSteering.jl") # Make sure this is the generic version

# Define default constants if needed elsewhere, but prefer passing via params
const EARTH_RADIUS_DEFAULT = 6378.1e3 # meters (Float64 ok)
const J2_VALUE_DEFAULT = 1.0826e-3    # Dimensionless (Float64 ok)
const TOL_ZERO_NORM = 1e-12           # Float64 ok for runtime checks

# ==========================================
# Helper Functions (Generic Adjustments)
# ==========================================

"""
    is_thrust_forbidden(t::Float64, forbidden_intervals::Vector{Tuple{Float64, Float64}}) -> Bool
    (Operates on time only, Float64 is usually sufficient)
"""
function is_thrust_forbidden(t::Float64, forbidden_intervals::Vector{Tuple{Float64, Float64}})
    if isempty(forbidden_intervals); return false; end
    # Basic interval check
    for interval in forbidden_intervals
        # Ensure interval bounds are ordered correctly if needed
        # lower = min(interval[1], interval[2])
        # upper = max(interval[1], interval[2])
        # if lower <= t <= upper; return true; end
        if interval[1] <= t <= interval[2]; return true; end # Assuming ordered
    end
    return false
end

"""
    is_eclipsed(r_eci::SVector{3, T}, t::Float64, NuES0::T, Re::T) where {T<:Real} -> Bool

Checks if the position `r_eci` (type T) is in Earth's cylindrical shadow.
Requires `rESECI` to return SVector{3, T}.
Parameters NuES0, Re are type T. Time t is Float64.
"""
function is_eclipsed(r_eci, t::Tt, NuES0::T, Re::T) where {T<:Real, Tt<:Real}
    zero_T = zero(T)
    if norm(r_eci) < Re; return true; end # Check if already inside Earth

    # Assume rESECI is generic and returns SVector{3, T}
    local rES_unit::SVector{3, T}
    try
        rES_unit = SVector{3, T}(rESECI(t, NuES0)) # Ensure SVector{3, T}
    catch e
        error("is_eclipsed requires a generic function 'rESECI(t, NuES0)' returning SVector{3, T}. Error: $e")
    end

    # Check if satellite is on the night side
    sun_angle_check = dot(r_eci, rES_unit) < zero_T
    if !sun_angle_check; return false; end # On sunlit side

    # Check if within the cylindrical shadow radius
    # dist_to_sun_line_sq = norm(r_eci - dot(r_eci, rES_unit) * rES_unit)^2 # More direct projection
    dist_to_sun_line_sq = norm(cross(r_eci, rES_unit))^2 # Equivalent for unit rES_unit
    shadow_check = dist_to_sun_line_sq <= Re^2

    return shadow_check
end


# ==========================================
# Core Dynamics & GVEs (LOCAL COPIES - Made Generic)
# ==========================================

"""
    gve(equin::AbstractVector{T}, mu::T; L::T = kepler(equin[6], equin[2], equin[3])) where {T<:Real}
        -> Tuple{Matrix{T}, T, T}

Computes Gaussian Variation Equations matrix G (6x3), mean motion n, and true longitude L.
Made generic for type T.
Requires local generic `kepler`.
"""
      
# In MinTimeDynamics.jl (inside gve function)
# Original: function gve(equin::AbstractVector{T}, mu::T; L::T = kepler(equin[6], equin[2], equin[3])) where {T<:Real}
# Change to:
function gve(equin::AbstractVector{EqT_gve}, mu_param::MuParamT_gve; L_arg = nothing) where {EqT_gve<:Real, MuParamT_gve<:Real}
    # Determine promoted type for internal calculations
    PromotedT = promote_type(EqT_gve, MuParamT_gve)

    if length(equin) != 6; error("gve: Input must have 6 elements."); end
    a  = equin[1]; ey = equin[2]; ex = equin[3]; qy = equin[4]; qx = equin[5]; l_mean = equin[6]

    # Use PromotedT for constants and intermediate variables where types might mix
    _zero = zero(PromotedT); _one = one(PromotedT); _two = convert(PromotedT, 2.0); 
    _three = convert(PromotedT, 3.0); _half = convert(PromotedT, 0.5)

    # Resolve L with the correct types
    # kepler(l::T, ey::T, ex::T) where T<:Real.
    # If ey, ex are EqT_gve (could be Dual) and l_mean is EqT_gve, then T in kepler becomes EqT_gve.
    local L_true::PromotedT
    if L_arg === nothing
        # Ensure kepler's inputs are consistent or kepler is fully generic itself.
        # If kepler expects all inputs to be same type T, and ey,ex,l_mean are EqT_gve, this is fine.
        L_true = kepler(l_mean, ey, ex) # kepler returns EqT_gve
    else
        L_true = convert(PromotedT, L_arg)
    end
    
    # Check 'a' using its own type or a safely converted zero
    if ForwardDiff.value(a) <= zero(ForwardDiff.value(a))
         error("gve: Semi-major axis 'a' must be positive. Got a = $a")
    end

    e2     = ex^2 + ey^2; # Type EqT_gve (or Dual of)
    if ForwardDiff.value(e2) >= one(ForwardDiff.value(e2)) # Compare value part
        error("gve: Eccentricity squared >= 1. e2 = $e2");
    end
    
    # All these will use promotion:
    n_motion      = sqrt(mu_param / a^_three);
    eta    = sqrt(_one - e2);
    h_ang_mom      = sqrt(mu_param * a) * eta;
    sL     = sin(L_true);
    cL     = cos(L_true);

    por    = _one + ey * sL + ex * cL;
    if abs(ForwardDiff.value(por)) <= eps(typeof(ForwardDiff.value(por)))
         error("gve: Denominator 'por' is near zero. por = $por");
    end
    roh    = h_ang_mom / mu_param / por;

    # GVE components - ensure SVector elements are of PromotedT
    dadt_rsw   = SVector{3, PromotedT}(ex * sL - ey * cL, por, _zero)
    dadt   = _two * a^2 / h_ang_mom * dadt_rsw

    deydt_rsw  = SVector{3, PromotedT}(- por * cL,
                                 ey + (_one + por) * sL,
                                 - ex * (qy * cL - qx * sL))
    deydt = roh * deydt_rsw

    dexdt_rsw  = SVector{3, PromotedT}(por * sL,
                                 ex + (_one + por) * cL,
                                 ey * (qy * cL - qx * sL))
    dexdt = roh * dexdt_rsw
    
    q_term_factor = (_one + qy^2 + qx^2) * _half
    dqydt_rsw  = SVector{3, PromotedT}(_zero, _zero, q_term_factor * sL)
    dqydt = roh * dqydt_rsw

    dqxdt_rsw  = SVector{3, PromotedT}(_zero, _zero, q_term_factor * cL)
    dqxdt = roh * dqxdt_rsw

    eta_p1 = _one + eta
    if abs(ForwardDiff.value(eta_p1)) <= eps(typeof(ForwardDiff.value(eta_p1)))
        error("gve: Denominator '1+eta' is near zero. eta_p1 = $eta_p1");
    end
    term_l1_factor = (por * (ey * sL + ex * cL) / eta_p1 + _two * eta)
    term_l2_factor = (_one + por) * (ey * cL - ex * sL) / eta_p1
    term_l3_factor = (qy * cL - qx * sL)
    dldt_rsw   = SVector{3, PromotedT}(-term_l1_factor, -term_l2_factor, -term_l3_factor)
    dldt = roh * dldt_rsw

    # Construct Matrix{PromotedT}
    G_mat = hcat(dadt, deydt, dexdt, dqydt, dqxdt, dldt)'

    return G_mat, n_motion, L_true;
end

    

"""
    kepler(l::T, ey::T, ex::T; Niter::Int = 20, tol::T = T(1e-12)) where {T<:Real} -> T

Solves Kepler's equation for true longitude `L`. Accepts generic Real type T.
(Local generic version - consider using Conversions.jl version)
"""
function kepler(l::T, ey::T, ex::T; Niter::Int = 20, tol::T = T(1e-12)) where {T<:Real}
    # This is the same generic implementation as provided for Conversions.jl
    e_sq = ex^2 + ey^2
    one_T = one(T); zero_T = zero(T)
    if e_sq >= one_T; error("Kepler solver: Eccentricity squared >= 1."); end
    eta = sqrt(one_T - e_sq)
    eta_p1 = one_T + eta
    if eta_p1 <= eps(T); error("Kepler solver: Denominator 1+eta near zero."); end
    AA = one_T - ex^2 / eta_p1; BB = ey * ex / eta_p1; CC = one_T - ey^2 / eta_p1
    K = l; f = zero_T; df = one_T; converged = false
    for j in 1:Niter
        cK = cos(K); sK = sin(K)
        f = K + ey * cK - ex * sK - l
        df = one_T - ey * sK - ex * cK
        if abs(f) < tol; converged = true; break; end
        if abs(df) < eps(T); error("Kepler solver: Derivative df near zero."); end
        K -= f / df
    end
    if !converged; error("Kepler solver: Did not converge."); end
    cK = cos(K); sK = sin(K)
    sL_num = (AA * sK + BB * cK - ey); cL_num = (BB * sK + CC * cK - ex)
    L = atan(sL_num, cL_num)
    return L # Returns type T
end

"""
    rv(EQ::SVector{6, T}, Mu::T) where {T<:Real} -> SVector{6, T}

Computes Cartesian state [r; v] (SVector) from equinoctial state (SVector).
(Local version - consider using Conversions.jl version if it exists)
Requires generic `equin2cart`.
"""
function rv(EQ::SVector{6, T}, Mu::T) where {T<:Real}
   # Assumes equin2cart from Conversions.jl is generic and returns SVector{6, T}
   cart_vec = equin2cart(EQ, Mu)
   return cart_vec
end

"""
    rNoT(EQ::AbstractVector{EQT}, Mu::MuT) where {EQT<:Real, MuT<:Real} -> SVector{3, EQT} or SVector{3, Dual}

Computes position vector r (SVector) from equinoctial state.
MUST be generic and handle Dual numbers if EQT is a Dual type.
Requires a generic `equin2cart` function.
"""
function rNoT(EQ::AbstractVector{EQT}, Mu::MuT) where {EQT<:Real, MuT<:Real}
    # Assuming equin2cart works generically and returns something indexable
    # with element type matching its input (e.g., EQT if Mu is not Dual, Dual if EQT is Dual)
    cart_vec = equin2cart(EQ, Mu) # This MUST handle Dual numbers in EQ and Mu

    # Extract elements. The resulting SVector's element type will depend on cart_vec's.
    # If cart_vec has Duals, this SVector will have Duals.
    # Use SA[] constructor for StaticArrays
    # Return type is inferred based on cart_vec's element type.
    return SA[cart_vec[1], cart_vec[2], cart_vec[3]]
end


"""
    ∂r∂yNoT(Mu::T, yp0::SVector{6, T}) where {T<:Real} -> SMatrix{3, 6, T, 18}

Computes Jacobian ∂r/∂y using ForwardDiff.
Requires generic `rNoT` and `equin2cart`.
"""
function ∂r∂yNoT(Mu::T, yp0::SVector{6, T}) where {T<:Real}
    # Define the closure WITHOUT the input type constraint on `y`.
    # Let ForwardDiff pass SVector{6, Dual{...}} into it.
    # The captured Mu::T is fine.
    r_func = y -> rNoT(y, Mu)

    # Calculate the Jacobian. ForwardDiff uses Dual numbers internally.
    # The resulting Jacobian matrix J will have element type T (e.g., Float64).
    J = ForwardDiff.jacobian(r_func, yp0)

    # Assert the final type and size AFTER the calculation.
    return J::SMatrix{3, 6, T, 18}
end

# In MinTimeDynamics.jl (around line 222)

"""
    Flow(x::SVector{6, XT}, t::TimeT, Mu::MuT) where {XT<:Real, TimeT<:Real, MuT<:Real} -> SVector{6, <PromotedType>}

Computes Keplerian flow (updates mean longitude).
Allows x, t, and Mu to be different Real subtypes (including Duals).
"""
function Flow(x::AbstractVector, t::Real, Mu::MuT) where {MuT<:Real}
    a = x[1] # Type XT or Dual{...,XT,...}

    # Determine a common type for calculations if needed, or let promotion work.
    # For constants, use the most general type or promote explicitly.
    # Example: zero_promoted = zero(promote_type(XT, TimeT, MuT))
    # However, often direct operations work fine due to Julia's promotion.

    # Ensure constants like 3.0, 2.0, pi are converted to a suitable type if they interact
    # with Dual numbers. Using `oneunit` or direct conversion can help.
    # For example, if `a` is Dual, `a^3` is fine.
    # If `Mu` is Float64 and `a` is Dual, `Mu / a^3` will promote correctly.

    _zero = zero(a) # Use type of 'a' for zero, or promote further if needed
    _three = convert(typeof(a), 3.0) # Convert 3.0 to type of 'a' (or a promoted type)
    _two = convert(typeof(t), 2.0)   # Convert 2.0 to type of 't'
    _pi = convert(typeof(t), pi)     # Convert pi to type of 't'


    if a <= _zero; error("Flow: Semi-major axis 'a' must be positive."); end
    
    # term_under_sqrt will be of a promoted type (XT, MuT)
    term_under_sqrt = Mu / a^_three # Mu is MuT, a is XT (or Dual of XT)
    
    # Check for negative before sqrt if types can be Real (not just Dual of positive)
    # If XT can be a non-Dual Real, this check is important.
    # For Duals, if the value part is negative, sqrt will error appropriately.
    if !(typeof(term_under_sqrt) <: ForwardDiff.Dual) && term_under_sqrt < zero(term_under_sqrt) # Check only if not Dual
        error("Flow: Mu/a^3 is negative ($term_under_sqrt). Mu=$Mu, a=$a")
    end
    
    n = sqrt(term_under_sqrt) # n will be of promoted type, handles Duals

    # Update mean longitude - t is TimeT (or Dual of TimeT)
    # n*t will promote correctly.
    #l_new = mod(x[6] + n * t, _two * _pi)
    l_new = x[6] + n * t

    # The result SVector's element type will be promoted from XT, TimeT, MuT
    return SVector(x[1], x[2], x[3], x[4], x[5], l_new)
end

"""
    f0(x::SVector{6, T}, Mu::T) where {T<:Real} -> SVector{6, T}

Keplerian drift rate vector [0,0,0,0,0,n].
(Local generic version - consider using Conversions.jl version)
"""
function f0(x, Mu::T) where {T<:Real}
    a = x[1];
    zero_T = zero(T); three_T = T(3.0)
    if a <= zero_T; error("f0: Semi-major axis 'a' must be positive."); end
    term_under_sqrt = Mu / a^three_T
    if term_under_sqrt < zero_T; error("f0: Mu/a^3 is negative."); end
    n = sqrt(term_under_sqrt)
    return SA[zero_T, zero_T, zero_T, zero_T, zero_T, n] # StaticArray literal
end
"""
    Phi0(x::SVector{6, T}, t::T, Mu::T) where {T<:Real} -> SMatrix{6, 6, T, 36}

State transition matrix for Keplerian drift (∂Flow(x,t)/∂x).
TIME 't' IS NOW TYPE T.
(Local generic version - consider using Conversions.jl version)
"""
function Phi0(x::SVector{6, XT}, t_arg::TimeT, mu_arg::MuParamT) where {XT<:Real, TimeT<:Real, MuParamT<:Real}
    a = x[1]
    #zero_T = zero(T); five_T = T(5.0)
    PromotedElementType = promote_type(XT, TimeT, MuParamT)
    _zero = zero(PromotedElementType)
    _one = one(PromotedElementType)
    _five = convert(PromotedElementType, 5.0)
    _minus_one_point_five = convert(PromotedElementType, -1.5)
    # Check 'a' using its own type or a safely converted zero
    if ForwardDiff.value(a) <= zero(ForwardDiff.value(a)) # Compare value part for safety with Duals
         error("Phi0: Semi-major axis 'a' must be positive. Got a = $a")
    end
    term_under_sqrt = mu_arg / (a^_five)
    val_term_under_sqrt = (term_under_sqrt isa ForwardDiff.Dual) ? ForwardDiff.value(term_under_sqrt) : term_under_sqrt
    if val_term_under_sqrt < zero(val_term_under_sqrt)
        error("Phi0: Argument to sqrt (mu/a^5) is negative: $term_under_sqrt. Mu=$mu_arg, a=$a")
    end
    sqrt_term = sqrt(term_under_sqrt) # Promoted type
    val = _minus_one_point_five * sqrt_term * t_arg
    mat = SMatrix{6, 6, PromotedElementType, 36}(
        _one,   _zero,  _zero,  _zero,  _zero,  val,
        _zero,  _one,   _zero,  _zero,  _zero,  _zero,
        _zero,  _zero,  _one,   _zero,  _zero,  _zero,
        _zero,  _zero,  _zero,  _one,   _zero,  _zero,
        _zero,  _zero,  _zero,  _zero,  _one,   _zero,
        _zero,    _zero,  _zero,  _zero,  _zero,  _one
    )
    return mat
end

"""
    G(x_gve::AbstractVector, mu_gve, t_g, nu_es0_gve, odb_gve::Int)

Computes the G matrix for the Gaussian Variational Equations.
This version is generic and compatible with ForwardDiff.
It relies on `gve` and `YawThrust` also being AD-compatible.
"""
function G(x_gve::AbstractVector, mu_gve, t_g, nu_es0_gve, odb_gve::Int)
    
    # 1. Call the (assumed to be generic) gve function.
    # The result `GMat_raw_dense_tuple` will contain Duals if x_gve does.
    GMat_raw_dense_tuple = gve(x_gve, mu_gve)
    GMat_raw_dense = GMat_raw_dense_tuple[1]

    # 2. Use the generic SMatrix constructor. It will automatically infer the
    #    correct element type (e.g., Float64 or Dual{...}).
    GMat_raw = SMatrix{6, 3}(GMat_raw_dense)

    if odb_gve == 1 # BangBang/Yaw control magnitude only
        
        # 3. Call the (assumed to be generic) YawThrust function.
        # `thrust_vector_part` will contain Duals if x_gve or t_g does.
        thrust_vector_part, _ = YawThrust(x_gve, t_g, nu_es0_gve, mu_gve)
        
        # 4. Convert to SVector using the generic constructor.
        # This ensures thrust_dir_rsw is a static vector with the correct element type.
        thrust_dir_rsw = SVector{3}(thrust_vector_part)

        # The matrix-vector product is AD-compatible out of the box.
        GMat_vec = GMat_raw * thrust_dir_rsw
        return GMat_vec

    else # odb_gve == 0
        return GMat_raw
    end
end

"""
    DMatL(tdif::Float64, Mu::T, yp0::SVector{6, T}, t::Float64, NuES0::T, ODB::Int) where {T<:Real}
        -> Union{SVector{6, T}, SMatrix{6, 3, T, 18}}

Sensitivity mapping D = Φ(t, t-tdif) * G(t-tdif).
Requires generic `Flow`, `Phi0`, `G`.
"""
function DMatL(tdif::Float64, Mu::T, yp0::SVector{6, T}, t::Float64, NuES0::T, ODB::Int) where {T<:Real}
    # State at time t - tdif (relative to reference state yp0 at t=0)
    FlowVal = Flow(yp0, - tdif, Mu); # State at time control is applied
    # STM from (t-tdif) to t
    PhiVal = Phi0(FlowVal, tdif, Mu) # Propagates effect over duration tdif
    # G evaluated at state FlowVal at time (t)
    G_at_t = G(FlowVal, Mu, t, NuES0, ODB)
    
    # Result is 6x1 SVector or 6x3 SMatrix
    DMatL_val = PhiVal * G_at_t
    return DMatL_val
end

"""
    U(Mu::T, yp0::SVector{6, T}, Δvtca::SVector{3, T},
      OTE::SMatrix{3, 6, T, 18}, t::Float64, NuES0::T, ODB::Int) where {T<:Real}
        -> Union{SVector{3, T}, SMatrix{3, 3, T, 9}}

Computes U matrix/vector: Maps control u(t) to B-plane deviation rate ∂(δr_B)/∂u.
Requires generic `safe_normalize`, `DMatL`.
"""
function U(Mu::T, yp0::SVector{6, T}, Δvtca::SVector{3, T},
           OTE::SMatrix{3, 6, T, 18}, t::Float64, NuES0::T, ODB::Int) where {T<:Real}

    # Assumes safe_normalize from Conversions.jl is generic
    B3U = safe_normalize(Δvtca, "Δvtca in U")
    # Build projection matrix using StaticArrays operations (type T)
    Proj_Bplane = one(SMatrix{3, 3, T}) - B3U * transpose(B3U)

    # Calculate DMatL (Result is 6x1 SVector or 6x3 SMatrix)
    # Control applied at time t, effect measured at t=0 (TCA), so tdif = 0 - t = -t
    DMatL_val = DMatL(-t, Mu, yp0, t, NuES0, ODB)

    # Calculate UE = Proj_Bplane * OTE * DMatL_val
    # Result UE is 3x1 SVector or 3x3 SMatrix
    UE = Proj_Bplane * OTE * DMatL_val

    # Check for near-zero result (use eps(T))
    if norm(UE) < eps(T) # Check norm for vector or matrix
         # @warn "U matrix/vector is near zero at time t=$t. Norm: $(norm(UE)). Adding small perturbation."
         # Use T literals for perturbation
         one_T = one(T)
         perturb_val = sqrt(eps(T))
         if ODB == 1 # UE is SVector{3}
             UE = UE .+ perturb_val * SA[one_T, one_T, one_T]
         else # UE is SMatrix{3,3}
             UE = UE .+ perturb_val * one(SMatrix{3, 3, T})
         end
     end

    return UE
end


# ==========================================
# ODE Function for Backward Integration (Generic)
# ==========================================
"""
    IntδΔr!(dX::Vector{T}, X::Vector{T}, p::NamedTuple, t::Float64) where {T<:Real}

ODE right-hand side for backward integration of δΔr.
State X is Vector{T}, parameters p contain generic types T. Time t is Float64.
Requires generic `Flow`, `equin2cart`, `is_eclipsed`, `is_thrust_forbidden`, `U`.
"""
function IntδΔr!(dX::Vector{T}, X::Vector{T}, p::NamedTuple, t::Float64) where {T<:Real}

    # Unpack parameters by name from NamedTuple 'p' (expecting type T where appropriate)
    Mu         = p.Mu::T
    yp0        = p.yp0::SVector{6, T}

    Δvtca      = p.Δvtca::SVector{3, T}
    NTheta     = p.NTheta::SVector{3, T} # ECI target direction on B-plane unit sphere

    OTE        = p.OTE::SMatrix{3, 6, T, 18}
    NuES0      = p.NuES0::T
    ODB        = p.ODB::Int
    ECB        = p.EclipseBool::Int
    LXB        = p.LuxBool::Int
    TuningP    = p.TuningP::T # Assume TuningP might be differentiated w.r.t.
    Re_const   = p.Re_const::T
    forbidden_intervals = p.forbidden_intervals::Vector{Tuple{Float64, Float64}}
   zero_T = zero(T); one_T = one(T)

    # Current state deviation (X is the state vector [δΔr_x, δΔr_y, δΔr_z])
    # No need to convert X if it's just the 3-element state.

    # --- Get current ECI position for eclipse check ---
    # Uses nominal Keplerian path from yp0 at time t
    FlowVal = Flow(yp0, t, Mu) # SVector{6, T}
    # Assumes equin2cart is generic
    r_eci_t = equin2cart(FlowVal, Mu)[SA[1, 2, 3]] # SVector{3, T}
    # --- Check Constraints ---
    # Assumes is_eclipsed and is_thrust_forbidden are generic/correct
    AlwaysOn = false
    if ECB == 1
        is_ecl = is_eclipsed(r_eci_t, t, NuES0, Re_const)
        AlwaysOn = false
    else
        is_ecl = false # No eclipse check if ECB == 0
        AlwaysOn = false
    end
    is_forbid = is_thrust_forbidden(t, forbidden_intervals)

    # --- Calculate Dynamics Matrix U ---
    # Assumes U is generic
    UE = U(Mu, yp0, Δvtca, OTE, t, NuES0, ODB) # Returns SVector{3, T} or SMatrix{3, 3, T, 9}

    # --- Determine Optimal Control u ---
    uOpt_val = zero_T # Default for ODB=1 scalar
    uOpt_vec = SA[zero_T, zero_T, zero_T] # Default for ODB=0 vector

    if !(is_ecl || is_forbid) # Only calculate non-zero control if thrust allowed
        if ODB == 1 # u is scalar magnitude [0, 1]
            if LXB == 1  # Luxurious case, two thrusters LXB == 1 || 
                InnerProduct = dot(UE, NTheta) # UE::SVector{3,T}, NTheta::SVector{3,T} -> T
                # Simple bang-bang (use eps(T) for comparison)
                if InnerProduct > eps(T)
                    uOpt_val = one_T # Apply thrust
                else
                    #uOpt_val = zero_T # Or small value if needed
                    uOpt_val = -one_T # Small positive value
                    #uOpt_val = one_T
                end
            else 
                InnerProduct = dot(UE, NTheta) # UE::SVector{3,T}, NTheta::SVector{3,T} -> T

                # Simple bang-bang (use eps(T) for comparison)
                    
                if InnerProduct > eps(T)
                    uOpt_val = one_T # Apply thrust 
                else# Current state deviation
                    if AlwaysOn 
                        uOpt_val = one_T
                    # uOpt_val = zero_T # Or small value if needed
                    else
                        uOpt_val = -one_T#T(1e-8) # Small positive value
                    #uOpt_val = one_T
                    end
                end
            end
            # Apply smoothing if needed (using TuningP)
            # uOpt_val = one_T / (one_T + exp(-InnerProduct / TuningP))
        else # ODB == 0, u is 3x1 vector in RSW (U maps RSW control to ECI B-plane deviation)
            # We want control u that maximizes <NTheta, U*u> subject to ||u||=1
            # Optimal u direction is proportional to U' * NTheta
            uOpt_dir = transpose(UE) * NTheta # UE::SMatrix{3,3,T}, NTheta::SVector{3,T} -> SVector{3,T}
            norm_uOpt_dir = norm(uOpt_dir)
            if norm_uOpt_dir > eps(T) # Use eps(T)
                uOpt_vec = uOpt_dir / norm_uOpt_dir # Normalized optimal direction
            end
            # else uOpt_vec remains zero
        end
    end

    # --- Calculate State Derivative dδΔr/dt ---
    local dδΔr::SVector{3, T} # Declare type
    if ODB == 1
        dδΔr = UE * uOpt_val # SVector{3,T} * T -> SVector{3,T}
    else
        dδΔr = UE * uOpt_vec # SMatrix{3,3,T} * SVector{3,T} -> SVector{3,T}
    end

    # --- Update the standard Vector dX used by the solver ---
    dX .= dδΔr # Assign components from SVector{3, T} to Vector{T}
end

# ==========================================
# Function to Get Control Input (Generic)
# ==========================================
"""
    GetCon(p::NamedTuple, t::Float64) where {T<:Real}
        -> Tuple{SVector{3, T}, T}

Calculates the optimal control vector `u_opt` (RSW frame, type T) and the
inner product `<UE, NTheta>` (type T) at a given time `t` (Float64).
Parameters `p` contain generic types T.
Requires generic `Flow`, `equin2cart`, `is_eclipsed`, `is_thrust_forbidden`, `U`, `YawThrust`.
"""
function GetCon(p::NamedTuple, t::Float64)  # Infer T from p
    # Unpack parameters by name (expecting type T where appropriate)
    _T         = eltype(p.yp0) # Infer T inside if needed
    Mu         = p.Mu::_T
    yp0        = p.yp0::SVector{6, _T}
    Δvtca      = p.Δvtca::SVector{3, _T}
    NTheta     = p.NTheta::SVector{3, _T}
    OTE        = p.OTE::SMatrix{3, 6, _T, 18}
    NuES0      = p.NuES0::_T
    ODB        = p.ODB::Int
    ECB        = p.EclipseBool::Int
    LXB        = p.LuxBool::Int
    TuningP    = p.TuningP::_T
    Re_const   = p.Re_const::_T
    forbidden_intervals = p.forbidden_intervals::Vector{Tuple{Float64, Float64}}

    zero_T = zero(_T); one_T = one(_T)
    nan_T = _T(NaN)

    # --- Check constraints ---
    FlowVal = Flow(yp0, t, Mu); # SVector{6, T}
    r_eci_t = equin2cart(FlowVal, Mu)[SA[1, 2, 3]] # SVector{3, T}
    if ECB == 1
        is_ecl = is_eclipsed(r_eci_t, t, NuES0, Re_const)
    else
        is_ecl = false # No eclipse check if ECB == 0
    end
    is_forbid = is_thrust_forbidden(t, forbidden_intervals)

    # --- Calculate U ---
    UE = U(Mu, yp0, Δvtca, OTE, t, NuES0, ODB) # SVector{3, T} or SMatrix{3, 3, T, 9}

    # --- Calculate Control ---
    InnerProduct = nan_T # Default value
    accDir_rsw = SA[zero_T, zero_T, zero_T] # Default zero acceleration in RSW
    

    local thrust_dir_rsw::SVector{3, _T}
    if !(is_ecl || is_forbid)
        if ODB == 1
            if LXB == 1 # Luxurious case, two thrusters
                InnerProduct = dot(UE, NTheta) # T
                # uMag = one_T / (one_T + exp(-InnerProduct / TuningP)) # Smoothing
                uMag = (InnerProduct > eps(_T)) ? one_T : -1*one_T#_T(1e-8) # Bang-bang
                
                try
                    thrust_dir_rsw_tuple_part, _ = YawThrust(FlowVal, t, NuES0, Mu)
                    thrust_dir_rsw = thrust_dir_rsw_tuple_part # This is SVector
                catch e
                    error("GetCon requires a generic function 'YawThrust(x, t, NuES0, Mu)' returning SVector{3, T}. Error: $e")
                end
                accDir_rsw = thrust_dir_rsw * uMag # Applied thrust vector in RSW
            else
                InnerProduct = dot(UE, NTheta) # T
                # uMag = one_T / (one_T + exp(-InnerProduct / TuningP)) # Smoothing
                uMag = (InnerProduct > eps(_T)) ? one_T : _T(1e-8) # Bang-bang
                try
                    thrust_dir_rsw_tuple_part, _ = YawThrust(FlowVal, t, NuES0, Mu)
                    thrust_dir_rsw = thrust_dir_rsw_tuple_part # This is SVector
                catch e
                    error("GetCon requires a generic function 'YawThrust(x, t, NuES0, Mu)' returning SVector{3, T}. Error: $e")
                end
                accDir_rsw = thrust_dir_rsw * uMag # Applied thrust vector in RSW
            end
        else # ODB == 0
            uOpt_dir = transpose(UE) * NTheta # SVector{3, T}
            norm_uOpt_dir = norm(uOpt_dir)
            if norm_uOpt_dir > eps(_T) # Use eps(T)
                accDir_rsw = uOpt_dir / norm_uOpt_dir # Unit thrust vector in RSW
            end
            # InnerProduct could be norm or zero
            InnerProduct = norm_uOpt_dir # Return the norm for ODB=0 case
        end
    end

    # Return applied thrust vector (SVector{3, T} in RSW) and inner product (T)
    return accDir_rsw, InnerProduct
end

"""
    GetConMag(p::NamedTuple, t::Float64) where {T<:Real} -> T

Calculates the optimal control magnitude (smoothed bang-bang or similar) at time `t`.
Refactored to use NamedTuple `p` containing generic types T.
Requires generic `U`.
"""
function GetConMag(p::NamedTuple, t::Float64) # Infer T from p
    _T         = eltype(p.yp0) # Infer T inside if needed
    # Unpack parameters
    Mu         = p.Mu::_T
    yp0        = p.yp0::SVector{6, _T}
    Δvtca      = p.Δvtca::SVector{3, _T}
    NTheta     = p.NTheta::SVector{3, _T}
    OTE        = p.OTE::SMatrix{3, 6, _T, 18}
    NuES0      = p.NuES0::_T
    ODB        = p.ODB::Int
    ECB        = p.EclipseBool::Int
    LXB        = p.LuxBool::Int
    TuningP    = p.TuningP::_T # Assume TuningP is passed in p
    # Constraints needed if control is zero during eclipse/forbidden
    Re_const   = p.Re_const::_T
    forbidden_intervals = p.forbidden_intervals::Vector{Tuple{Float64, Float64}}

    zero_T = zero(_T); one_T = one(_T)

    # --- Check constraints ---
    FlowVal = Flow(yp0, t, Mu); # SVector{6, T}
    r_eci_t = equin2cart(FlowVal, Mu)[SA[1, 2, 3]] # SVector{3, T}
    if ECB == 1
        is_ecl = is_eclipsed(r_eci_t, t, NuES0, Re_const)
    else
        is_ecl = false # No eclipse check if ECB == 0
    end
    is_forbid = is_thrust_forbidden(t, forbidden_intervals)

    uOptN = zero_T # Default control magnitude

    if !(is_ecl || is_forbid)
        # Calculate U (only makes sense for ODB=1 if we want a scalar magnitude)
        if ODB != 1
            @warn "GetConMag is typically used for ODB=1 (scalar magnitude control)."
            # Fallback or error? For now, calculate as if ODB=1
        end
        # Force calculation as if ODB=1 to get UE as SVector
        UE_vec = U(Mu, yp0, Δvtca, OTE, t, NuES0, 1) # Calculate 3x1 UE

        InnerProduct = dot(UE_vec, NTheta) # T

        # Smoothed version (ensure TuningP is positive and non-zero)
        if TuningP <= eps(_T)
             @warn "GetConMag: TuningP is near zero or negative ($TuningP). Using bang-bang."
             uOptN = (InnerProduct > eps(_T)) ? one_T : _T(1e-8) # Bang-bang fallback
        else
             uOptN = one_T / (one_T + exp(-InnerProduct / TuningP))
        end

        # Simple Bang-Bang version:
        # uOptN = (InnerProduct > eps(T)) ? one_T : T(1e-8) # Or zero_T
    end

    return uOptN # Return scalar control magnitude (type T)
end

# ==========================================
# Constraint Dynamics (g, ∂g∂y, ∂g∂tCA) - Generic
# ==========================================

"""
    rvs(t::Float64, p::NamedTuple) where {T<:Real} -> SVector{6, T}

Computes Cartesian state of secondary object at time `t` using Keplerian flow.
Assumes `p` contains `p.Mu::T` and `p.ys0::SVector{6, T}`.
Requires generic `Flow`, `equin2cart`.
"""
function rvs(t::Float64, p::NamedTuple) # Infer T from p
   # T = eltype(p.yp0) # Infer T inside if needed
   # Uses Flow which takes SVector{6,T}, Float64, T; returns SVector{6,T}
   # Uses equin2cart which takes SVector{6,T}, T; returns SVector{6,T}
   ys_at_t = Flow(p.ys0, t, p.Mu) # p.ys0 must be SVector{6, T}
   return equin2cart(ys_at_t, p.Mu) # Returns SVector{6, T}
end

"""
    g(y::SVector{6, YT}, tCA::TCAT, p::NamedTuple) where {YT<:Real, TCAT<:Real} -> <PromotedType>

Computes constraint g = Δr ⋅ Δv = 0.
Allows y's elements (YT) and tCA (TCAT) to be different Real subtypes (including Duals).
`p` must contain `p.Mu` (Real) and `p.ys0` (SVector{6, Real}).
"""
function g(y::AbstractVector, tCA::Real, p::NamedTuple) #where {TCAT<:Real}
    # Evaluate states at tCA using Keplerian flow
    # Flow now accepts time tCA::T directly
    yp_tca = Flow(y, tCA, p.Mu)     
    ys_tca = Flow(p.ys0, tCA, p.Mu) 

    # Get Cartesian states

    cart_p = equin2cart(yp_tca, p.Mu)
    cart_s = equin2cart(ys_tca, p.Mu)

    # Use static indexing
    delta_r = cart_p[SA[1, 2, 3]] - cart_s[SA[1, 2, 3]]
    delta_v = cart_p[SA[4, 5, 6]] - cart_s[SA[4, 5, 6]]

    return dot(delta_r, delta_v)
end
# --- Gradient functions using ForwardDiff (Generic) ---

"""
    ∂g∂y(y::SVector{6, T}, tCA::T, p::NamedTuple) where {T<:Real} -> SVector{6, T}

Computes gradient ∇y g(y, tCA) using ForwardDiff.
Requires generic `g`.
"""
function ∂g∂y(y::AbstractVector, tCA::Real, p::NamedTuple) #where {TCAT<:Real}
    # Closure captures tCA and p
    g_func = y_state -> g(y_state, tCA, p)
    ∇g = ForwardDiff.gradient(g_func, y)
    return ∇g
end


"""
    ∂g∂tCA(y::SVector{6, T}, tCA::T, p::NamedTuple) where {T<:Real} -> T

Computes partial derivative ∂g/∂tCA using ForwardDiff.
Requires generic `g` (which now handles tCA::T).
"""
function ∂g∂tCA(y::AbstractVector, tCA::Real, p::NamedTuple) #where {TCAT<:Real}
    # Closure captures y and p
    # The input tca_val to this closure will be Dual if differentiating wrt tCA
    g_func = tca_val -> g(y, tca_val, p) # Calls corrected g
    # ForwardDiff.derivative passes Dual tca_val to g_func
    ∂g∂t = ForwardDiff.derivative(g_func, tCA) # tCA is the point of evaluation
    return ∂g∂t
end

# ==========================================
# B Matrix (Sensitivity of y, tca to Control u) - Generic
# ==========================================


# In MinTimeDynamics.jl

"""
    DMat(tdif::TDIF_dm, y::SVector{6, YT_dm}, t_current::TCURR_dm, p::NamedTuple) where {TDIF_dm<:Real, YT_dm<:Real, TCURR_dm<:Real}
        -> Union{SVector{6, <PromotedElType>}, SMatrix{6, 3, <PromotedElType>}}

Computes sensitivity D = Φ(t_current, t_current-tdif) * G(t_current-tdif).
Allows tdif, y elements, and t_current to be different Real subtypes.
Requires generic `Flow`, `Phi0`, `G`.
`p` must contain `p.Mu`, `p.NuES0`, `p.ODB`.
"""
function DMat(tdif::TDIF_dm, y::AbstractVector, t_current::TCURR_dm, p::NamedTuple) where {TDIF_dm<:Real, TCURR_dm<:Real}
    # Determine the time at which control is applied.
    # This will be of a promoted type, potentially Dual.
    time_control_applied = t_current - tdif

    # Calculate the state at the time control is applied.
    # Flow(y::SVector{6,YT_dm}, time_control_applied::Promoted(TCURR,TDIF), p.Mu::TypePMu)
    # FlowVal will have elements of type promote_type(YT_dm, TCURR_dm, TDIF_dm, typeof(p.Mu))
    FlowVal = Flow(y, time_control_applied, p.Mu)

    # Calculate the State Transition Matrix for the drift.
    # Phi0(FlowVal, tdif::TDIF_dm, p.Mu::TypePMu)
    # PhiVal will have elements of type promote_type(eltype(FlowVal), TDIF_dm, typeof(p.Mu))
    PhiVal = Phi0(FlowVal, tdif, p.Mu)

    G_at_t = G(FlowVal, p.Mu, time_control_applied, p.NuES0, p.ODB)

    # The result DMat_val will have elements of a type promoted from PhiVal and G_at_t.
    DMat_val = PhiVal * G_at_t
    return DMat_val
end

"""
    By(t::Float64, tCA::T, y::SVector{6, T}, p::NamedTuple) where {T<:Real}
        -> Union{SVector{6, T}, SMatrix{6, 3, T, 18}}

Computes By = (I - (f0 * ∂g∂y') / (∂g∂y' * f0 + ∂g∂tCA)) * DMat sensitivity.
Requires generic `f0`, `∂g∂y`, `∂g∂tCA`, `DMat`. Parameters `p` contain type T.
"""
function By(t::Float64, tCA::Real, y::AbstractVector, p::NamedTuple) # where {TCAT_by<:Real}
    # ... (calculations for f_drift, ∇g_y, PgPyE, ∂g∂t, gaccE are fine) ...
    f_drift = f0(y, p.Mu) # SVector{6, T}
    ∇g_y = ∂g∂y(y, tCA, p) # SVector{6, T}
    PgPyE = transpose(∇g_y) # Transpose{T, SVector{6, T}} (1x6 row vector)
    ∂g∂t = ∂g∂tCA(y, tCA, p) # T
    gaccE = dot(∇g_y, f_drift) + ∂g∂t # T

    val_gaccE = (gaccE isa ForwardDiff.Dual) ? ForwardDiff.value(gaccE) : gaccE
    if abs(val_gaccE) < eps(typeof(val_gaccE)) # Compare value part
        error("By calculation: gaccE is near zero ($gaccE). Inputs: tCA=$tCA, y=$y")
    end

    # Let tdif_val be tCA - t_float. Its type will be promoted from TCAT_by and Float64.
    tdif_val = tCA - t # If TCAT_by is Dual, tdif_val is Dual. If TCAT_by is Float64, tdif_val is Float64.
    
    # y is SVector{6, YT_by}
    # t_float is Float64

    # DMat maps control at time 't' to state change at time 'tCA'
    # tdif = tCA - t (duration from control application to TCA, type T)
    # t_current = t (time control is applied, type Float64)
    DMat_val = DMat(tdif_val, y, t, p)

        # Element type for the intermediate calculations before multiplying by DMat_val
    # This will be a promotion of typeof(gaccE), eltype(f_drift), eltype(PgPyE)
    Term2ElType = promote_type(typeof(gaccE), eltype(f_drift), eltype(PgPyE))
    
    # I6 should be compatible with Term2ElType for subtraction, or use a common promoted type.
    # The most general element type for (I6 - term2) will be Term2ElType.
    I6 = one(SMatrix{6, 6, Term2ElType})

    # Ensure one(gaccE) is used if gaccE can be Dual.
    # If gaccE is scalar, one(gaccE) is its multiplicative identity.
    term2_factor = (one(gaccE) / gaccE)
    term2_matrix_part = f_drift * PgPyE # Resulting matrix eltype also promotes from f_drift and PgPyE

    # term2 will have eltype promoted from term2_factor and term2_matrix_part
    term2 = term2_factor * term2_matrix_part

    # By_val = (I6 - term2) * DMat_val
    # The element type of (I6 - term2) is Term2ElType.
    # The element type of DMat_val is eltype(DMat_val).
    # The final By_val eltype will be promote_type(Term2ElType, eltype(DMat_val)).
    By_val = (I6 - term2) * DMat_val


    
    return By_val
end

"""
    BtCA(t::Float64, tCA::T, y::SVector{6, T}, p::NamedTuple) where {T<:Real}
        -> Union{T, Transpose{T, SVector{3, T}}}

Computes BtCA = (- ∂g∂y' / (∂g∂y' * f0 + ∂g∂tCA)) * DMat sensitivity.
Requires generic `f0`, `∂g∂y`, `∂g∂tCA`, `DMat`. Parameters `p` contain type T.
"""
function BtCA(t_float::Float64, tCA::Real, y::AbstractVector, p::NamedTuple) #where {TCAT_btca<:Real}
    # f_drift will have eltype promoted from YT_btca and typeof(p.Mu)
    f_drift = f0(y, p.Mu)

    # ∇g_y will have eltype promoted from YT_btca, TCAT_btca, and types in p
    ∇g_y = ∂g∂y(y, tCA, p)
    PgPyE = transpose(∇g_y) # PgPyE is a 1x6 row vector (Adjoint or Transpose)

    # ∂g∂t will have type promoted from YT_btca, TCAT_btca, and types in p
    ∂g∂t = ∂g∂tCA(y, tCA, p)

    # gaccE will be a promoted scalar type (potentially Dual)
    gaccE = dot(∇g_y, f_drift) + ∂g∂t

    # Check gaccE's value part for being near zero
    val_gaccE = (gaccE isa ForwardDiff.Dual) ? ForwardDiff.value(gaccE) : gaccE
    if abs(val_gaccE) < eps(typeof(val_gaccE)) # Compare value part against its own type's epsilon
        error("BtCA calculation: gaccE is near zero ($gaccE). Inputs: t_float=$t_float, tCA=$tCA, y=$y")
    end

    # Calculate tdif for DMat:
    # t_float is Float64.
    # tCA is TCAT_btca (can be Dual or Float64).
    # tdif_val will be promoted type (e.g., Dual if tCA is Dual, Float64 otherwise).
    tdif_val = tCA - t_float

    # Call DMat:
    # DMat(tdif::Promoted, y::SVector{6,YT_btca}, t_current::Float64, p)
    # DMat_val can have Dual elements.
    DMat_val = DMat(tdif_val, y, t_float, p)

    # Calculate BtCA_val:
    # -(one(gaccE) / gaccE) is a scalar (potentially Dual).
    # PgPyE is 1x6 (elements potentially Dual).
    # DMat_val is 6x1 (SVector, elements potentially Dual) or 6x3 (SMatrix, elements potentially Dual).
    # The result of PgPyE * DMat_val will be scalar or 1x3.
    # Then multiplied by the scalar factor.

    # one(gaccE) ensures the '1' is of a type compatible with gaccE (e.g., Dual if gaccE is Dual)
    factor = -(one(gaccE) / gaccE)
    
    # PgPyE * DMat_val:
    # If DMat_val is SVector (6x1), result is scalar.
    # If DMat_val is SMatrix (6x3), result is 1x3 (e.g., a Transpose of an SVector).
    product_term = PgPyE * DMat_val

    BtCA_val = factor * product_term
    
    return BtCA_val
end


"""
    B(t::Float64, tCA::T, Y::SVector{6, T}, p::NamedTuple) where {T<:Real}
        -> Union{SVector{7, T}, SMatrix{7, 3, T, 21}}

Combines By and BtCA into the full sensitivity B = [By; BtCA].
Requires generic `By`, `BtCA`. Parameters `p` contain type T.
"""
function B(t::Float64, tCA::Real, Y::AbstractVector, p::NamedTuple) #where {TCAT_b<:Real}
    ByE = By(t, tCA, Y, p)     # Returns SVector{6,T} or SMatrix{6,3,T}
    BtCAE = BtCA(t, tCA, Y, p) # Returns T or Transpose{T, SVector{3,T}}

    # Combine using StaticArrays methods
    if p.ODB == 1 # ByE is 6x1, BtCAE is T
        # Ensure BtCAE is treated as a 1-element SVector for vcat
        BEvalV = vcat(ByE, SA[BtCAE])
        return BEvalV #::SVector{7, T}
    else # ByE is 6x3, BtCAE is 1x3 (Transpose)
        # vcat works directly for SMatrix and Transpose
        BEvalM = vcat(ByE, BtCAE)
        return BEvalM #::SMatrix{7, 3, T, 21}
    end
end

# ==========================================
# Hamiltonian (Generic)
# ==========================================
"""
    Ham(t::Float64, Lambda_vec::Vector{T}, Z_vec::Vector{T}, p::NamedTuple) where {T<:Real} -> T

Computes the maximised Hamiltonian H*.
Requires state Z=[Y;tCA] and costate Lambda (Vector{T}). Time t is Float64.
Uses parameter tuple `p` containing generic types T.
Requires generic `equin2cart`, `is_eclipsed`, `is_thrust_forbidden`, `B`.
"""
function Ham(t::Float64, Lambda_arg::AbstractVector, Z_arg::AbstractVector, p::NamedTuple) #where {LT<:Real, ZT<:Real}
    if length(Lambda_arg) != 7 || length(Z_arg) != 7
        error("Ham: Input Lambda and Z arguments must have 7 elements. Got $(length(Lambda_arg)) and $(length(Z_arg))")
    end
    # Convert input Vectors to StaticArrays for internal use
    Lambda = SVector{7}(Lambda_arg)
    Z      = SVector{7}(Z_arg)

    # Unpack state Z = [Y; tCA]
    # Y_equin will be SVector{6, ZT}
    Y_equin = SVector{6}(ntuple(i -> Z[i], 6))
    tCA     = Z[7] # Type ZT


    # Parameters from p needed here (expecting type T where appropriate)
    Eps_param      = p.Eps       # Type of p.Eps
    ODB_param      = p.ODB       # Int
    Mu_param       = p.Mu        # Type of p.Mu
    NuES0_param    = p.NuES0     # Type of p.NuES0
    Re_const_param = p.Re_const  # Type of p.Re_const
    #TuningP  = p.TuningP
    forbidden_intervals = p.forbidden_intervals   # ::Vector{Tuple{Float64, Float64}}

    #zero_T = zero(T); one_T = one(T)
    cart_state_Y = equin2cart(Y_equin, Mu_param) # Result elements are promoted from ZT, typeof(p.Mu)
    

    # --- Check Constraints ---
    # Hamiltonian is evaluated at time 't' and state Z(t) = [Y(t); tCA(t)]
    r_eci_t = cart_state_Y[SA[1, 2, 3]]
    is_ecl = is_eclipsed(r_eci_t, t, NuES0_param, Re_const_param)
    is_forbid = is_thrust_forbidden(t, p.forbidden_intervals)

    # --- Calculate B matrix/vector ---
    # Needs state Y (at time t), tCA (at time t), and time t
    B_val = B(t, tCA, Y_equin, p) # Returns SVector{7,T} or SMatrix{7,3,T}

    # --- Determine Optimal Control u and Hamiltonian Value ---
    # uMag = zero_T # Default scalar control magnitude
    # uVec = SA[zero_T, zero_T, zero_T] # Default vector control

    # This is the most robust and explicit way to initialize the variable.
    # It determines the "highest" numeric type needed to store the result.
    PromotedType = promote_type(eltype(Lambda_arg), eltype(Z_arg), typeof(Eps_param))
    Hamiltonian_val = zero(PromotedType)

    if !(is_ecl || is_forbid) # Only apply control if allowed
        if ODB_param == 1 # u is scalar magnitude
            InnerProd = dot(Lambda, B_val) # Lambda::SVector{7,T}, B_val::SVector{7,T} -> T
            local uMag
            # Smoothed version:
            # if TuningP <= eps(T); error("Ham: TuningP must be positive for smoothing."); end
            # uMag = one_T / (one_T + exp(-InnerProd / TuningP))

            # Bang-bang version:
            val_InnerProd = (InnerProd isa ForwardDiff.Dual) ? ForwardDiff.value(InnerProd) : InnerProd
            uMag = (val_InnerProd > eps(typeof(val_InnerProd))) ? one(InnerProd) : convert(typeof(InnerProd), 1e-8)
            Hamiltonian_val = convert(typeof(Hamiltonian_val), Eps_param) * InnerProd * uMag
        else # ODB == 0, u is 3x1 vector
            uOpt_dir = transpose(B_val) * Lambda 
            norm_uOpt_dir = norm(uOpt_dir)
            
            val_norm_uOpt_dir = (norm_uOpt_dir isa ForwardDiff.Dual) ? ForwardDiff.value(norm_uOpt_dir) : norm_uOpt_dir
            if val_norm_uOpt_dir > eps(typeof(val_norm_uOpt_dir))
                # Hamiltonian_val = Eps_param * norm_uOpt_dir
                Hamiltonian_val = convert(typeof(Hamiltonian_val), Eps_param) * norm_uOpt_dir
            end
            # else: uVec remains zero, H remains zero
        end
    end
    # else: u remains zero, Hamiltonian_val remains zero

    return Hamiltonian_val # Return the maximised Hamiltonian value H*
end

# ==========================================
# Optional Helper Functions (Generic)
# ==========================================

# In MinTimeDynamics.jl

"""
    Δr(z_vec::AbstractVector, p::NamedTuple) -> SVector{3}

Computes relative position Δr at time tCA. Input z = [y; tCA] (AbstractVector).
This version is generic and compatible with ForwardDiff.
Requires `Flow` and `equin2cart` to be AD-compatible.
`p` must contain `p.Mu` and `p.ys0`.
"""
function Δr(z_vec::AbstractVector, p::NamedTuple)
    # The check can be kept, it's good practice.
    if length(z_vec) != 7
        error("Δr: Input z vector must have 7 elements, got length $(length(z_vec))")
    end
    
    # Let the constructor infer the element type. 
    # If z_vec contains Duals, equin will be an SVector of Duals.
    # A more idiomatic way to do this is SVector(z_vec[1:6]), but ntuple also works.
    equin = SVector(z_vec[1], z_vec[2], z_vec[3], z_vec[4], z_vec[5], z_vec[6])

    # tCA will also be a Dual number during differentiation
    tCA = z_vec[7]

    # These functions must now be able to handle Dual numbers.
    # This means all operations inside Flow and equin2cart must be generic.
    yp_tca = Flow(equin, tCA, p.Mu)
    ys_tca = Flow(p.ys0, tCA, p.Mu)

    cart_p = equin2cart(yp_tca, p.Mu)
    cart_s = equin2cart(ys_tca, p.Mu)

    delta_r_pos = cart_p[SA[1, 2, 3]] - cart_s[SA[1, 2, 3]]

    # REMOVED: ::SVector{3, T}
    # The return type will be SVector{3, Dual{...}} when used by ForwardDiff,
    # and SVector{3, T} when called directly with numbers.
    return delta_r_pos
end

"""
    ∂Δr∂z(z_vec::AbstractVector{T}, p::NamedTuple) where {T<:Real} -> SMatrix{3, 7, T, 21}

Computes Jacobian ∂Δr/∂z using ForwardDiff.
Input z_vec can be Vector{T} or SVector{7,T}.
Requires generic `Δr`. `p` must contain `p.Mu` and `p.ys0`.
"""
function ∂Δr∂z(z_vec::AbstractVector{T}, p::NamedTuple) where {T<:Real}
    # Closure captures p
    # Δr_func will be called by ForwardDiff with a z of the same type as z_vec 
    # (e.g., SVector of Duals if z_vec is SVector)
    Δr_func = z -> Δr(z, p) 
    
    # ForwardDiff.jacobian works fine with SVector inputs.
    # The result J will be a standard Matrix if z_vec is Vector,
    # or an SMatrix if z_vec is SVector.
    J = ForwardDiff.jacobian(Δr_func, z_vec)
    
    # Ensure the output is always SMatrix{3, 7, T, 21} for type stability downstream
    return SMatrix{3, 7, T, 21}(J)
end











########### NOT SURE IF REALLY NEEDED THIS WAY ############
# You might want to define this struct in a common `Types.jl` file or similar
struct BPlaneTrajectory
    time_grid_dl::Vector{Float64} # Dimensionless time grid
    b1_scaled::Vector{Float64}    # Scaled B1 coordinate
    b2_scaled::Vector{Float64}    # Scaled B2 coordinate
end

# New helper function for optimal control and yaw calculation.
# This is called *during* forward integration.
function calculate_u_rsw_and_yaw(
    current_x::SVector{6, T}, # Current equinoctial state from forward integration
    t_current::Float64,       # Current physical time (or dimensionless time)
    sim_p_orig::NamedTuple,   # Original simulation parameters (contains Mu, NuES0, EclipseBool, etc.)
    NTheta_eci::SVector{3, T}, # Target direction in ECI (derived from theta_f)
    mode_odb::Int             # 1 for 1D, 0 for 3D
) where T<:Real
    # Extract relevant parameters from sim_p_orig
    Mu = sim_p_orig.Mu
    NuES0 = sim_p_orig.NuES0
    EclipseBool = sim_p_orig.EclipseBool
    LuxBool = sim_p_orig.LuxBool # Eclipse check enabled/disabled
    Re_const = sim_p_orig.Re_const
    forbidden_intervals = sim_p_orig.forbidden_intervals
    Epsilon = sim_p_orig.Epsilon # Thrust magnitude

    # OTE and Δvtca are fixed at the initial conditions of the backward propagation,
    # they don't change with current_x or t_current for the U matrix calculation.
    # OTE was pre-calculated and added to sim_p_orig (as `sim_p_orig.OTE`)
    OTE = sim_p_orig.OTE
    Δvtca = sim_p_orig.Δvtca

    zero_T = zero(T); one_T = one(T)
    
    # 1. Check for eclipse/forbidden thrust
    r_eci_t = equin2cart(current_x, Mu)[SA[1, 2, 3]]
    is_ecl = (EclipseBool == 1) ? is_eclipsed(r_eci_t, t_current, NuES0, Re_const) : false
    is_forbid = is_thrust_forbidden(t_current, forbidden_intervals)

    u_magnitude = zero_T
    thrust_dir_rsw = SA[zero_T, zero_T, zero_T]
    yaw_angle_rad = zero_T # Default, only for 1D
    switching_function_val = zero_T # dot(UE, NTheta_eci) or similar, for plotting

    if true
        # Calculate U matrix/vector. This is sensitive to the `t` argument,
        # but the `yp0` and `OTE` within `U` should be consistent with the nominal path
        # at t=0 that was used for backward integration.
        UE = U(Mu, sim_p_orig.yp0, Δvtca, OTE, t_current, NuES0, mode_odb) 

        if mode_odb == 1 # 1D Yaw-Steering Control
            if LuxBool ==   1 # Luxurious two-thruster case 
                switching_function_val = dot(UE, NTheta_eci) # Scalar
                u_magnitude = -1 
                if switching_function_val < eps(T)
                    u_magnitude = 1
                end 
                if  (is_ecl || is_forbid)
                    u_magnitude = 0
                end
                # YawThrust needs the *current* state `current_x`
                # For this comparison, assume no rate limiting for yaw.
                _, yaw_angle_rad = YawThrust(current_x, t_current, NuES0, Mu, use_rate_limit=false)

                # The yaw thrust direction in RSW
                thrust_dir_rsw = SA[zero_T, cos(yaw_angle_rad), sin(yaw_angle_rad)]
            else 
                switching_function_val = dot(UE, NTheta_eci) # Scalar
                u_magnitude = (switching_function_val > eps(T)) ? one_T : T(1e-8) # Bang-bang
                if  (is_ecl || is_forbid)
                    u_magnitude = 0
                end
                # YawThrust needs the *current* state `current_x`
                # For this comparison, assume no rate limiting for yaw.
                _, yaw_angle_rad = YawThrust(current_x, t_current, NuES0, Mu, use_rate_limit=false)

                # The yaw thrust direction in RSW
                thrust_dir_rsw = SA[zero_T, cos(yaw_angle_rad), sin(yaw_angle_rad)]
            end
        else # 0D (3D) Unconstrained Control
            # Optimal u direction is proportional to U' * NTheta
            uOpt_dir_rsw = transpose(UE) * NTheta_eci # SVector{3,T}
            norm_uOpt_dir = norm(uOpt_dir_rsw)
            
            switching_function_val = norm_uOpt_dir # Use norm as proxy for switching function for 3D
            u_magnitude = one_T # Full thrust magnitude for 3D

            if  (is_ecl || is_forbid)
                u_magnitude = 0
            end

            if norm_uOpt_dir > eps(T)
                thrust_dir_rsw = uOpt_dir_rsw / norm_uOpt_dir # Normalized optimal direction
            end
            # yaw_angle_rad is not applicable for 3D, remains zero_T
        end
    end
    
    return thrust_dir_rsw, u_magnitude, yaw_angle_rad, switching_function_val
end

# New ODE right-hand side for forward integration of the *full* state `x`
function ode_forward_dynamics!(dx::Vector{T}, x::Vector{T}, p_fwd::NamedTuple, t_fwd::Float64) where T<:Real
    # Unpack from p_fwd
    Mu = p_fwd.Mu::T
    Epsilon = p_fwd.Epsilon::T # Thrust magnitude
    sim_p_orig_for_control = p_fwd.sim_p_orig::NamedTuple # All original sim params needed for control
    NTheta_eci = p_fwd.NTheta_eci::SVector{3, T}
    mode_odb = p_fwd.mode_odb::Int

    # Convert x to SVector for type consistency with other functions
    x_svec = SVector{6,T}(x)

    # 1. Keplerian drift
    f0_val = f0(x_svec, Mu) # SVector{6, T}

    # 2. Optimal control acceleration component
    thrust_dir_rsw_unit_vec, u_magnitude_scalar, _, _ = calculate_u_rsw_and_yaw(
        x_svec, t_fwd, sim_p_orig_for_control, NTheta_eci, mode_odb
    )
    
    # Convert thrust_dir_rsw from RSW to ECI
    # RSW2ECI requires r and v in ECI
    R_p, V_p = equin2cart(x_svec, Mu)[SA[1,2,3]], equin2cart(x_svec, Mu)[SA[4,5,6]]
    thrust_dir_eci = R2EU(thrust_dir_rsw_unit_vec, R_p, V_p)

    # Acceleration in ECI frame: Epsilon * u_magnitude_scalar * thrust_dir_eci
    accel_eci = Epsilon * u_magnitude_scalar * thrust_dir_eci
    
    # Convert ECI acceleration to equinoctial elements variation (using GVEs)
    G_mat, _, _ = gve(x_svec, Mu) # G_mat is 6x3
    accel_equin_gve = G_mat * accel_eci # Result is 6x1 SVector

    # Total dynamics: Keplerian drift + GVEs contribution from thrust
    total_dx_dt = f0_val + accel_equin_gve
    
    # Assign to dx
    dx .= total_dx_dt
end


# New function to reconstruct the B-plane trajectory (Ar_B_scaled) for a given t0 and theta_f
function get_bplane_trajectory(
    t0_dl::Float64, theta_f_rad::Float64, 
    fft_coeffs_b1_t::Matrix{Complex{Float64}}, fft_coeffs_b2_t::Matrix{Complex{Float64}}, 
    common_t_grid_dc::Vector{Float64}, num_theta_fft::Int;
    num_points_on_traj::Int = 200
)
    if isempty(common_t_grid_dc) || length(common_t_grid_dc) < 2
        @warn "Time grid for B-plane trajectory reconstruction is too short or empty. Returning empty trajectory."
        return BPlaneTrajectory(Float64[], Float64[], Float64[])
    end

    # Create interpolators for p_k1(t) and p_k2(t) (re-using logic from DC.jl's animate function)
    # The `common_t_grid_dc` is ordered from 0 down to t0_physical (most negative).
    # For interpolation, we need an increasing range.
    t_interp_range = range(common_t_grid_dc[end], stop=common_t_grid_dc[1], length=length(common_t_grid_dc))
    
    # The FFT coefficients `fft_coeffs_b1_t` are (NumTheta x NumPoints),
    # where NumPoints corresponds to `common_t_grid_dc`.
    # For cubic_spline_interpolation, data `real.(vec(fft_coeffs_b1_t[k_idx, end:-1:1]))` means
    # taking the k-th row (for a specific wavenumber), reversing its time order to match `t_interp_range`.
    itp_pk1_real_list = [cubic_spline_interpolation(t_interp_range, real.(vec(fft_coeffs_b1_t[k_idx, end:-1:1])), extrapolation_bc=Throw()) for k_idx in 1:num_theta_fft]
    itp_pk1_imag_list = [cubic_spline_interpolation(t_interp_range, imag.(vec(fft_coeffs_b1_t[k_idx, end:-1:1])), extrapolation_bc=Throw()) for k_idx in 1:num_theta_fft]
    itp_pk2_real_list = [cubic_spline_interpolation(t_interp_range, real.(vec(fft_coeffs_b2_t[k_idx, end:-1:1])), extrapolation_bc=Throw()) for k_idx in 1:num_theta_fft]
    itp_pk2_imag_list = [cubic_spline_interpolation(t_interp_range, imag.(vec(fft_coeffs_b2_t[k_idx, end:-1:1])), extrapolation_bc=Throw()) for k_idx in 1:num_theta_fft]

    wavenumbers_k = get_wavenumbers(num_theta_fft)

    # Helper functions to get all coeffs for a given time from the interpolators
    function get_all_pk1_at_t(t_phys)
        return SVector{num_theta_fft, Complex{Float64}}(itp_pk1_real_list[i](t_phys) + im * itp_pk1_imag_list[i](t_phys) for i in 1:num_theta_fft)
    end
    function get_all_pk2_at_t(t_phys)
        return SVector{num_theta_fft, Complex{Float64}}(itp_pk2_real_list[i](t_phys) + im * itp_pk2_imag_list[i](t_phys) for i in 1:num_theta_fft)
    end

    # Time grid for plotting the trajectory (from t0_dl up to 0.0)
    time_grid_dl = LinRange(t0_dl, 0.0, num_points_on_traj)
    b1_scaled_path = zeros(num_points_on_traj)
    b2_scaled_path = zeros(num_points_on_traj)

    for i in 1:num_points_on_traj
        t_current = time_grid_dl[i]
        pk1_coeffs = get_all_pk1_at_t(t_current)
        pk2_coeffs = get_all_pk2_at_t(t_current)
        b1_scaled_path[i] = reconstruct_from_fft_coeffs(theta_f_rad, pk1_coeffs, wavenumbers_k)
        b2_scaled_path[i] = reconstruct_from_fft_coeffs(theta_f_rad, pk2_coeffs, wavenumbers_k)
    end

    return BPlaneTrajectory(collect(time_grid_dl), b1_scaled_path, b2_scaled_path)
end

# New function to integrate forward the full dynamics to get Δr and control profiles
function integrate_full_dynamics_for_comparison(
    results::Dict{String, Any}, time_dl_ref::Float64; num_points_fwd_int::Int = 200
)
    # This will store the outputs
    # Added SVector{3,Float64} for rsw_direction_profile_1D (actual thrust vector)
    fwd_results_1D = (Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{SVector{3,Float64}}()) # time, dr_norm, u_mag, yaw_angle, switching_func, rsw_dir_vec
    fwd_results_3D = (Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{SVector{3,Float64}}()) # time, dr_norm, u_mag, rsw_dir_vec

    for (mode_str, data) in results
        t0_dl = data["t0_dl"]
        theta_f_rad = data["theta_f_rad"]
        sim_p_current_mode = data["sim_params"] # sim_params specific to this mode (e.g., ODB)

        Mu = sim_p_current_mode.Mu
        Epsilon = sim_p_current_mode.Epsilon
        
        x_init_fwd = Flow(sim_p_current_mode.yp0, t0_dl, Mu)

        N_B_vec_T = SVector{3,eltype(x_init_fwd)}(cos(theta_f_rad), sin(theta_f_rad), 0.0)
        NTheta_eci = B2EU(N_B_vec_T, data["B1Ref"], data["B3U"])
        
        tspan_fwd_ode = (t0_dl, 0.0)
        
        p_fwd_ode = (
            Mu = Mu,
            Epsilon = Epsilon,
            sim_p_orig = sim_p_current_mode,
            NTheta_eci = NTheta_eci,
            mode_odb = sim_p_current_mode.ODB
        )
        
        prob = ODEProblem(ode_forward_dynamics!, collect(x_init_fwd), tspan_fwd_ode, p_fwd_ode)
        
        t_save_fwd = LinRange(t0_dl, 0.0, num_points_fwd_int)
        sol = solve(prob, Tsit5(), abstol=1e-9, reltol=1e-9, saveat=t_save_fwd, dense=false, dtmin=1e-6)

        if sol.retcode != SciMLBase.ReturnCode.Success
            @warn "Forward ODE for $mode_str failed with retcode: $(sol.retcode). Skipping data."
            # Fill with NaNs for consistency in array sizes if plotting will iterate over them
            nan_val = NaN
            nan_svec = SVector{3,Float64}(NaN,NaN,NaN)
            if mode_str == "1D"  || mode_str == "Eclipse" || mode_str == "One" || mode_str == "High"
                fwd_results_1D = (sol.t, fill(nan_val, length(sol.t)), fill(nan_val, length(sol.t)), fill(nan_val, length(sol.t)), fill(nan_val, length(sol.t)), fill(nan_svec, length(sol.t)))
            else
                fwd_results_3D = (sol.t, fill(nan_val, length(sol.t)), fill(nan_val, length(sol.t)), fill(nan_svec, length(sol.t)))
            end
            continue
        end

        sol_times = sol.t
        sol_states = sol.u

        delta_r_norm_profile = Float64[]
        u_mag_profile = Float64[]
        yaw_angle_profile_rad = Float64[]
        rsw_direction_profile_vecs = SVector{3,Float64}[] # Full RSW vector for either 1D or 3D
        switching_function_profile = Float64[]

        for (i, t_fwd_curr) in enumerate(sol_times)
            x_fwd_curr = SVector{6,eltype(x_init_fwd)}(sol_states[i])

            cart_p = equin2cart(x_fwd_curr, Mu)
            cart_s = equin2cart(Flow(sim_p_current_mode.ys0, t_fwd_curr, Mu), Mu) 
            delta_r_curr = cart_p[SA[1,2,3]] - cart_s[SA[1,2,3]]
            push!(delta_r_norm_profile, norm(delta_r_curr))

            # Calculate control magnitude, yaw angle, and switching function at this point
            thrust_dir_rsw_unit_vec, u_mag_curr, yaw_angle_curr, switching_func_curr = 
                calculate_u_rsw_and_yaw(x_fwd_curr, t_fwd_curr, sim_p_current_mode, NTheta_eci, sim_p_current_mode.ODB)
            
            push!(u_mag_profile, u_mag_curr)
            push!(switching_function_profile, switching_func_curr)

            # Store the *actual* thrust vector in RSW, scaled by magnitude
            push!(rsw_direction_profile_vecs, thrust_dir_rsw_unit_vec * u_mag_curr)
            
            if mode_str == "1D" || mode_str == "Eclipse" || mode_str == "One" || mode_str == "High"
                push!(yaw_angle_profile_rad, yaw_angle_curr)
            end
        end

        if mode_str == "1D" || mode_str == "Eclipse"|| mode_str == "One" || mode_str == "High"
            fwd_results_1D = (sol_times, delta_r_norm_profile, u_mag_profile, yaw_angle_profile_rad, switching_function_profile, rsw_direction_profile_vecs)
        else
            fwd_results_3D = (sol_times, delta_r_norm_profile, u_mag_profile, rsw_direction_profile_vecs)
        end
    end # end for (mode_str, odb_mode)

    return fwd_results_1D, fwd_results_3D
end