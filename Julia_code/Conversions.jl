# File: elements.jl (or Conversions.jl as used later)
# Contains functions for converting between orbital element sets (Keplerian, Equinoctial)
# and Cartesian coordinates, as well as various reference frame transformations (ECI, RSW, NTW, B-Plane, OPF).
# Refactored to use StaticArrays for performance and generic types for AD compatibility.

using LinearAlgebra # Required for I, cross, dot, norm
using StaticArrays # Required for SVector, SMatrix, SA, etc.
using Test # Used for internal checks like isapprox

# --- Helper Constants ---
# These are often used for runtime checks, keeping as Float64 might be acceptable.
# If AD fails due to these, they might need to become generic or passed as arguments.
const RTOL_NORMALIZATION = 1e-8 # Relative tolerance for norm checks
const ATOL_ORTHOGONALITY = 1e-10 # Absolute tolerance for dot product checks
const TOL_ZERO_NORM = 1e-12 # Tolerance for checking near-zero norms before normalizing

# ==========================================
# Kepler's Equation Solver
# ==========================================

"""
    kepler(l::T, ey::T, ex::T; Niter::Int = 20, tol::T = T(1e-12)) where {T<:Real} -> T

Solves Kepler's equation for true longitude `L`. Accepts generic Real type T.
"""
function kepler(l::T, ey::T, ex::T; Niter::Int = 20, tol::T = T(1e-12)) where {T<:Real}
    e_sq = ex^2 + ey^2
    one_T = one(T)
    zero_T = zero(T)
    if e_sq >= one_T; error("Kepler solver: Eccentricity squared >= 1."); end
    eta = sqrt(one_T - e_sq)
    eta_p1 = one_T + eta
    AA = one_T - ex^2 / eta_p1; BB = ey * ex / eta_p1; CC = one_T - ey^2 / eta_p1
    K = l; f = zero_T; df = one_T; converged = false
    for j in 1:Niter
        cK = cos(K); sK = sin(K)
        f = K + ey * cK - ex * sK - l
        df = one_T - ey * sK - ex * cK
        if abs(f) < tol; converged = true; break; end
        # Use eps(T) for checking near-zero derivative
        if abs(df) < eps(T); error("Kepler solver: Derivative df near zero."); end
        K -= f / df
    end
    if !converged; error("Kepler solver: Did not converge."); end
    cK = cos(K); sK = sin(K)
    sL_num = (AA * sK + BB * cK - ey); cL_num = (BB * sK + CC * cK - ex)
    L = atan(sL_num, cL_num)
    return L # Returns type T
end

# ==========================================
# Element Conversions: Equinoctial <=> Cartesian (StaticArrays)
# ==========================================

# CORRECTED SIGNATURE
function equin2cart(equin::AbstractVector{EqT}, mu::MuT; L::Union{EqT, Nothing} = nothing) where {EqT<:Real, MuT<:Real}
    # Implementation needs to handle EqT and MuT (promotion usually works)
    # Example: Use zero(EqT), one(EqT), etc. if needed
    if length(equin) != 6; error("equin2cart: Input must have 6 elements."); end
    a = equin[1]; ey = equin[2]; ex = equin[3]; qy = equin[4]; qx = equin[5]; l = equin[6]

    # Use types of inputs where appropriate
    zero_EqT = zero(EqT); one_EqT = one(EqT); two_EqT = EqT(2.0)

    if a <= zero_EqT; error("equin2cart: Semi-major axis 'a' must be positive."); end
    e_sq = ex^2 + ey^2
    if e_sq >= one_EqT; error("equin2cart: Eccentricity squared must be < 1."); end

    # Calculate L if not provided, ensure it matches the state type EQT
    L_calc::EqT = (L === nothing) ? kepler(l, ey, ex) : L
    sL = sin(L_calc); cL = cos(L_calc)

    p = a * (one_EqT - e_sq) # Result type depends on a and e_sq (EqT)
    if p <= zero(p); error("equin2cart: Semi-latus rectum 'p' must be positive."); end
    # Denominator involves EqT types
    den_rNorm = one_EqT + ex * cL + ey * sL
    if den_rNorm <= eps(EqT); error("equin2cart: Denominator for rNorm near zero."); end
    rNorm = p / den_rNorm # Type EqT

    # vHat involves MuT and p (EqT). Result type depends on promotion.
    # sqrt should handle Duals correctly.
    vHat = sqrt(mu / p)
    xPF = rNorm * cL; yPF = rNorm * sL # Type EqT
    # VxPF/VyPF involve vHat (promoted type) and EqT terms
    VxPF = -vHat * (ey + sL); VyPF = vHat * (ex + cL)

    # q terms are EqT
    qxx = qx^2; qyy = qy^2; qxy = two_EqT * qx * qy # Use EqT(2.0)
    eta_q = one_EqT + qxx + qyy

    # Frame components involve EqT types
    f_hat_x = (one_EqT - qyy + qxx) / eta_q; f_hat_y = qxy / eta_q
    g_hat_x = qxy / eta_q; g_hat_y = (one_EqT + qyy - qxx) / eta_q
    w_hat_x = -two_EqT * qy / eta_q; w_hat_y = two_EqT * qx / eta_q

    # ECI components result from promoted types
    r_x = f_hat_x * xPF + g_hat_x * yPF; v_x = f_hat_x * VxPF + g_hat_x * VyPF
    r_y = f_hat_y * xPF + g_hat_y * yPF; v_y = f_hat_y * VxPF + g_hat_y * VyPF
    r_z = w_hat_x * xPF + w_hat_y * yPF; v_z = w_hat_x * VxPF + w_hat_y * VyPF

    # Resulting SVector element type is determined by promotion of r_x, v_x etc.
    # This should naturally become Dual if EqT was Dual.
    cart = SA[r_x, r_y, r_z, v_x, v_y, v_z]
    return cart # Return type is SVector{6, <Promoted Type>}
end

# CORRECTED SIGNATURE
function cart2equin(cart::SVector{6, CartT}, mu::MuT) where {CartT<:Real, MuT<:Real}
    # Implementation needs to handle CartT and MuT
    # Ensure use of zero(CartT), one(CartT), eps(CartT) etc.
    r = cart[SA[1, 2, 3]] # SVector{3, CartT}
    v = cart[SA[4, 5, 6]] # SVector{3, CartT}

    zero_CartT = zero(CartT); one_CartT = one(CartT);
    half_CartT = CartT(0.5); two_CartT = CartT(2.0); pi_CartT = CartT(pi)

    rNorm = norm(r) # Type CartT
    vNorm = norm(v) # Type CartT
    if rNorm < eps(CartT); error("cart2equin: Position vector magnitude is near zero."); end
    if vNorm < eps(CartT); error("cart2equin: Velocity vector magnitude is near zero."); end
    rDir = r / rNorm # SVector{3, CartT}

    h = cross(r, v) # SVector{3, CartT}
    hNorm = norm(h) # Type CartT
    if hNorm < eps(CartT); error("cart2equin: Angular momentum magnitude is near zero."); end

    # e_vec involves CartT and MuT, result promoted
    e_vec = cross(v, h) / mu - rDir
    e_sq = dot(e_vec, e_vec) # Promoted type T_e = promote_type(CartT, MuT)
    if e_sq >= one(e_sq) - eps(e_sq); @warn("cart2equin: Orbit e_sq = $e_sq >= 1.0."); end
    e = sqrt(max(zero(e_sq), e_sq)) # Type T_e

    # energy involves CartT and MuT
    energy = half_CartT * vNorm^2 - mu / rNorm
    if abs(energy) < eps(typeof(energy)) && e_sq >= one(e_sq) - eps(e_sq); error("cart2equin: Parabolic orbit."); end
    a = -mu / (two_CartT * energy) # Promoted type
    if a <= zero(a) && e_sq < one(e_sq) - eps(e_sq); error("cart2equin: Non-positive 'a' ($a) for apparently elliptic orbit."); end

    zEq = h / hNorm # SVector{3, CartT}
    zECI = SA[zero_CartT, zero_CartT, one_CartT] # SVector{3, CartT}
    N = cross(zECI, zEq) # SVector{3, CartT}
    NNorm = norm(N) # Type CartT

    den_q = one_CartT + zEq[3]
    qx = zero_CartT; qy = zero_CartT;
    if den_q < eps(CartT)
        @warn("cart2equin: Inclination near 180 deg. qx, qy unstable.")
        if NNorm < eps(CartT); error("cart2equin: Retrograde equatorial orbit. qx, qy ill-defined."); end
        qx = N[1] / den_q; qy = N[2] / den_q # Type CartT
    else
        qx = N[1] / den_q; qy = N[2] / den_q # Type CartT
    end

    eta_q = one_CartT + qx^2 + qy^2

    # f-vector (xEq)
    f_hat_x = (one_CartT + qx^2 - qy^2) / eta_q
    f_hat_y = (2 * qx * qy) / eta_q
    f_hat_z = (-2 * qy) / eta_q
    xEq = SA[f_hat_x, f_hat_y, f_hat_z]

    # g-vector (yEq)
    g_hat_x = (2 * qx * qy) / eta_q
    g_hat_y = (one_CartT - qx^2 + qy^2) / eta_q
    g_hat_z = (2 * qx) / eta_q
    yEq = SA[g_hat_x, g_hat_y, g_hat_z]

    # 4. NOW calculate ex and ey using the correct frame
    ey = dot(e_vec, yEq)
    ex = dot(e_vec, xEq)

    # L involves CartT
    rNorm = norm(r)
    r_xEq = dot(r, xEq); r_yEq = dot(r, yEq)
    sL = r_yEq / rNorm; cL = r_xEq / rNorm
    L = atan(sL, cL) # Type CartT

    l = zero(L) # Initialize l with type CartT
    if e < eps(typeof(e))
        l = L
    else
        # Calculations for K and l involve promoted types (from e_sq, a, etc.)
        eta = sqrt(max(zero(e_sq), one(e_sq) - e_sq))
        if abs(eta) < eps(typeof(eta)); error("cart2equin: eta is near zero (e near 1)."); end
        eta_p1 = one(eta) + eta
        # AA, BB, CC involve ex, ey (promoted type) and eta
        AA = one(eta) - ex^2 / eta_p1; BB = ey * ex / eta_p1; CC = one(eta) - ey^2 / eta_p1
        den_r_over_a = one(L) + ex * cL + ey * sL
        if den_r_over_a <= eps(typeof(den_r_over_a)); error("cart2equin: Denominator for r/a calculation near zero or negative."); end
        # r_over_a involves e_sq and den_r_over_a (promoted types)
        r_over_a = (one(e_sq) - e_sq) / den_r_over_a

        # sLmod, cLmod involve r_over_a, sL, cL, ey, ex (promoted types)
        sLmod = r_over_a * sL + ey; cLmod = r_over_a * cL + ex
        det = eta # Type of eta
        # sK, cK involve AA, BB, CC, sLmod, cLmod, det (promoted types)
        sK = (CC * sLmod - BB * cLmod) / det; cK = (-BB * sLmod + AA * cLmod) / det
        norm_scK = sqrt(sK^2 + cK^2) # Promoted type
        K = (norm_scK < eps(norm_scK)) ? zero(norm_scK) : atan(sK / norm_scK, cK / norm_scK) # Promoted type
        # l involves K, ey, ex (promoted types)
        l = mod(K + ey * cos(K) - ex * sin(K), two_CartT * pi_CartT) # Use CartT for 2pi
    end

    # Final elements will have promoted types
    equin = SA[a, ey, ex, qy, qx, l]
    # Assert final type if needed, but StaticArrays should handle promotion
    # return equin::SVector{6, promote_type(CartT, MuT)}? Maybe too complex.
    # Let return type be inferred.
    return equin
end
# ==========================================
# Element Conversions: Equinoctial <=> Keplerian
# ==========================================

"""
    equin2kepl(equin::AbstractVector{T}) where {T<:Real} -> Vector{T}

Converts equinoctial elements (AbstractVector of type T) to Keplerian elements (Vector of type T).
"""
function equin2kepl(equin::AbstractVector{T}) where {T<:Real}
    if length(equin) != 6; error("equin2kepl: Input must have 6 elements."); end
    a=equin[1]; ey=equin[2]; ex=equin[3]; qy=equin[4]; qx=equin[5]; l=equin[6]

    zero_T = zero(T); two_T = T(2.0); pi_T = T(pi)

    e=sqrt(ex^2+ey^2)
    q_sq=qy^2+qx^2
    # Use max to avoid sqrt of small negative number due to precision
    i=two_T*atan(sqrt(max(zero_T, q_sq)))
    # Use eps(T) for zero check
    Ω=(q_sq < eps(T)) ? zero_T : atan(qy,qx)
    Ω=mod(Ω, two_T * pi_T)
    # Use eps(T) for zero check
    pomega=(e < eps(T)) ? Ω : atan(ey,ex)
    ω=mod(pomega-Ω, two_T * pi_T)
    M=mod(l-pomega, two_T * pi_T)
    # Output is standard Vector, but elements are type T
    return T[a, e, i, Ω, ω, M]
end

"""
    kepl2equin(kepl::AbstractVector{T}) where {T<:Real} -> SVector{6, T}

Converts Keplerian elements (AbstractVector of type T) to equinoctial elements (SVector of type T).
"""
function kepl2equin(kepl::AbstractVector{T}) where {T<:Real}
     if length(kepl) != 6; error("kepl2equin: Input must have 6 elements."); end
     a=kepl[1]; e=kepl[2]; i=kepl[3]; Ω=kepl[4]; ω=kepl[5]; M=kepl[6]

     zero_T = zero(T); one_T = one(T); two_T = T(2.0); pi_T = T(pi)

     if e < zero_T; error("kepl2equin: Eccentricity 'e' cannot be negative."); end
     if e >= one_T; @warn("kepl2equin: Eccentricity 'e' >= 1."); end
     # Allow slight tolerance for pi due to floating point representation?
     # Or rely on caller to ensure i is valid. Sticking to strict check for now.
     if i < zero_T || i > pi_T; error("kepl2equin: Inclination 'i' must be in [0, π]."); end

     pomega=mod(Ω+ω, two_T * pi_T)
     ey=e*sin(pomega); ex=e*cos(pomega)
     ti2=tan(i/two_T) # Use generic two
     qy=ti2*sin(Ω); qx=ti2*cos(Ω)
     l=mod(M+pomega, two_T * pi_T)
     # Output SVector with elements of type T
     return SA[a, ey, ex, qy, qx, l]
end

# ==========================================
# Reference Frame Transformations (StaticArrays)
# ==========================================

"""
    safe_normalize(v::SVector{N, T}, vec_name::String) where {N, T<:Real} -> SVector{N, T}

Normalizes an SVector of type T, checking for zero norm using TOL_ZERO_NORM (Float64).
"""
function safe_normalize(v::SVector{N, T}, vec_name::String) where {N, T<:Real}
    n = norm(v)
    # Runtime check using Float64 tolerance. If AD fails here, make tolerance generic.
    if n < TOL_ZERO_NORM
        error("Cannot normalize near-zero vector: $vec_name = $v (norm=$n)")
    end
    return v / n # Returns SVector{N, T}
end

# --- B-Plane <=> ECI ---
""" B2EU: Transforms B-Plane SVector{3, T} to ECI SVector{3, T}. """
function B2EU(BVec::SVector{3, BtT}, B1Ref::SVector{3, B1tT}, B3::SVector{3, B3tT}) where {BtT<:Real, B1tT<:Real, B3tT<:Real}
    B3U = safe_normalize(B3, "B3")
    B1U_unnormalized = B1Ref - dot(B1Ref, B3U) * B3U
    B1U = safe_normalize(B1U_unnormalized, "Projected B1Ref")
    B2U = cross(B3U, B1U)
    # Build SMatrix column-wise (will be SMatrix{3, 3, T, 9})
    RMat = SMatrix{3, 3}(B1U..., B2U..., B3U...)
    EVec = RMat * BVec
    return EVec # Returns SVector{3, T}
end

""" E2BU: Transforms ECI SVector{3, T} to B-Plane SVector{3, T}. """
function E2BU(EVec::SVector{3, EtT}, B1Ref::SVector{3, B1tT}, B3::SVector{3, B3tT}) where {EtT<:Real, B1tT<:Real, B3tT<:Real}
    B3U = safe_normalize(B3, "B3")
    B1U_unnormalized = B1Ref - dot(B1Ref, B3U) * B3U
    B1U = safe_normalize(B1U_unnormalized, "Projected B1Ref")
    B2U = cross(B3U, B1U)
    # Build transpose of RMat directly (will be SMatrix{3, 3, T, 9})
    RMat_inv = transpose(SMatrix{3, 3}(B1U..., B2U..., B3U...))
    BVec = RMat_inv * EVec
    return BVec # Returns SVector{3, T}
end

# --- RSW <=> ECI ---
""" R2EU: Transforms RSW SVector{3, T} to ECI SVector{3, T}. """
function R2EU(RVec::SVector{3, RVtT}, Rp::SVector{3, RtT}, Vp::SVector{3, VtT}) where {RVtT<:Real, RtT<:Real, VtT<:Real}
    RU = safe_normalize(Rp, "Rp")
    h = cross(Rp, Vp)
    WU = safe_normalize(h, "h = Rp x Vp")
    SU = cross(WU, RU)
    RMat = SMatrix{3, 3}(RU..., SU..., WU...) # Column-wise

    EVec = RMat * RVec

    return EVec # Returns SVector{3, T}
end

""" E2RU: Transforms ECI SVector{3, T} to RSW SVector{3, T}. """
function E2RU(EVec::SVector{3, EtT}, Rp::SVector{3, RtT}, Vp::SVector{3, VtT}) where {EtT<:Real, RtT<:Real, VtT<:Real}
    RU = safe_normalize(Rp, "Rp")
    h = cross(Rp, Vp)
    WU = safe_normalize(h, "h = Rp x Vp")
    SU = cross(WU, RU)
    RMat_inv = transpose(SMatrix{3, 3}(RU..., SU..., WU...))
    RVec = RMat_inv * EVec
    return RVec # Returns SVector{3, T}
end

# --- NTW <=> ECI ---
""" N2EU: Transforms NTW SVector{3, T} to ECI SVector{3, T}. """
function N2EU(NVec::SVector{3, NtT}, Rp::SVector{3, RtT}, Vp::SVector{3, VtT}) where {NtT<:Real, RtT<:Real, VtT<:Real}
    TU = safe_normalize(Vp, "Vp")
    h = cross(Rp, Vp)
    HU = safe_normalize(h, "h = Rp x Vp") # W-axis
    NU = cross(TU, HU)
    RMat = SMatrix{3, 3}(NU..., TU..., HU...) # Column-wise
    EVec = RMat * NVec
    return EVec # Returns SVector{3, T}
end

""" E2NU: Transforms ECI SVector{3, T} to NTW SVector{3, T}. """
function E2NU(EVec::SVector{3, EtT}, Rp::SVector{3, RtT}, Vp::SVector{3, VtT}) where {EtT<:Real, RtT<:Real, VtT<:Real}
    TU = safe_normalize(Vp, "Vp")
    h = cross(Rp, Vp)
    HU = safe_normalize(h, "h = Rp x Vp") # W-axis
    NU = cross(TU, HU)
    RMat_inv = transpose(SMatrix{3, 3}(NU..., TU..., HU...))
    NVec = RMat_inv * EVec
    return NVec # Returns SVector{3, T}
end

function LVLHCorrection(LVLHThalesVec)
    RU                  = SVector{3}(0, 0, -1)
    SU                  = SVector{3}(1, 0,  0)
    WU                  = SVector{3}(0, -1, 0)
    RMat                = SMatrix{3, 3}(RU..., SU..., WU...) 
    NewVec              = RMat*LVLHThalesVec;
    return NewVec
end

# --- OPF <=> ECI ---
""" E2OU: Transforms ECI SVector{3, T} to OPF SVector{3, T}. """
function E2OU(EVec::SVector{3, EtT}, h::SVector{3, HtT}, rES_unit::SVector{3, RtT}) where {EtT<:Real, HtT<:Real, RtT<:Real}
    # Runtime check using Float64 tolerance.
    if !isapprox(norm(rES_unit), 1.0, rtol=RTOL_NORMALIZATION); @warn "E2OU: Input rES_unit should be normalized."; end
    O3U = safe_normalize(h, "h")
    O2U_unnormalized = cross(O3U, rES_unit)
    O2U = safe_normalize(O2U_unnormalized, "O3 x rES_unit")
    O1U = cross(O2U, O3U)
    RMat_inv = transpose(SMatrix{3, 3}(O1U..., O2U..., O3U...))
    OVec = RMat_inv * EVec

    return OVec # Returns SVector{3, T}
end

""" O2EU: Transforms OPF SVector{3, T} to ECI SVector{3, T}. """
function O2EU(OVec::SVector{3, OtT}, h::SVector{3, HtT}, rES_unit::SVector{3, RtT}) where {OtT<:Real, HtT<:Real, RtT<:Real}
     # Runtime check using Float64 tolerance.
    if !isapprox(norm(rES_unit), 1.0, rtol=RTOL_NORMALIZATION); @warn "O2EU: Input rES_unit should be normalized."; end
    O3U = safe_normalize(h, "h")
    O2U_unnormalized = cross(O3U, rES_unit)
    O2U = safe_normalize(O2U_unnormalized, "O3 x rES_unit")
    O1U = cross(O2U, O3U)
    RMat = SMatrix{3, 3}(O1U..., O2U..., O3U...) # Column-wise
    EVec = RMat * OVec
    return EVec # Returns SVector{3, T}
end

# ==========================================
# B-Plane Related Helpers (StaticArrays)
# ==========================================
"""
    B1B3(Δvtca::SVector{3, T}, yp0_eq::AbstractVector{T}, Mu::T) where {T<:Real}
        -> Tuple{SVector{3, T}, SVector{3, T}}

Calculates B-plane basis vectors B1 (unit) and B3 (unit) using generic type T.
"""
function B1B3(Δvtca::SVector{3, VecT}, yp0_eq::AbstractVector{EqT}, Mu::MuT) where {VecT<:Real, EqT<:Real, MuT<:Real}
    # Implementation calls safe_normalize, equin2cart (corrected), dot
    # Type T for output will be promotion of VecT, EqT, MuT
    B3U = safe_normalize(Δvtca, "Δvtca") # SVector{3, VecT}
    # equin2cart takes EqT, MuT -> returns SVector{6, promote(EqT,MuT)}
    yp0_cart = equin2cart(yp0_eq, Mu)
    # Element type of yp0_cart determines vy0's type
    vy0 = yp0_cart[SA[4, 5, 6]]
    # dot result type depends on vy0 and B3U
    dot_prod = dot(vy0, B3U)
    B1Ref_unnormalized = vy0 - dot_prod * B3U
    B1Ref = safe_normalize(B1Ref_unnormalized, "Projected vy0")
    # Return types determined by promotion
    return B1Ref, B3U
end

"""
    Theta2NThetaE(Δvtca::SVector{3, T}, yp0_eq::AbstractVector{T},
                  Mu::T, Theta::T) where {T<:Real} -> SVector{3, T}

Calculates ECI direction vector for a B-plane angle Theta using generic type T.
"""
function Theta2NThetaE(Δvtca::SVector{3, VecT}, yp0_eq::AbstractVector{EqT},
                       Mu::MuT, Theta::ThetaT) where {VecT<:Real, EqT<:Real, MuT<:Real, ThetaT<:Real}
    # Calls B1B3 (corrected), cos, sin, SA, B2EU
    # Output type T will be promotion of VecT, EqT, MuT, ThetaT
    B1Ref, B3U = B1B3(Δvtca, yp0_eq, Mu) # Types depend on promotion inside B1B3
    cosT = cos(Theta); sinT = sin(Theta) # Type ThetaT
    zero_ThetaT = zero(ThetaT)
    NThetaB = SA[cosT, sinT, zero_ThetaT] # SVector{3, ThetaT}
    # B2EU takes NThetaB(ThetaT), B1Ref(Promoted), B3U(VecT)
    # -> Output type determined by promotion inside B2EU
    NThetaE = B2EU(NThetaB, B1Ref, B3U)
    return NThetaE
end


# ==========================================
# Cartesian <=> Polar Coordinates (2D)
# ==========================================
""" c2p(x::T, y::T) where {T<:Real} -> Tuple{T, T} """
function c2p(x::T, y::T) where {T<:Real}
    zero_T = zero(T); two_T = T(2.0); pi_T = T(pi)
    r = sqrt(x^2 + y^2); θ = atan(y, x)
    if θ < zero_T; θ += two_T * pi_T; end # Use generic 2pi
    return r, θ # Returns Tuple{T, T}
end

""" p2c(r::T, θ::T) where {T<:Real} -> Tuple{T, T} """
function p2c(r::T, θ::T) where {T<:Real}
    zero_T = zero(T)
    if r < zero_T; error("p2c: Radius 'r' cannot be negative."); end
    x = r * cos(θ); y = r * sin(θ)
    return x, y # Returns Tuple{T, T}
end

# ==========================================
# Initialization Helper (StaticArrays)
# ==========================================
"""
    ys0FromInit(yp0_tca::AbstractVector{T}, Mu::T,
                Δvtca::SVector{3, T}) where {T<:Real} -> SVector{6, T}

Calculates the initial equinoctial state (`ys0_tca`) of the secondary at TCA using generic type T.
"""
# CORRECTED SIGNATURE
function ys0FromInit(yp0_tca::AbstractVector{EqT}, Mu::MuT,
                     Δvtca::SVector{3, VecT}) where {EqT<:Real, MuT<:Real, VecT<:Real}
    # Calls equin2cart (corrected), SA, vcat, cart2equin (corrected)
    # Output type T will be promotion of EqT, MuT, VecT
    # equin2cart takes EqT, MuT -> SVector{6, promote(EqT,MuT)}
    xp_cart_tca = equin2cart(yp0_tca, Mu)
    rp_tca = xp_cart_tca[SA[1, 2, 3]]
    vp_tca = xp_cart_tca[SA[4, 5, 6]]
    # rs_tca has same type as rp_tca
    rs_tca = rp_tca
    # vs_tca involves types of vp_tca and Δvtca (VecT) -> promoted type
    vs_tca = vp_tca - Δvtca
    # xs_cart_tca element type is promotion of rs_tca and vs_tca types
    xs_cart_tca = vcat(rs_tca, vs_tca)
    # cart2equin takes xs_cart_tca and Mu (MuT)
    # -> Output type determined by promotion inside cart2equin
    ys0_tca = cart2equin(xs_cart_tca, Mu)
    return ys0_tca
end