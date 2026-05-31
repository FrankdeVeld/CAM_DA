# Geometries.txt
# Defines different orbital geometry scenarios for rendezvous/collision analysis.
# Refactored to use StaticArrays where appropriate.

# Assumes conversions.jl provides StaticArray compatible versions
include("Conversions.jl")
using LinearAlgebra
using StaticArrays

"""
    GeometryPicker(Scn::Int, Mu::Float64;
                   aC::Union{Float64, Nothing} = nothing,
                   eC::Union{Float64, Nothing} = nothing,
                   iC::Union{Float64, Nothing} = nothing,
                   OmC::Union{Float64, Nothing} = nothing,
                   omC::Union{Float64, Nothing} = nothing,
                   nuC::Union{Float64, Nothing} = nothing,
                   deltaVcMag::Union{Float64, Nothing} = nothing,
                   deltaVAngA::Union{Float64, Nothing} = nothing,
                   deltaVAngB::Union{Float64, Nothing} = nothing)
                   -> Tuple{SVector{6, Float64}, SVector{6, Float64}, SVector{3, Float64}}

Scenario Definitions:

Scn = 1: Primary defined by Keplerian elements (aC...nuC) at TCA.
Relative velocity defined by magnitude (deltaVcMag) and angles
(deltaVAngA, deltaVAngB) assumed to define the NTW components
[N, T, W] as deltaVcMag * [sin(A)sin(B), cos(B), cos(A)sin(B)].
Requires aC...nuC, deltaVcMag, deltaVAngA, deltaVAngB keywords.

Scn = 2: Toy model; circular orbits, specific Cartesian setup. Ignores keywords.

Scn = 3: ISS-like orbits defined by specific Cartesian states. Ignores keywords.

Scn = 4: Lamberto validation case; primary defined by Equinoctial elements.
Relative velocity by magnitude/angles (NTW assumed as in Scn 1). Ignores keywords.

Scn = 5: Circular primary, specific Cartesian setup with different deltaV. Ignores keywords.

Scn = 6: Circular primary, specific Cartesian setup with different deltaV. Ignores keywords.

Selects and computes initial orbital states (StaticArrays) for primary (yp0) and
secondary (ys0) objects at TCA, and the relative velocity vector (deltaV_ECI)
based on the chosen scenario (Scn).

Ellipticity Check:
The function verifies that both the resulting primary (yp0) and secondary (ys0)
orbits are elliptic (e < 1). If either is parabolic or hyperbolic (e >= 1),
it throws an error. This often indicates the specified deltaVcMag is too large
for the given primary orbit, leading to an escape trajectory for the secondary.

Arguments

    Scn::Int: Scenario identifier (1-6).

    Mu::Float64: Gravitational parameter.

    aC, eC, iC, OmC, omC, nuC: Keplerian elements for primary (required for Scn=1). Angles in radians.

    deltaVcMag: Magnitude of the relative velocity (required for Scn=1, Scn=4 uses fixed value).

    deltaVAngA, deltaVAngB: Angles defining deltaV direction in NTW frame (radians, interpretation specific to scenario).

# Returns
- `Tuple{SVector{6, Float64}, SVector{6, Float64}, SVector{3, Float64}}`:
    - `yp0`: Equinoctial elements of the primary at TCA (StaticArray).
    - `ys0`: Equinoctial elements of the secondary at TCA (StaticArray).
    - `deltaV_ECI`: Relative velocity vector (Primary - Secondary) in ECI (StaticArray).

# Throws

    ArgumentError: If required keyword arguments for Scn=1 are missing.

    ArgumentError: If resulting primary or secondary orbit eccentricity is >= 1.0.

    Errors from underlying conversion functions (e.g., invalid elements).
"""
function GeometryPicker(Scn::Int, Mu::Float64;
                        aC::Union{Float64, Nothing}=nothing,
                        eC::Union{Float64, Nothing}=nothing,
                        iC::Union{Float64, Nothing}=nothing,
                        OmC::Union{Float64, Nothing}=nothing,
                        omC::Union{Float64, Nothing}=nothing,
                        nuC::Union{Float64, Nothing}=nothing,
                        deltaVcMag::Union{Float64, Nothing}=nothing,
                        deltaVAngA::Union{Float64, Nothing}=nothing,
                        deltaVAngB::Union{Float64, Nothing}=nothing)

    # Initialize output variables (use zero StaticArrays)
    yp0::SVector{6, Float64} = zero(SVector{6, Float64})
    ys0::SVector{6, Float64} = zero(SVector{6, Float64})
    deltaV_ECI::SVector{3, Float64} = zero(SVector{3, Float64})
    xStca::SVector{6, Float64} = zero(SVector{6, Float64})
    xS_vel_tca::SVector{3, Float64} = zero(SVector{3, Float64})

    xptca::SVector{6, Float64} = zero(SVector{6, Float64})
    xstca::SVector{6, Float64} = zero(SVector{6, Float64})
    R_EARTH = 6.378e6
    if Scn == 1
        req_vars = (aC, eC, iC, OmC, omC, nuC, deltaVcMag, deltaVAngA, deltaVAngB)
        if any(isnothing, req_vars); error("GeometryPicker Scn=1 requires keywords..."); end

        E_C = 2 * atan(sqrt(max(0.0,(1.0 - eC) / (1.0 + eC))) * tan(nuC / 2.0)) # Added max(0,...) for robustness
        M_C = mod(E_C - eC * sin(E_C), 2pi)

        xCKep::Vector{Float64} = [aC, eC, iC, OmC, omC, M_C]
        yp0 = kepl2equin(xCKep) # Returns SVector{6}
        xCCart::SVector{6, Float64} = equin2cart(yp0, Mu) # Returns SVector{6}

        A = deltaVAngA; B = deltaVAngB
        # deltaVcNTW components (standard vector first)
        deltaVcNTW_vec = [sin(A) * sin(B), cos(B), cos(A) * sin(B)] * deltaVcMag
        # Convert deltaV from NTW frame to ECI (using StaticArray inputs/outputs)
        # deltaV_ECI = N2EU(SVector{3}(deltaVcNTW_vec), xCCart[SA[1,2,3]], xCCart[SA[4,5,6]])
        deltaV_ECI = R2EU(SVector{3}(deltaVcNTW_vec), xCCart[SA[1,2,3]], xCCart[SA[4,5,6]])

        # Calculate secondary state at TCA
        xS_vel_tca = xCCart[SA[4,5,6]] - deltaV_ECI
        xStca  = vcat(xCCart[SA[1,2,3]], xS_vel_tca)
        ys0 = cart2equin(xStca, Mu) # Returns SVector{6}

    elseif Scn == 2
        xp_pos_tca = SA[7e6, 0.0, 0.0] # Use SA literal
        xs_pos_tca = xp_pos_tca
        r_norm = norm(xp_pos_tca)
        v_mag_circ = sqrt(Mu / r_norm)
        xp_vel_tca = SA[0.0, v_mag_circ / sqrt(2.0), v_mag_circ / sqrt(2.0)]
        deltaV_ECI_init = SA[0.0, 0.0, 10000.0] # Large deltaV
        xs_vel_tca = xp_vel_tca - deltaV_ECI_init
        xptca = vcat(xp_pos_tca, xp_vel_tca)
        xstca = vcat(xs_pos_tca, xs_vel_tca)
        yp0 = cart2equin(xptca, Mu) # Returns SVector{6}
        ys0 = cart2equin(xstca, Mu) # Returns SVector{6}
        # deltaV_ECI calculated at end

    elseif Scn == 3
        xstca_vec = [-1843631.11; -6438727.45; -1034886.61; 4872.312236; -450.317177; -5911.324142]
        xptca_vec = [-1843631.11; -6438727.45; -1034886.61; 4872.312236;  450.317177; -5911.324142]
        xstca  = SVector{6}(xstca_vec)
        xptca  = SVector{6}(xptca_vec)
        yp0 = cart2equin(xptca, Mu) # Returns SVector{6}
        ys0 = cart2equin(xstca, Mu) # Returns SVector{6}
        deltaV_ECI = xptca[SA[4,5,6]] - xstca[SA[4,5,6]] # StaticArray subtraction

    elseif Scn == 4
        yp0_vec = [1.0 * (6.378 + 0.5) * 1e6; 0.0; 0.05; 0.0; 0.0; deg2rad(30.0)]
        yp0 = SVector{6}(yp0_vec) # Primary state as SVector
        xptca = equin2cart(yp0, Mu) # Returns SVector{6}

        alpha = deg2rad(40.0 + 180.0); beta = deg2rad(140.0); magDeltaV = 5e3
        deltaVcNTW_vec = [sin(alpha) * sin(beta), cos(beta), cos(alpha) * sin(beta)] * magDeltaV
        #deltaV_ECI = N2EU(SVector{3}(deltaVcNTW_vec), xptca[SA[1,2,3]], xptca[SA[4,5,6]])
        deltaV_ECI = R2EU(SVector{3}(deltaVcNTW_vec), xptca[SA[1,2,3]], xptca[SA[4,5,6]])

        xs_vel_tca = xptca[SA[4,5,6]] - deltaV_ECI
        xstca = vcat(xptca[SA[1,2,3]], xs_vel_tca)
        ys0 = cart2equin(xstca, Mu) # Returns SVector{6}

    elseif Scn == 5
        xp_pos_tca = SA[7e6, 0.0, 0.0]; xs_pos_tca = xp_pos_tca
        r_norm = norm(xp_pos_tca); v_mag_circ = sqrt(Mu / r_norm)
        xp_vel_tca = SA[0.0, v_mag_circ, 0.0]
        deltaV_ECI_init = SA[0.0, 5e3, 5e3]
        xs_vel_tca = xp_vel_tca - deltaV_ECI_init
        xptca = vcat(xp_pos_tca, xp_vel_tca)
        xstca = vcat(xs_pos_tca, xs_vel_tca)
        yp0 = cart2equin(xptca, Mu); ys0 = cart2equin(xstca, Mu)
        # deltaV_ECI calculated at end

    elseif Scn == 6
        xp_pos_tca = SA[7e6, 0.0, 0.0]; xs_pos_tca = xp_pos_tca
        r_norm = norm(xp_pos_tca); v_mag_circ = sqrt(Mu / r_norm)
        xp_vel_tca = SA[0.0, v_mag_circ / sqrt(2.0) * 1.01, v_mag_circ / sqrt(2.0)]
        deltaV_ECI_init = SA[0.0, 0.0, 8e3]
        xs_vel_tca = xp_vel_tca - deltaV_ECI_init
        xptca = vcat(xp_pos_tca, xp_vel_tca)
        xstca = vcat(xs_pos_tca, xs_vel_tca)
        yp0 = cart2equin(xptca, Mu); ys0 = cart2equin(xstca, Mu)
        # deltaV_ECI calculated at end
    elseif Scn == 7 # ISS-like Orbit
        # Typical ISS parameters: a ~ 6793 km (415km alt), e ~ 0.007, i ~ 51.6 deg
        # nuC will be passed via keyword for variation.
        # RAAN (OmC), AOP (omC) can be set to 0 or other values if needed.
        if isnothing(aC) || isnothing(eC) || isnothing(iC) || isnothing(nuC)
            error("GeometryPicker Scn=7 (ISS) requires aC, eC, iC, nuC to be specified. Pass them via keywords from typical ISS values.")
        end
        
        _OmC = isnothing(OmC) ? 0.0 : OmC
        _omC = isnothing(omC) ? 0.0 : omC
        _deltaVcMag = isnothing(deltaVcMag) ? 10.0 : deltaVcMag # 10 m/s relative speed
        _deltaVAngA = isnothing(deltaVAngA) ? deg2rad(0.0) : deltaVAngA   # Mostly along-track collision
        _deltaVAngB = isnothing(deltaVAngB) ? deg2rad(90.0) : deltaVAngB                                                         
        _deltaVAngA_actual = isnothing(deltaVAngA) ? 0.0 : deltaVAngA # Angle in R-W plane for component perpendicular to T
        _deltaVAngB_actual = isnothing(deltaVAngB) ? 0.0 : deltaVAngB # Angle determining T component vs R-W component magnitude. Let cos(B) be T-component factor.

        E_C = 2 * atan(sqrt(max(0.0,(1.0 - eC) / (1.0 + eC))) * tan(nuC / 2.0))
        M_C = mod(E_C - eC * sin(E_C), 2pi)

        xCKep = [aC, eC, iC, _OmC, _omC, M_C]
        yp0 = kepl2equin(xCKep)
        xCCart = equin2cart(yp0, Mu)

        _deltaVAngA_param = isnothing(deltaVAngA) ? 0.0 : deltaVAngA # Default: No component in W direction if B is small
        _deltaVAngB_param = isnothing(deltaVAngB) ? deg2rad(5.0) : deltaVAngB # Small B for T-dominant deltaV

        deltaVcNTW_vec = [sin(_deltaVAngA_param) * sin(_deltaVAngB_param), cos(_deltaVAngB_param), cos(_deltaVAngA_param) * sin(_deltaVAngB_param)] * _deltaVcMag
        # deltaV_ECI = N2EU(SVector{3}(deltaVcNTW_vec), xCCart[SA[1,2,3]], xCCart[SA[4,5,6]])
        deltaV_ECI = R2EU(SVector{3}(deltaVcNTW_vec), xCCart[SA[1,2,3]], xCCart[SA[4,5,6]])

        xS_vel_tca = xCCart[SA[4,5,6]] - deltaV_ECI # Primary - Secondary = deltaV_ECI => V_sec = V_prim - deltaV_ECI
        xStca  = vcat(xCCart[SA[1,2,3]], xS_vel_tca) # Secondary pos = Primary pos at TCA

        ys0 = cart2equin(xStca, Mu)

        if equin2kepl(ys0)[1]*(1-equin2kepl(yp0)[2]) < R_EARTH + 200e3
            error("Perigee radius too small")
        end


    elseif Scn == 8 # Iridium NEXT-like Orbit
        # Typical Iridium: a ~ 7158 km (780km alt), e ~ 0.0001, i ~ 86.4 deg
        if isnothing(aC) || isnothing(eC) || isnothing(iC) || isnothing(nuC)
            error("GeometryPicker Scn=8 (Iridium) requires aC, eC, iC, nuC.")
        end
        _OmC = isnothing(OmC) ? 0.0 : OmC
        _omC = isnothing(omC) ? 0.0 : omC
        _deltaVcMag = isnothing(deltaVcMag) ? 10.0 : deltaVcMag
        _deltaVAngA_param = isnothing(deltaVAngA) ? 0.0 : deltaVAngA 
        _deltaVAngB_param = isnothing(deltaVAngB) ? deg2rad(5.0) : deltaVAngB

        E_C = 2 * atan(sqrt(max(0.0,(1.0 - eC) / (1.0 + eC))) * tan(nuC / 2.0))
        M_C = mod(E_C - eC * sin(E_C), 2pi)
        xCKep = [aC, eC, iC, _OmC, _omC, M_C]
        yp0 = kepl2equin(xCKep)
        xCCart = equin2cart(yp0, Mu)

        deltaVcNTW_vec = [sin(_deltaVAngA_param)*sin(_deltaVAngB_param), cos(_deltaVAngB_param), cos(_deltaVAngA_param)*sin(_deltaVAngB_param)] * _deltaVcMag
        # deltaV_ECI = N2EU(SVector{3}(deltaVcNTW_vec), xCCart[SA[1,2,3]], xCCart[SA[4,5,6]])
        deltaV_ECI = R2EU(SVector{3}(deltaVcNTW_vec), xCCart[SA[1,2,3]], xCCart[SA[4,5,6]])
        xS_vel_tca = xCCart[SA[4,5,6]] - deltaV_ECI
        xStca = vcat(xCCart[SA[1,2,3]], xS_vel_tca)
        ys0 = cart2equin(xStca, Mu)

    elseif Scn == 9 # Falcon 9 GTO-like Orbit
        # Typical GTO: perigee alt 400km, apogee alt 35786km, i ~ 27 deg
        # R_earth = 6378.1e3 m
        # rp = 6378.1e3 + 400e3; ra = 6378.1e3 + 35786e3
        # a_gto = (rp+ra)/2; e_gto = (ra-rp)/(ra+rp)
        # a_gto ~ 24382.1e3 m; e_gto ~ 0.73
        if isnothing(aC) || isnothing(eC) || isnothing(iC) || isnothing(nuC)
            error("GeometryPicker Scn=9 (GTO) requires aC, eC, iC, nuC.")
        end
        _OmC = isnothing(OmC) ? 0.0 : OmC
        _omC = isnothing(omC) ? 0.0 : omC
        _deltaVcMag = isnothing(deltaVcMag) ? 100.0 : deltaVcMag # Higher rel speed for GTO encounters
        _deltaVAngA_param = isnothing(deltaVAngA) ? 0.0 : deltaVAngA
        _deltaVAngB_param = isnothing(deltaVAngB) ? deg2rad(5.0) : deltaVAngB

        E_C = 2 * atan(sqrt(max(0.0,(1.0 - eC) / (1.0 + eC))) * tan(nuC / 2.0))
        M_C = mod(E_C - eC * sin(E_C), 2pi)
        xCKep = [aC, eC, iC, _OmC, _omC, M_C]
        yp0 = kepl2equin(xCKep)
        xCCart = equin2cart(yp0, Mu)
        deltaVcNTW_vec = [sin(_deltaVAngA_param)*sin(_deltaVAngB_param), cos(_deltaVAngB_param), cos(_deltaVAngA_param)*sin(_deltaVAngB_param)] * _deltaVcMag
        #deltaV_ECI = N2EU(SVector{3}(deltaVcNTW_vec), xCCart[SA[1,2,3]], xCCart[SA[4,5,6]])
        deltaV_ECI = R2EU(SVector{3}(deltaVcNTW_vec), xCCart[SA[1,2,3]], xCCart[SA[4,5,6]])
        xS_vel_tca = xCCart[SA[4,5,6]] - deltaV_ECI
        xStca = vcat(xCCart[SA[1,2,3]], xS_vel_tca)
        ys0 = cart2equin(xStca, Mu)
    else
        error("GeometryPicker: Invalid Scenario number Scn = $Scn.")
    end

    # --- Final Ellipticity Check ---
    e_p_sq = yp0[2]^2 + yp0[3]^2
    e_s_sq = ys0[2]^2 + ys0[3]^2
    if e_p_sq >= 1.0; e_p = sqrt(e_p_sq); error("GeometryPicker Scn=$Scn: Primary orbit non-elliptic (e_p=$e_p)."); end
    if e_s_sq >= 1.0; e_s = sqrt(e_s_sq); error("GeometryPicker Scn=$Scn: Secondary orbit non-elliptic (e_s=$e_s)."); end

    return yp0::SVector{6, Float64}, ys0::SVector{6, Float64}, deltaV_ECI::SVector{3, Float64}
end