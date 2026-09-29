using NonlinearSolve
import Base.Threads.@threads

export wave_U_ofz


############################################
"""
    wave_U_ofz(IA, IF, IP, z, h, θ, β, α; method=:velocity, incl_trend=false, IMFs2incl=1:size(IA)[2], g = 9.80665)

Calculate the instantaneous surface gravity wave-induced along-beam velocity.

The matrix algebra only works one way so be careful what shape the inputs are.
All the trigonometry and geometry assumes the ADCP is negligably far from the surface. The projections would be off if the distance is significant.

...
# Arguments
- `IA`: instantaneous amplitude from the Hilbert-Huang transform.
- `IF`: instantaneous frequency from the Hilbert-Huang transform.
- `IP`: instantaneous phase from the Hilbert-Huang transform.
- `r`: along-beam range coordinate (from the surface).
- `h`: bottom depth.
- `θ`: angle of the beam relative to vertical toward the beam.
- `θ₀`: fixed beam angle.
- `β`: angle of the beam relative to vertical perpendicular the beam.
- `α`: angle between the wave vector and the beam axis projected onto the surface plane.
- `method`: what property to base the amplitude on. Options are `:velocity` and `:acceleration`.
- `incl_trend=true`: optionally include the mean trend of the IMF decomposition.
- `IMFs2incl`: the mode numbers to include in the calculation.
- `g`: optionally specify the acceleration due to gravity.
- `u0`: initial wavenumber guess for the numerical solver.
...
"""
function wave_U_ofz(IA, IF, IP, r, h, θ, θ₀, β, α; method=:velocity, incl_trend=false, IMFs2incl=1:size(IA)[2], g = 9.80665, u0=0.1)
    # Expects all IMFs and trend

    # Optionally remove the trend
    if !incl_trend
        IA = IA[:,1:end-1]
    end

    ## get the appropriate function for the amplitude input ##
    wave_func_name = Symbol("wave_from_", method)
    wave_func = getfield(DataMethods, wave_func_name)

    ## Get the data type function from IF to propagate it forward ##
    DataTypeFunc = getfield(Base, Symbol(eltype(IF)))

    wave_u_imf = Array{eltype(IA)}(undef, size(r, 1), size(IA,1), length(IMFs2incl))
    @threads for (imfidx, imf_n) in collect(enumerate(IMFs2incl))

        # Extract the amplitude, frequency, and phase of an imf
        A = IA[:,imf_n]
        ωₖ = DataTypeFunc.(2π*IF[:,imf_n])
        ϕₖ = IP[:,imf_n]

        # Find the corresponding wavenumber
        kₖ = Array{eltype(ωₖ)}(undef, size(ωₖ))
        @threads for (pos, freq) in enumerate(ωₖ)
            # Only index h if it is not a constant #
            if length(h) > 1
                h2use = h[pos]
            else
                h2use = h
            end
            disp_r(k, p) = k.*tanh.(k.*h2use) .- freq^2/g
            k_prob = NonlinearProblem(disp_r, u0)
            k_sol = solve(k_prob)
            kₖ[pos] = k_sol.u[1]
        end

        # Calculate the wave velocity
        wave_u_imf[:,:,imfidx] = wave_func(A, ωₖ, ϕₖ, r, h, kₖ, θ, θ₀, β, α)
    end

    return dropdims(sum(wave_u_imf, dims=3), dims=3)
end


"""
    wave_from_velocity(A, ω, ϕ, r, h, k, θ, β, α)

Calculate the projected along-beam velocity from a linear surface gravity wave with amplitude `A` in terms of velocity.

...
# Arguments
- `A`: amplitude in terms of velocity.
- `ω`: angular frequency.
- `ϕ`: wave phase.
- `r`: along-beam range coordinate (from the surface).
- `h`: bottom depth.
- `θ`: angle of the beam relative to vertical toward the beam.
- `θ₀`: fixed beam angle.
- `β`: angle of the beam relative to vertical perpendicular the beam.
- `α`: angle between the wave vector and the beam axis projected onto the surface plane.
...
"""
function wave_from_velocity(A, ω, ϕ, r, h, k, θ, θ₀, β, α)
    amp = A' ./ sinh.(k' .* h')
    u_term = unit_spherical_vector_projection_u.(θ, α, β)' .* cosh.(k' .*(h' .+ r.*rfactor2z.(θ, θ₀, β)')) .* cos.(r*k'.*rfactor2x.(θ, θ₀, β, α)' .- ϕ')
    w_term = -unit_spherical_vector_projection_w.(θ, β)' .* sinh.(k' .*(h' .+ r.*rfactor2z.(θ, θ₀, β)')) .* sin.(r*k'.*rfactor2x.(θ, θ₀, β, α)' .- ϕ')
    
    return amp .* (u_term .+ w_term)
end


"""
    wave_from_acceleration(A, ω, ϕ, r, h, k, θ, β, α)

Calculate the projected along-beam velocity from a linear surface gravity wave with amplitude `A` in terms of acceleration.

...
# Arguments
- `A`: amplitude in terms of acceleration.
- `ω`: angular frequency.
- `ϕ`: wave phase.
- `r`: along-beam range coordinate (from the surface).
- `h`: bottom depth.
- `θ`: angle of the beam relative to vertical toward the beam.
- `θ₀`: fixed beam angle.
- `β`: angle of the beam relative to vertical perpendicular the beam.
- `α`: angle between the wave vector and the beam axis projected onto the surface plane.
...
"""
function wave_from_acceleration(A, ω, ϕ, r, h, k, θ, θ₀, β, α)
    amp = - A' ./ ω' ./ sinh.(k' .* h')
    u_term = unit_spherical_vector_projection_u.(θ, α, β)' .* cosh.(k' .*(h' .+ r.*rfactor2z.(θ, θ₀, β)')) .* cos.(r*k'.*rfactor2x.(θ, θ₀, β, α)' .- (ϕ .+ π/2)')
    w_term = -unit_spherical_vector_projection_w.(θ, β)' .* sinh.(k' .*(h' .+ r.*rfactor2z.(θ, θ₀, β)')) .* sin.(r*k'.*rfactor2x.(θ, θ₀, β, α)' .- (ϕ .+ π/2)')
    
    return amp .* (u_term .+ w_term)
end


"""
    unit_spherical_vector_projection_u(θ, α, β)

Utility function to define the scale factor for projection of the horizontal wave component.
"""
function unit_spherical_vector_projection_u(θ, α, β)
    return sind(θ)*cosd(α)*cosd(β) + cosd(θ)*sind(α)*sind(β)
end


"""
    unit_spherical_vector_projection_w(θ, β)

Utility function to define the scale factor for projection of the vertical wave component.
"""
function unit_spherical_vector_projection_w(θ, β)
    return cosd(θ)*cosd(β)
end

"""
    rfactor2x(θ, θ₀, β, α)

Utility function to define the scale factor for the projection of the ADCP range onto the wave x-axis.
"""
function rfactor2x(θ, θ₀, β, α)
    return secd(θ₀)*unit_spherical_vector_projection_u(θ, α, β)
end

"""
    rfactor2z(θ, θ₀, β)

Utility function to define the scale factor for the projection of the ADCP range onto the wave z-axis.
"""
function rfactor2z(θ, θ₀, β)
    return -secd(θ₀)*unit_spherical_vector_projection_w(θ, β)
end