export calculate_SNR,
        HR_V_variance, HR_c_from_V_variance,
        PulseLag, IntrinsicVelocityRange, ExtendedVelocityRange,
        CorrelationFromTurbulence, TurbulenceFromCorrelation


############################################
"""
    calculate_SNR(As, An; mode=:log)

Calculate the signal-to-noise ratio of signal amplitude `As` relative to noise amplitude `An`.

Specify that `As` and `An` are in `:log`-space or `:lin`-space with `mode`.
"""
function calculate_SNR(As, An; mode=:log)
    if mode == :log
        return As - An

    elseif mode == :lin
        return 20*log10(As/An)
    end
end


####################
"""
    HR_V_variance(c, ds)

Calculate the velocity variance of sample of correlation `c` and parameters in NCDataset `ds`.
"""
function HR_V_variance(c, ds)

    Δr = ds.group["Config"].attrib["bursthr_cellSize"]
    cs = Array(ds.group["Data"].group["IBurstHR"]["SpeedOfSound"]) # (m/s) Speed of sound
    F₀ = 1000*ds.group["Config"].attrib["beamConfiguration5_frequency"] # (Hz) Center frequency
    ΔF = 1000*ds.group["Config"].attrib["beamConfiguration5_bandwidth"] # (Hz) Bandwidth
    τ = 2 * ds.group["Config"].attrib["bursthr_lag"] ./ cs # (s) lag between pings in ping pair

    return cs'.^3 .* (c.^(-2) .- 1) ./ (64 * π^2 * Δr * ΔF * F₀.^2 .* τ'.^2)
end

"""
    HR_V_variance(c, Δr, F₀, ΔF, τ; cs = 1500)

Calculate the velocity variance of sample of correlation `c`.

...
# Arguments
- `Δr`: cell size.
- `F₀`: center frequency.
- `ΔF`: frequency bandwidth.
- `τ`: time lag between pings in a ping pair.
- `cs`: speed of sound.
...
"""
function HR_V_variance(c, Δr, F₀, ΔF, τ; cs = 1500)
    return cs^3 * (c^(-2) - 1) / (64 * π^2 * Δr * ΔF * F₀^2 * τ^2)
end


"""
    HR_c_from_V_variance(var, ds)

Calculate the velocity correlation of sample of variance `var` and parameters in `ds`.
"""
function HR_c_from_V_variance(var, ds)

    Δr = ds.group["Config"].attrib["bursthr_cellSize"]
    cs = Array(ds.group["Data"].group["IBurstHR"]["SpeedOfSound"]) # (m/s) Speed of sound
    F₀ = 1000*ds.group["Config"].attrib["beamConfiguration5_frequency"] # (Hz) Center frequency
    ΔF = 1000*ds.group["Config"].attrib["beamConfiguration5_bandwidth"] # (Hz) Bandwidth
    τ = 2 * ds.group["Config"].attrib["bursthr_lag"] ./ cs # (s) lag between pings in ping pair

    return (var .* (64 * π^2 * Δr * ΔF * F₀.^2 * τ'.^2) ./ (cs'.^3) .+ 1).^(-1/2)
end


"""
    HR_c_from_V_variance(var, Δr, F₀, ΔF, τ; cs = 1500)

Calculate the velocity correlation of sample of variance `var`.

...
# Arguments
- `Δr`: cell size.
- `F₀`: center frequency.
- `ΔF`: frequency bandwidth.
- `τ`: time lag between pings in a ping pair.
- `cs`: speed of sound.
...
"""
function HR_c_from_V_variance(var, Δr, F₀, ΔF, τ; cs = 1500)
    return (var * (64 * π^2 * Δr * ΔF * F₀^2 * τ^2) / (cs.^3) + 1)^(-1/2)
end


####################
"""
    PulseLag(λ, c)

Calculate the pulse lag `τ` from lag distance `λ` and speed of sound `c`.
"""
function PulseLag(λ, c)
    return 2 * λ / c
end

"""
    IntrinsicVelocityRange(c, F₀, τ)

Calculate the intrinsic velocity range for a pulse-coherent ADCP with speed of sound `c`, carrier frequency `F₀`, and pulse lag `τ`.
"""
function IntrinsicVelocityRange(c, F₀, τ)
    return c / (4 * F₀ * τ)
end

"""
    ExtendedVelocityRange(c, F₀, τ; nmax=6)

Calculate the extended velocity range for a pulse-coherent ADCP with speed of sound `c`, carrier frequency `F₀`, 
pulse lag `τ`, and the number of side correlations considered `nmax`.

"""
function ExtendedVelocityRange(c, F₀, τ; nmax=6)
    return 2 * nmax * IntrinsicVelocityRange(c, F₀, τ)
end


####################
"""
    CorrelationFromTurbulence(ϵ, d, Vᵣ; C₁=0.707)

Calculate the maximum correlation for a pulse-coherent ADCP expected in a region with TKE dissipation rate `ϵ`.

...
# Arguments
- `ϵ`: TKE dissipation rate (m² s⁻³)
- `d`: characteristic measurement size (largest of either cell length or width) (m)
- `Vᵣ`: intrinsic velocity range (m s⁻¹).
...
"""
function CorrelationFromTurbulence(ϵ, d, Vᵣ; C₁=0.707)
    factor = -1/3 * (2π)^(-2/3) * C₁
    return exp(factor * π^2 * Vᵣ^(-2) * (ϵ*d)^(2/3))
end


"""
    TurbulenceFromCorrelation(cor, d, Vᵣ; C₁=0.707)

Calculate the TKE dissipation rate that would result in correlation `c` for a pulse-coherent ADCP.

...
# Arguments
- `cor`: correlation
- `d`: characteristic measurement size (largest of either cell length or width) (m)
- `Vᵣ`: intrinsic velocity range (m s⁻¹).
...
"""
function TurbulenceFromCorrelation(cor, d, Vᵣ; C₁=0.707)
    factor = -1/3 * (2π)^(-2/3) * C₁
    return (log(cor) / (factor * π^2) * Vᵣ^(2))^(3/2) / d
end