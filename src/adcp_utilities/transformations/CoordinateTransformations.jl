export align_z0, r_2_slant_r, z_2_brange, brange_2_z


############################################
"""
    align_z0(data, zs, z0, Δz; buffer=15)

Shift the columns in `data` and `zs` such that the point closest to `z0` in `zs` is aligned in a row.

This is useful for shifting bottom-referenced data to surface-referenced and vice versa, for example.

...
# Arguments
- `data`: the data to shift.
- `zs`: the coordinates of `data` (must be size(data)).
- `z0`: the value in `zs` to use for aligning.
- `Δz`: the range step to use to construct a reference coordinate vector.
- `buffer`: an arbitrary buffer to add to the size shifted array to ensure all the data fits (must be greater than the largest shift). The arrays are reduced to the correct size after.
...
"""
function align_z0(data, zs, z0, Δz; buffer=50)

    # Create the reference z axis, bigger than it needs to be to adjust later #
    half_idx = round(Int, (size(data)[1] + buffer) / 2)        # To make sure the array is bigger than the original data
    ref_z = collect(z0 - half_idx*Δz:Δz:z0 + half_idx*Δz)    # Evenly spaced z axis to reference to

    z0_idx = half_idx + 1                                    # Index of z0 in ref_z

    # Find the index of z closest to z0 #
    (vals, idx) = findmin(abs.(zs .- z0), dims=1)
    idx = vec(first.(Tuple.(idx)))

    # Shift each column of the matrix so that z0 is aligned in a row
    shift_zs = fill!(Array{eltype(zs)}(undef, size(ref_z)[1], size(zs)[2]), NaN)
    shift_data = fill!(Array{eltype(data)}(undef, size(ref_z)[1], size(zs)[2]), NaN)
    @threads for tt in 1:size(data)[2]
        idx_start = z0_idx - (idx[tt] - 1)
        idx_end = z0_idx + (size(zs)[1] - idx[tt])

        shift_zs[idx_start:idx_end, tt] = zs[:,tt]
        shift_data[idx_start:idx_end, tt] = data[:,tt]
    end

    # Determine which rows contain useful information by the number of nans #
    usefull_row = first.(Tuple.(findall(sum(isnan.(shift_zs), dims=2) .< size(shift_zs)[2])))

    # Resample the data to remove useless points #
    return (ref_z = ref_z[usefull_row], z_true = shift_zs[usefull_row,:], data = shift_data[usefull_row,:])
end


"""
    r_2_slant_r(range, θ)

Convert ADCP nominal `range` or cell size to the equivalent for slant beams at angle `θ` (in deg) to beam 5.
"""
function r_2_slant_r(range, θ)
    return range ./ cosd(θ)
end


"""
    z_2_brange(z, θ, ϕ, roll, pitch)

Convert `z`-coordinate to along beam range using the ADCP `roll` and `pitch` attitude.

`θ` is the slant beam angle in the ADCP x-direction from z and `ϕ` is the same in the y-direction.
"""
function z_2_brange(z, θ, ϕ, roll, pitch)
    return z .* sqrt.(1 .+ tand.(ϕ .+ roll).^2 .+ tand.(θ .+ pitch).^2)'
end


"""
    brange_2_z(brange, θ, ϕ, roll, pitch)

Convert along-beam range to z coordinate using the ADCP `roll` and `pitch` attitude.

`θ` is the slant beam angle in the ADCP x-direction from z and `ϕ` is the same in the y-direction.
"""
function brange_2_z(brange, θ, ϕ, roll, pitch)
    return brange ./ sqrt.(1 .+ tand.(ϕ .+ roll).^2 .+ tand.(θ .+ pitch).^2)'
end