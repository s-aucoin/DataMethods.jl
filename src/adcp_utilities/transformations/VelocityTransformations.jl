export transformation_matrix, heading_R_matrix, pitchroll_R_matrix, rotation_matrix, 
        beam2xyz, beam2ENU,
        RDI_pitch2pitch, RDI_rotation_matrix,
        RDI_ENU2XYZ, RDI_XYZ2ENU


# missing functions to go from ENU to XYZ or beam, and from beam to XYZ #
############################################
"""
    transformation_matrix(θ)

Make the transformation matrix for a 4-beam ADCP with beam angle `θ` in degrees.
"""
function transformation_matrix(θ)
    ϕ = π * θ / 180 # convert from degrees to radians
    # Define the transformation matrix #
    return 1/2 * [csc(ϕ) 0 -csc(ϕ) 0;
                    0 -csc(ϕ) 0 csc(ϕ);
                    sec(ϕ) 0 sec(ϕ) 0;
                    0 sec(ϕ) 0 sec(ϕ)]
end


"""
    heading_R_matrix(heading)

Make the heading rotation matrix for a 4-beam ADCP with magnetometer reading `heading`.
"""
function heading_R_matrix(heading)
    ntimes = size(heading)
    return cat( [cos(heading) sin(heading) zeros(ntimes)],
                [-sin(heading) cos(heading) zeros(ntimes)],
                [zeros(ntimes) zeros(ntimes) ones(ntimes)], dims=1)
end


"""
    pitchroll_R_matrix(roll, pitch)

Make the pitch/roll rotation matrix for a 4-beam ADCP with attitude readings `roll` and `pitch`.
"""
function pitchroll_R_matrix(roll, pitch)
    ntimes = size(roll)
    return cat( [cos(pitch)  -sin(pitch)*sin.(roll)  -cos.(roll)*sin.(pitch)],
                [zeros(ntimes)      cos(roll)              -sin(roll)],
                [sin(pitch)  sin.(roll)*cos.(pitch)  cos(pitch)*cos.(roll)], dims=1)
end


"""
    rotation_matrix(AHRSraw)

Make the full rotation matrix for a 4-beam ADCP using the data from an on-board AHRS sensor.
"""
function rotation_matrix(AHRSraw)

    ntimes = size(AHRSraw)[2]
    R = Array{Float64}(undef, 4, 4, ntimes) # the matrix needs to be modified for 4 beams
    R[1:3, 1:3, :] = permutedims(reshape(AHRSraw, (3,3,ntimes)), [2,1,3]) # copy over the 3x3 section

    # The modification is mostly copying rows/columns #
    R[4,1:4,:] .= R[3,1:4,:]             # copy the third row into the fourth row
    R[1:4,4,:] .= R[1:4,3,:]             # Same for column
    R[3,4,:] .= 0; R[4,3,:] .= 0         # replace the (3,4) and (4,3) positions with 0

    return R
end


"""
    rotation_matrix(heading, pitch, roll)

Make the full rotation matrix for a 4-beam ADCP using the attitude data (in degrees) `heading`, `pitch`, and `roll`.
"""
function rotation_matrix(heading, pitch, roll)
    # convert to radians first #
    heading_rad = π * (heading .- 90 )/180
    # raw heading is defined as 0 for x facing north, so subtracting 90 makes 0 defined facing y
    pitch_rad = π * pitch/180
    roll_rad = π * roll/180

    # Multiply the heading and pitch/roll matrices to get the full rotation matrix #
    R3x3 = reshape(reduce(hcat, heading_R_matrix.(heading_rad) .* pitchroll_R_matrix.(roll_rad, pitch_rad)), 3, 3, :)

    ## Need to reshape the matrix for the 4 beams ##
    ntimes = size(R3x3)[3]
    R = Array{Float64}(undef, 4, 4, ntimes) # the matrix needs to be modified for 4 beams
    R[1:3, 1:3, :] = R3x3 # copy over the 3x3 section

    # The modification is mostly copying rows/columns #
    R[4,1:4,:] .= R[3,1:4,:]             # copy the third row into the fourth row
    R[1:4,4,:] .= R[1:4,3,:]             # Same for column
    R[3,4,:] .= 0; R[4,3,:] .= 0         # replace the (3,4) and (4,3) positions with 0

    return R
end



#######################
## 4-beam ADCP transformations (signature 1000) ##
"""
    beam2xyz(V1, V2, V3, V4; T=nothing, θ=nothing)

Convert 4 beam-coordinate velocities `V1`, `V2`, `V3`, `V4` to xyz-referenced velocities.
Uses either a supplied transformation matrix `T`, or makes one from the beam angle `θ`.
"""
function beam2xyz(V1, V2, V3, V4; T=nothing, θ=nothing)
    if T == nothing .&& θ == nothing
        error("Must provide transformation matrix T or beam angle θ")

    elseif T == nothing .&& θ != nothing
        # Make the transformation matrix #
        T = transformation_matrix(θ)
    end

    V = [V1, V2, V3, V4] # put the beam velocities into a length 4 vector of each beam
    xyz = T * V

    return (u = xyz[1], v = xyz[2], w1 = xyz[3], w2 = xyz[4])
end


"""
    beam2ENU(V1, V2, V3, V4, AHRSraw; T=nothing, θ=nothing)

Convert 4 beam-coordinate velocities `V1`, `V2`, `V3`, `V4` to earth-referenced velocities using the data from an AHRS sensor.
Uses either a supplied transformation matrix `T`, or makes one from the beam angle `θ`.
"""
function beam2ENU(V1, V2, V3, V4, AHRSraw; T=nothing, θ=nothing)
    if T == nothing .&& θ == nothing
        error("Must provide transformation matrix T or beam angle θ")

    elseif T == nothing .&& θ != nothing
        # Make the transformation matrix #
        T = transformation_matrix(θ)
    end

    # Make the rotation matrix #
    R = rotation_matrix(AHRSraw)


    ENU = Array{eltype(V1)}(undef, size(V1)..., 4)
    @threads for (i, j) in collect(Iterators.product(axes(V1)[1], axes(V1)[2]))
        ENU[i, j, :] = R[:,:,j] * T * [V1[i,j]; V2[i,j]; V3[i,j]; V4[i,j]]
    end

    return (uE = ENU[:,:,1], vN = ENU[:,:,2], U1 = ENU[:,:,3], U2 = ENU[:,:,4])
end


"""
    beam2ENU(V1, V2, V3, V4, heading, pitch, roll; T=nothing, θ=nothing)

Convert 4 beam-coordinate velocities `V1`, `V2`, `V3`, `V4` to earth-referenced velocities using the attitude data.
Uses either a supplied transformation matrix `T`, or makes one from the beam angle `θ`.
"""
function beam2ENU(V1, V2, V3, V4, heading, pitch, roll; T=nothing, θ=nothing)
    if T == nothing .&& θ == nothing
        error("Must provide transformation matrix T or beam angle θ")

    elseif T == nothing .&& θ != nothing
        # Make the transformation matrix #
        T = transformation_matrix(θ)
    end

    # Make the rotation matrix #
    R = rotation_matrix(heading, pitch, roll)


    ENU = Array{eltype(V1)}(undef, size(V1)..., 4)
    for (i, j) in collect(Iterators.product(axes(V1)[1], axes(V1)[2]))
        ENU[i, j, :] = R[:,:,j] * T * [V1[i,j]; V2[i,j]; V3[i,j]; V4[i,j]]
    end

    return (uE = ENU[:,:,1], vN = ENU[:,:,2], U1 = ENU[:,:,3], U2 = ENU[:,:,4])
end



#######################
## RDI ADCP transformations ##

"""
    RDI_pitch2pitch(pitch, roll)

Convert the RDI instrument measured `pitch` to the gimbal pitch required for the rotation matrix.

Units are in degrees.
"""
function RDI_pitch2pitch(pitch, roll)
    return atand(tand(pitch)*cosd(roll))
end



"""
    RDI_rotation_matrix(pitch, roll, heading)

Make the rotation matrix for an RDI instrument defined with `pitch`, `roll`, and `heading`.

Units are in degrees.
"""
function RDI_rotation_matrix(pitch, roll, heading)
    return cat( [ (cosd(heading)*cosd(roll) + sind(heading)*sind(pitch)*sind(roll))  (sind(heading)*cosd(pitch))  (cosd(heading)*sind(roll) - sind(heading)*sind(pitch)*cosd(roll))],
                [(-sind(heading)*cosd(roll) + cosd(heading)*sind(pitch)*sind(roll))  (cosd(heading)*cosd(pitch))  (-sin(heading)*sind(roll) - cosd(heading)*sind(pitch)*cosd(roll))],
                [              (-cosd(pitch)*sind(roll))                                     (sind(pitch))                         (cosd(pitch)*cos(roll))], dims=1)
end



"""
    RDI_ENU2XYZ(uE, vN, U, pitch, roll, heading)

Convert earth-referenced to instrument coordinate velocities for RDI data.

Angle units are in degrees.
"""
function RDI_ENU2XYZ(uE, vN, U, pitch, roll, heading)

    # Make the inverse rotation matrix #
    R⁻¹ = inv.(RDI_rotation_matrix.(pitch, roll, heading))

    XYZ = Array{eltype(uE)}(undef, size(uE)..., 3)
    @threads for (i, j) in collect(Iterators.product(axes(uE)[1], axes(uE)[2]))
        XYZ[i, j, :] = R⁻¹[j] * [uE[i,j]; vN[i,j]; U[i,j]]
    end

    return (u = XYZ[:,:,1], v = XYZ[:,:,2], w = XYZ[:,:,3])

end


"""
    RDI_XYZ2ENU(u, v, w, pitch, roll, heading)

Convert instrument coordinates to earth-referenced for RDI data.

Angle units are in degrees.
"""
function RDI_XYZ2ENU(u, v, w, pitch, roll, heading)

    # Make the rotation matrix #
    R = RDI_rotation_matrix.(pitch, roll, heading)

    ENU = Array{eltype(u)}(undef, size(u)..., 3)
    @threads for (i, j) in collect(Iterators.product(axes(u)[1], axes(u)[2]))
        ENU[i, j, :] = R[j] * [u[i,j]; v[i,j]; w[i,j]]
    end

    return (uE = ENU[:,:,1], vN = ENU[:,:,2], U = ENU[:,:,3])

end