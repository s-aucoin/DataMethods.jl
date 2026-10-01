# Library

All functions exported by this package are listed below.


## Library Contents
```@contents
Pages = ["library.md"]
Depth = 3
```

## Index
```@index
```


## Functions
```@docs
IntervalAverage
interp_2d
```

### Calculus Utilities
```@docs
∂_cen
∂_back
int_def_trap
int_trap
```

### Signal Processing
```@docs
Welch
FT_params
detrend2d
dfiltfilt
findspikes
despike
```


### ADCP Utilities
```@docs
get_sig_fnames
```

#### Signal Analysis
```@docs
calculate_SNR
HR_V_variance
HR_c_from_V_variance
PulseLag
IntrinsicVelocityRange
ExtendedVelocityRange
CorrelationFromTurbulence
TurbulenceFromCorrelation
```

#### Wave Analysis
```@docs
wave_U_ofz
```

##### Unexported
```@docs
DataMethods.wave_from_velocity
DataMethods.wave_from_acceleration
DataMethods.unit_spherical_vector_projection_u
DataMethods.unit_spherical_vector_projection_w
DataMethods.rfactor2x
DataMethods.rfactor2z
```

#### Transformations
```@docs
rm_sidelobe_contam
```
##### Coordinate Tranformations
```@docs
align_z0
r_2_slant_r
z_2_brange
brange_2_z
```

##### Velocity Transformations
```@docs
transformation_matrix
heading_R_matrix
pitchroll_R_matrix
rotation_matrix
beam2xyz
beam2ENU
RDI_pitch2pitch
RDI_rotation_matrix
RDI_ENU2XYZ
RDI_XYZ2ENU
```


### CTD Utilities
```@docs
find_ctd_profiles
correct_CT_lag
remove_loops
FindLoops
```


### GPS Utilities
```@docs
add_ms
clean_gps_data
mean_gps_positions
clean_and_mean
index_by_time
```