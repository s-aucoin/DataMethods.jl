# ADCP Utilities

This is a collection of functions that are helpful for analysing ADCP data. Some of them are specifically for pulse-coherent ADCPs, but some are more generally useful.

## General
```@docs; canonical=false
get_sig_fnames
```


## Signal Analysis
These functions relate to analysing the pulse and acoustic signal itself for pulse coherent ADCPs.

```@docs; canonical=false
calculate_SNR
HR_V_variance
HR_c_from_V_variance
PulseLag
IntrinsicVelocityRange
ExtendedVelocityRange
CorrelationFromTurbulence
TurbulenceFromCorrelation
```

## Transformations
```@docs; canonical=false
rm_sidelobe_contam
```

### Coordinate Tranformations

```@docs; canonical=false
align_z0
r_2_slant_r
z_2_brange
brange_2_z
```

### Velocity Transformations
Developed for the Signature 1000, but useful more generally.
```@docs; canonical=false
transformation_matrix
heading_R_matrix
pitchroll_R_matrix
rotation_matrix
beam2xyz
beam2ENU
```

Developed for the RDI Workhorse 600
```@docs; canonical=false
RDI_pitch2pitch
RDI_rotation_matrix
RDI_ENU2XYZ
RDI_XYZ2ENU
```



## Wave Analysis
This function is for wave analysis from the motion of a floating ADCP. The method is described in Aucoin et al. 2026 (under review).

```@docs; canonical=false
wave_U_ofz
```