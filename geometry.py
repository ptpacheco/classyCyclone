#!/usr/bin/env python
# geometry.py:
import numpy as np

## Geometry ID: ###
# Stairmand, D = 100 mm
# Dipleg: D = 37 mm, H = 120 mm
# Bin: D = 200 mm, H = 300 mm
# Units: [mm]

## Input geometry data ###
R_cyc = 50
L_a = 20
H_b = 50
R_vx = 25
wallThickness_vx = 2
R_b = 18.5
H_vx = 50
H_cyl = 150
H_cone = 250
R_dip = 0 # Set to 0 for no dipleg
H_dip = 0 # Set to 0 for no dipleg
R_bin = 0 # Set to 0 for no bin
H_bin = 0 # Set to 0 for no bin

outletLengthToR_ratio = 10
inletLengthToRH_ratio = 6.5

## Input cell size data
coreCellLength = 1
rad_c2c = 1.05
bin_rad_c2c = 1.1
bin_z_c2c = 1.05

## Other geometrical parameters
diagonal_ratio = 0.7 # For disks
side_ratio = 0.62 # For disks
###############################

## Calculated geometry data ###
# Z-levels
H_out = R_vx*outletLengthToR_ratio
H0 = 0
H1 = H0 + H_bin
H2 = H1 + H_dip # if h_dip = 0, h1 = h2
H3 = H2 + H_cone
if H_vx > H_b:
    H4 = H3 + H_cyl - H_vx
    H5 = H3 + H_cyl - H_b
    H6 = H5 + H_b
elif H_vx < H_b:
    H4 = H3 + H_cyl - H_b
    H5 = H3 + H_cyl - H_vx
    H6 = H4 + H_b
else:
    H4 = H3 + H_cyl - H_b
    H5 = H4 + H_b

# Inlet length
L_in = inletLengthToRH_ratio*4*H_b*L_a/(2*(H_b + L_a))

# R-levels
R1 = R_vx
R2 = R1 + wallThickness_vx
R3 = R_cyc
coneBottomContraction = R_b/R_cyc
R1b = R1*coneBottomContraction
R2b = R2*coneBottomContraction
R3b = R_b
if R_dip > R_b: R4b = R_dip
if H_bin != 0: R5b = R_bin

## Angular levels
# Note that angular direction is CCW on X-Y plane, starting from 0Y
Angle_inlet = np.arcsin((R_cyc - L_a)/R_cyc) # angle from 0X to radial position of inlet side wall (towards 0X)
if Angle_inlet > np.pi/4: alpha = np.pi/2 - Angle_inlet
else: alpha = Angle_inlet
Angles = []
if Angle_inlet == np.pi/4:
    for j in [0,1,2,3]:
        Angles.append(j*np.pi/2)
        Angles.append(Angles[-1] + np.pi/4)
else:
    for j in [0,1,2,3]:
        Angles.append(j*np.pi/2)
        Angles.append(Angles[-1] + alpha)
        Angles.append(Angles[-2] + np.pi/4)
        Angles.append(Angles[-3] + np.pi/2 - alpha)

Angles.append(2*np.pi)


## Chopping
z_c2c_cone = (1 - coreCellLength*(1 - coneBottomContraction)/H_cone)
z_size = coreCellLength*2
z_count = np.abs(np.ceil(np.log(coneBottomContraction) / np.log(z_c2c_cone) ))/2
tan_count_core = np.ceil(np.pi*R1*diagonal_ratio/(4*coreCellLength))
rad_count_core = np.ceil(R1*(1 - diagonal_ratio)/coreCellLength)
if rad_count_core%2 == 1:
    rad_count_core += 1
rad_cellN_core = (R1*(1-diagonal_ratio)/rad_count_core) / 2 # controls out of core radial chopping, refers to cell side length at the outermost cell
rad_count_vx = np.ceil(np.log(1 + (R2b - R1b)*(rad_c2c - 1)/(rad_c2c*(rad_cellN_core)) )/np.log(rad_c2c))
rad_count_outer = np.ceil(np.log(1 + (R3b - R2b)*(rad_c2c - 1)/(rad_c2c*rad_cellN_core) )/np.log(rad_c2c)) #- rad_count_vx
if rad_count_outer%2 == 1:
    rad_count_outer += 1
rad_cellN_outer = rad_cellN_core*rad_c2c**(rad_count_outer+1)
if R_dip > R_b:
    rad_count_rim = np.ceil(np.log(1 + (R4b - R3b)*(rad_c2c - 1)/(rad_c2c*rad_cellN_outer) )/np.log(rad_c2c))
    if rad_count_rim%2 == 1:
        rad_count_rim += 1
    rad_cellN_rim = rad_cellN_outer*rad_c2c**(rad_count_rim+1)
    rad_count_rim2 = np.ceil(np.log(1 + (R5b - R4b)*(rad_c2c - 1)/(2*rad_c2c*rad_cellN_rim) )/np.log(rad_c2c))
else:
    rad_count_rim2 = np.ceil(np.log(1 + (R5b - R4b)*(rad_c2c - 1)/(2*rad_c2c*rad_cellN_outer) )/np.log(rad_c2c))


innerSquareTanCount = (tan_count_core*np.sin(alpha))
if innerSquareTanCount%2 == 1:
    innerSquareTanCount += 1

outerSquareTanCount = tan_count_core*(1 - np.sin(alpha))/2
if outerSquareTanCount%2 == 1:
    outererSquareTanCount += 1

inlet_size = rad_cellN_core*3*rad_c2c**(rad_count_outer+1)
