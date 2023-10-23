#!/usr/bin/env python
import numpy as np

from classy_blocks.mesh import *
from classy_blocks.construct.point import *
from classy_blocks.construct.edges import *
from classy_blocks.construct.flat.face import *
from classy_blocks.construct.operations.loft import Loft
from classy_blocks.construct.operations.extrude import Extrude
from classy_blocks.construct.operations.revolve import Revolve
from classy_blocks.construct.shapes.cylinder import Cylinder
from classy_blocks.construct.shapes.rings import ExtrudedRing
from classy_blocks.construct.shapes.frustum import Frustum
from classy_blocks.base.transforms import *
from classy_blocks.modify.optimizer import Optimizer


########################################################################
##########################  INPUTS  ####################################
########################################################################
## Geometry ID: ###
# Stairmand, D = 100 mm
# Dipleg: D = 37 mm, H = 120 mm
# Bin: D = 200 mm, H = 300 mm

## Units: [mm]
# If you want to use other units, edit the Allrun.sh script and edit the scaling factors in the transformPoints command

## INPUT: Cyclone geometry data ##
R_cyc = 50
L_a = 20
H_b = 50
R_vx = 25
wallThickness_vx = 2
R_b = 18.5
H_vx = 50
H_cyl = 150
H_cone = 250

## INPUT: Dipleg and bin geometry data ##
R_dip = 23 # Set to 0 for no dipleg
H_dip = 120.0 # Set to 0 for no dipleg
R_bin = 100 # Set to 0 for no bin
H_bin = 300.0 # Set to 0 for no bin

## INPUT: Inlet and outlet length to radius ratios ##
outletLengthToR_ratio = 10
inletLengthToRH_ratio = 6.5

## INPUT: Cell grading information ##
core_radSize = 1.8
core_tanSize = core_radSize*1.25
core_axSize = core_radSize*2

dipleg_radC2C = 1.05
bin_radC2C = dipleg_radC2C**2
bin_axC2C = 1.01
########################################################################
####################  AUXILIARY FUNCTIONS  #############################
########################################################################

def endSize_fromStartSize_andC2C(length, start_size, c2c_ratio):
    count = np.log(1 + length / start_size * (1 - (1/c2c_ratio)))/np.log(1/c2c_ratio)
    return start_size*(1/c2c_ratio)**int(count)

########################################################################
##################  RADIAL AND AXIAL LEVELS  ###########################
########################################################################

### CALCULATE DIMENSIONS ###
# Define height levels
Zlvl = [0]
if H_bin > 0: Zlvl.append(Zlvl[-1] + H_bin)
if H_dip > 0: Zlvl.append(Zlvl[-1] + H_dip)
Zlvl.append(Zlvl[-1] + H_cone)
if H_vx > H_b:
    Zlvl.append(Zlvl[-1] + H_cyl - H_vx)
    Zlvl.append(Zlvl[-2] + H_cyl - H_b)
    Zlvl.append(Zlvl[-3] + H_cyl)
elif H_vx < H_b:
    Zlvl.append(Zlvl[-1] + H_cyl - H_b)
    Zlvl.append(Zlvl[-2] + H_cyl - H_vx)
    Zlvl.append(Zlvl[-3] + H_cyl)
else:
    Zlvl.append(Zlvl[-1] + H_cyl - H_b)
    Zlvl.append(Zlvl[-2] + H_cyl)
Zlvl.append(Zlvl[-1] + R_vx * outletLengthToR_ratio)

# Define radial levels on cylinder region
Rlvl = [0]
Rlvl.append(R_vx)
Rlvl.append(Rlvl[-1] + wallThickness_vx)
Rlvl.append(R_cyc)

# Define radial levels on cone bottom region
coneBottomContraction = R_b/R_cyc
RBlvl = []
for i in range(len(Rlvl)): RBlvl.append(Rlvl[i]*coneBottomContraction)
if R_dip > R_b: RBlvl.append(R_dip)
if H_bin != 0: RBlvl.append(R_bin)

# Define angular levels (revolve direction is CCW in X-Y plane)
inletAngle = np.arcsin((R_cyc - L_a)/R_cyc) # angle from 0X to radial position of inlet side wall
subdivAngle = 0.5*(np.pi/2 - inletAngle)
##############################

########################################################################
#####################  MESH CONSTRUCTION  ##############################
########################################################################

### INITIALIZE MESH ###
mesh = Mesh()
##############################

### BUILD BIN, IF IT EXISTS ###
if H_bin != 0:
    bin = [ Cylinder( [0, 0, Zlvl[0]], [0, 0, Zlvl[1]], [0, RBlvl[2], Zlvl[0]] ) ]
    for i in range(3, len(RBlvl)):
        bin.append( ExtrudedRing.expand( bin[-1], RBlvl[i] - RBlvl[i-1] ) )
for i in range(len(bin)-1):
    bin[i].set_end_patch("slave_dipBin")


for i in range(len(bin)):
    mesh.add(bin[i])
##############################

### BUILD MAIN BODY OF CYCLONE ###
# Build dipleg, if it exists
if H_dip != 0:
    dipleg = [ Cylinder( [0, 0, Zlvl[1]], [0, 0, Zlvl[2]], [0, RBlvl[1], Zlvl[1]] ) ]
    for i in range(2, len(RBlvl) - 1):
        dipleg.append( ExtrudedRing.expand( dipleg[-1], RBlvl[i] - RBlvl[i-1] ) )
for i in range(len(dipleg)):
    dipleg[i].set_start_patch("master_dipBin")
    mesh.add(dipleg[i])

# Build cone
cone = [ Frustum( [0, 0, Zlvl[2]], [0, 0, Zlvl[3]], [RBlvl[1], 0, Zlvl[2]], Rlvl[1] ) ]
for i in range(2, len(Rlvl)):
    revFace = Face( [[RBlvl[i-1], 0, Zlvl[2]] , [Rlvl[i-1], 0, Zlvl[3]] , [Rlvl[i], 0, Zlvl[3]] , [RBlvl[i], 0, Zlvl[2]]] )
    for j in range(8):
        cone.append( Revolve( revFace.copy().rotate(j*np.pi/4, [0,0,1], [0,0,0]), np.pi/4, [0,0,1], [0,0,0] ) )
for i in range(len(cone)):
    mesh.add(cone[i])
if H_bin != 0 and H_dip == 0: # If there is a bin but no dipleg, set master patches for merging
    cone[0].set_start_patch("master_dipBin")
    for i in range(len(cone)-1):
        cone[i+1].set_patch("right","master_dipBin")

# Build bottom section of cylindrical region
cyl1 = [ Cylinder( [0, 0, Zlvl[3]], [0, 0, Zlvl[4]], [0, Rlvl[1], Zlvl[3]] ) ]
for i in range(2, len(Rlvl)):
    cyl1.append( ExtrudedRing.expand( cyl1[-1], Rlvl[i] - Rlvl[i-1] ) )
for i in range(len(cyl1)):
    mesh.add(cyl1[i])

# If H_b != H_vx, build extra cyl section
if H_vx != H_b:
    cyl2 = []
    if H_b > H_vx:
        cyl2.append( Cylinder([0, 0, Zlvl[4]], [0, 0, Zlvl[5]], [0, Rlvl[1], Zlvl[4]]) )
        for i in range(2, len(Rlvl)):
            revFace = Face( [[Rlvl[i-1], 0, Zlvl[4]] , [Rlvl[i-1], 0, Zlvl[5]] , [Rlvl[i], 0, Zlvl[5]] , [Rlvl[i], 0, Zlvl[4]]] )
            if i == 2: a = 8
            else: a = 6
            for j in range(a):
                cyl2.append( Revolve( revFace.copy().rotate((j+2)*np.pi/4, [0,0,1], [0,0,0]), np.pi/4, [0,0,1], [0,0,0] ) )
    else:
        for i in range(3,len(Rlvl)):
            revFace = Face( [[Rlvl[i-1], 0, Zlvl[4]] , [Rlvl[i-1], 0, Zlvl[5]] , [Rlvl[i], 0, Zlvl[5]] , [Rlvl[i], 0, Zlvl[4]]] )

            for j in range(2,8):
                cyl2.append( Revolve( revFace.copy().rotate(j*np.pi/4, [0,0,1], [0,0,0]), np.pi/4, [0,0,1], [0,0,0] ) )
    for i in range(len(cyl2)):
        mesh.add(cyl2[i])

# Build topmost cyl section
cylF = []
revFace = Face( [[Rlvl[-2], 0, Zlvl[-3]] , [Rlvl[-2], 0, Zlvl[-2]] , [Rlvl[-1], 0, Zlvl[-2]] , [Rlvl[-1], 0, Zlvl[-3]]] )
for i in range(6):
    cylF.append( Revolve( revFace.copy().rotate((i+2)*np.pi/4, [0,0,1], [0,0,0]), np.pi/4, [0,0,1], [0,0,0] ) )
for i in range(len(cylF)):
    mesh.add(cylF[i])
cylF[0].set_patch("bottom", "master_topStitch2")
cylF[-1].set_patch("top","master_topStitch1")
# Build outlet
if H_vx > H_b:
    outlet = Cylinder( [0, 0, Zlvl[-4]], [0, 0, Zlvl[-1]], [R_vx, 0, Zlvl[-4]] )
else:
    outlet = Cylinder( [0, 0, Zlvl[-3]], [0, 0, Zlvl[-1]], [R_vx, 0, Zlvl[-3]] )
outlet.set_end_patch("outlet")
mesh.add(outlet)
##############################

### BUILD INLET AND MERGE REGION ###
P_bottomIn = Point([Rlvl[2] , 0, Zlvl[-3]])
P_topIn = Point([Rlvl[2] , 0, Zlvl[-2]])
P_topOut = Point([Rlvl[3] , 0, Zlvl[-2]])
P_bottomOut = Point([Rlvl[3] , 0, Zlvl[-3]])
P_subdivInBottom = Point( [np.cos(inletAngle + subdivAngle)*Rlvl[2], np.sin(inletAngle + subdivAngle)*Rlvl[2], Zlvl[-3]] )
P_subdivInTop = Point( [np.cos(inletAngle + subdivAngle)*Rlvl[2], np.sin(inletAngle + subdivAngle)*Rlvl[2], Zlvl[-2]] )
P_subdivBottom = Point([ np.sin(subdivAngle) * Rlvl[3]/np.cos(subdivAngle), Rlvl[3], Zlvl[-3]])
P_subdivTop = Point([np.sin(subdivAngle) * Rlvl[3]/np.cos(subdivAngle), Rlvl[3], Zlvl[-2]])

face0X = Face( [P_bottomIn.position, P_topIn.position, P_topOut.position, P_bottomOut.position] )
faceSD1 = face0X.copy().rotate(inletAngle, [0,0,1], [0,0,0])
faceSD2 = Face( [P_subdivInBottom.position, P_subdivInTop.position, P_subdivTop.position, P_subdivBottom.position])
face0Y = face0X.copy().rotate(np.pi/2, [0,0,1], [0,0,0])

inlet = []
inlet.append( Revolve(face0X, inletAngle, [0,0,1], [0,0,0]) )
inlet[-1].set_patch("bottom", "slave_topStitch1")
inlet.append( Loft(faceSD1, faceSD2) )
inlet[-1].add_side_edge(0, Origin([0,0,Zlvl[-3]]))
inlet[-1].add_side_edge(1, Origin([0,0,Zlvl[-2]]))
inlet.append( Loft(faceSD2, face0Y) )
inlet[-1].add_side_edge(0, Origin([0,0,Zlvl[-3]]))
inlet[-1].add_side_edge(1, Origin([0,0,Zlvl[-2]]))
inlet[-1].set_patch("top", "slave_topStitch2")

if H_b > H_vx:
    for i in [0,1,2]:
        extrFace = inlet[i].get_face("left").copy().translate([0,0, -(H_b-H_vx)])
        if i == 0: extrFace.add_edge(2, Origin([0,0, Zlvl[-4]]))
        extrFace.add_edge(0, Origin([0,0, Zlvl[-4]]))
        inlet.append( Extrude( extrFace, [0, 0, H_b-H_vx ] ) )

inletLength = inletLengthToRH_ratio * 2*(L_a*H_b)/(L_a + H_b)

P_inletTopOut = [inletLength, Rlvl[-1], Zlvl[-2]]
P_inletTopIn = [inletLength, Rlvl[-1] - L_a, Zlvl[-2]]
P_inletBottomOut = [inletLength, Rlvl[-1], Zlvl[-3]]
P_inletBottomIn = [inletLength, Rlvl[-1] - L_a, Zlvl[-3]]

if H_b > H_vx:
    for i in [1, 4]:
        loftFace = inlet[i].get_face("back").copy() # FIX THIS TO GET IT WORKING FOR OTHER CASES ###################################
        #inlet[i].set_patch("back", "slave_inlet")
else:
    loftFace = inlet[1].get_face("back").copy()
    #inlet[1].set_patch("back", "slave_inlet")
    loftFace2 = Face( [P_inletBottomOut, P_inletTopOut, P_inletTopIn, P_inletBottomIn ] )
    inletExtr = Loft(loftFace, loftFace2)

inletExtr.set_patch("top", "inlet")

inletExtr.set_patch("bottom", "master_inletStitch")
inlet[1].set_patch("back", "slave_inletStitch")
if H_b <= H_vx:
    for i in range(len(inlet)):
        inlet[i].set_patch("left","master_bottomIn")
else:
    for i in range(3, len(inlet)):
        inlet[i].set_patch("bottom","master_bottomIn")

for i in range(len(inlet)):
    mesh.add(inlet[i])
mesh.add(inletExtr)
##############################

########################################################################
#########################  GRADINGS  ###################################
########################################################################
### DEFINE GRADINGS IN MAIN BODY ###
tanCount = np.ceil((np.pi/4)*Rlvl[1]/core_tanSize) # for quarter circle
radCount1 = np.ceil(Rlvl[1]/(2*core_radSize)) # Core
radCount2 = np.ceil((Rlvl[2]-Rlvl[1])/core_radSize) # VX
radCount3 = np.ceil((Rlvl[3]-Rlvl[2])/core_radSize) # Outer
radCount4 = np.ceil((RBlvl[4] - RBlvl[3])/(coneBottomContraction*core_radSize)) # Dipleg rim
bin_radSize = endSize_fromStartSize_andC2C(RBlvl[4] - RBlvl[3], 2*coneBottomContraction*core_radSize, dipleg_radC2C)
if tanCount % 2 != 0: tanCount += 1 # Make counts divisible by 2 for better poly faces in dipleg - bin transition
if radCount1 % 2 != 0: radCount1 += 1
if radCount3 % 2 != 0: radCount3 += 1
if radCount4 % 2 != 0: radCount4 += 1
cone_axialC2C = (1 - core_axSize*(1 - coneBottomContraction)/H_cone)

outlet.chop_tangential(count = tanCount) # Defines tan chop for all cyclone except bin and cylF
outlet.chop_radial(count = radCount1)#core_radSize*(1/0.7)) # Defines rad chop for core (outlet to dip)
outlet.chop_axial(start_size = core_axSize) # Defines axial chop for outlet

cyl1[0].chop_axial(start_size = core_axSize) # Defines axial chop for cyl1
cyl1[1].chop_radial(count = radCount2) # Defines rad chop for VX (outlet to dip)
cyl1[2].chop_radial(count = radCount3) # Defines rad chop for outer radius (outlet to dip)

cone[0].chop_axial(start_size = core_axSize, c2c_expansion = cone_axialC2C, invert = True) # Defines axial chop for cone

cone_axialEndSize = endSize_fromStartSize_andC2C(H_cone, core_axSize, cone_axialC2C)
if H_dip != 0:
    dipleg[0].chop_axial(start_size = cone_axialEndSize) # Defines axial chop for dipleg
    if R_dip > R_b: dipleg[3].chop_radial(count = radCount4, c2c_expansion = dipleg_radC2C) #int((RBlvl[3]-RBlvl[2])/core_radSize)) # Defines rad chop for dipleg rim

if H_vx < H_b: # Define axial chop for cyl2, if it exists
    cyl2[0].chop_axial(start_size = core_axSize)
elif H_vx > H_b:
    cyl2[0].chop(0, start_size = core_axSize)

cylF[0].chop(0, start_size = core_axSize)

### DEFINE GRADINGS IN BIN ###
if H_bin != 0:
    # Axial chopping
    bin[0].chop_axial(start_size = cone_axialEndSize*2, c2c_expansion = bin_axC2C, invert = True)
    if R_dip > R_b: bin[-1].chop_axial(start_size = cone_axialEndSize*2, c2c_expansion = bin_axC2C, invert = True)

    # Tangential chopping
    bin[0].chop_tangential(count = tanCount/2)

    # Radial chopping
    bin[0].chop_radial(count = np.ceil(radCount1/2 + radCount2/2))
    bin[1].chop_radial(count = radCount3/2)
    bin[2].chop_radial(count = radCount4/2, c2c_expansion = bin_radC2C)
    if R_dip > R_b: bin[3].chop_radial(start_size = bin_radSize, c2c_expansion = bin_radC2C) # c2c_expansion is not a typo

### DEFINE GRADINGS IN INLET MERGE REGION ###
tanCountIn1 = np.ceil(inletAngle*Rlvl[1]/(core_tanSize))
tanCountIn2 = np.ceil((subdivAngle)*Rlvl[1]/(core_tanSize))
if tanCountIn1 % 2 != 0: tanCountIn1 += 1 # Make counts divisible by 2 to match rest of mesh
if tanCountIn2 % 2 != 0: tanCountIn2 += 1
if H_b > H_vx:
    #inlet[0].chop(0, start_size = core_axSize)
    #inlet[-1].chop(2, start_size = core_axSize) # Axial
    #inlet[0].chop(1, count = radCount3) # Radial
    inlet[0].chop(2, count = tanCountIn1) # Tangential
    inlet[1].chop(2, count = tanCountIn2)
    inlet[2].chop(2, count = tanCountIn2)
else:
    inlet[-1].chop(0, start_size = core_axSize) # Axial
    #inlet[-1].chop(1, count = radCount3) # Radial
    inlet[-3].chop(2, count = tanCountIn1) # Tangential
    inlet[-2].chop(2, count = tanCountIn2)
    inlet[-1].chop(2, count = tanCountIn2)

inletExtr.chop(0, start_size = core_axSize)
inletExtr.chop(2, start_size = core_radSize, c2c_expansion = 1.02)
##############################

### MERGE PATCHES ###
if H_bin != 0: mesh.merge_patches("slave_dipBin", "master_dipBin")

for i in [-1, -2]:
    cyl1[2].shell[i].set_patch("top", "slave_bottomIn")
mesh.merge_patches("slave_bottomIn", "master_bottomIn")
#mesh.merge_patches("master_inlet", "slave_inlet")
##############################

### OPTIMIZE AND WRITE MESH ###
#mesh.assemble()
mesh.set_default_patch("walls", "wall")
#optimizer = Optimizer(mesh) # some problems in topmost cyl region? dunno why
#optimizer.optimize()
mesh.write('bodyMesh/system/blockMeshDict', debug_path="./body.vtk")
