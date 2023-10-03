#!/usr/bin/env python
from typing import List
import numpy as np

from classy_blocks.mesh import *
from classy_blocks.construct.point import *
from classy_blocks.construct.edges import *
from classy_blocks.construct.flat.face import Face
from classy_blocks.construct.flat.sketches.disk import HalfDisk, QuarterDisk
from classy_blocks.construct.shapes.shape import SketchedShape
from classy_blocks.construct.operations.loft import Loft
from classy_blocks.construct.operations.extrude import Extrude
from classy_blocks.construct.operations.revolve import Revolve
from classy_blocks.base.transforms import *

from importlib import resources as impresources
import geometry
from geometry import *

print("Making blockMeshDicts...")

# Initialize mesh
mesh = Mesh()
if H_bin != 0:
    binmesh = Mesh()

# Build core in conic region, along with dipleg ###
if Angle_inlet == np.pi/4:
    # In case inlet angle is exactly equal to 45deg, normal H-grid is fine
    # But I don't have time for this right now so here's a placeholder:
    print("Inlet angle equal to pi/4, aborting.")
    exit()

    #coreSketch_QIN_top = QuarterDisk([0, 0, H3], [R1, 0, H3], [0, 0, 1])
    #coreSketch_QIN_bottom = coreSketch_QIN_top.copy().scale(ratio = coneBottomContraction).translate(displacement = [0,0,-H_cone])
    #core_QIn = []
    #for i in range(len(coreSketch_QIN_bottom.faces)):
    #    core_QIn.append(Loft(coreSketch_QIN_bottom.faces[i], coreSketch_QIN_top[i]))
    #    mesh.add(core_QIn[i])
else:
    # Define quarter cylinder of core (same quarter where inlet is)
    # Inner circle is subdivided into four sections
    # First, points:
    P0 = Point([0, 0, H3])
    P1 = Point([R1*diagonal_ratio,0,H3]) #side_ratio, 0, H3])
    P2 = Point([0, R1*diagonal_ratio, H3]) #side_ratio, H3])
    P3 = Point([R1, 0, H3])
    P4 = Point([0, R1, H3])
    D1 = Point([diagonal_ratio*R1*np.sin(alpha), diagonal_ratio*R1*np.sin(alpha), H3])
    D2 = Point([R1*diagonal_ratio, 0, H3]).rotate(np.pi/4, [0,0,1], [0,0,0])
    D3 = P3.copy().rotate(np.pi/4, [0,0,1], [0,0,0])
    SH1 = D1.copy().rotate(-np.pi/4, [0,0,1], [0,0,0])#Point([diagonal_ratio*R1*np.sin(Angle_inlet), 0, H3])
    SH2 = D1.copy().rotate(np.pi/4, [0,0,1], [0,0,0])#Point([0, diagonal_ratio*R1*np.sin(Angle_inlet), H3])
    S1 = Point([diagonal_ratio*R1*np.cos(alpha),diagonal_ratio*R1*np.sin(alpha),H3])
    S2 = Point([diagonal_ratio*R1*np.sin(alpha),diagonal_ratio*R1*np.cos(alpha),H3])
    S3 = P3.copy().rotate(alpha, [0,0,1], [0,0,0])
    S4 = D3.copy().rotate((np.pi/4 - alpha), [0,0,1], [0,0,0])

    # Then, edge points
    #E1 = P3.copy().rotate(alpha/2, [0,0,1], [0,0,0])
    #E2 = S3.copy().rotate( (np.pi/4 - alpha)/2, [0,0,1], [0,0,0])
    #E3 = D3.copy().rotate( (np.pi/4 - alpha)/2, [0,0,1], [0,0,0])
    #E4 = S4.copy().rotate( alpha/2, [0,0,1], [0,0,0])

    # Then, faces
    core_faces = []
    core_faces.append( Face( [P0.position, SH2.position, D1.position, SH1.position]) )
    core_faces.append( Face( [SH1.position, D1.position, S1.position, P1.position]) )
    core_faces.append( Face( [SH2.position, P2.position, S2.position, D1.position]) )
    core_faces.append( Face( [D1.position, S2.position, D2.position, S1.position]))
    core_faces.append( Face( [P1.position, S1.position, S3.position, P3.position], [None, None, Origin([0,0,H3]), None] ) ) #Arc(E1)
    core_faces.append( Face( [S1.position, D2.position, D3.position, S3.position], [None, None, Origin([0,0,H3]), None] ) ) #Arc(E2)
    core_faces.append( Face( [D2.position, S2.position, S4.position, D3.position], [None, None, Origin([0,0,H3]), None] ) ) #Arc(E3)
    core_faces.append( Face( [S2.position, P2.position, P4.position, S4.position], [None, None, Origin([0,0,H3]), None] ) ) #Arc(E4)

    # Finally, make blocks
    core_cone = []
    if H_dip != 0: core_dipleg = []
    if H_bin != 0: core_bin = []
    for j in [0,1,2,3]:
        for i in range(len(core_faces)):
            if j == 0:  topFace = core_faces[i]
            else: topFace = core_faces[i].copy().rotate(j*np.pi/2, [0,0,1], [0,0,0])
            bottomFace = topFace.copy().scale(ratio = coneBottomContraction, origin = [0,0,H3]).translate(displacement = [0,0,-H_cone])
            core_cone.append( Loft( topFace, bottomFace) )
            if H_dip != 0: core_dipleg.append( Extrude(bottomFace, (H2-H1) ) )
            if H_bin != 0:
                core_bin.append( Extrude(bottomFace.copy().translate([0, 0, -H2+H1]), H1) )

    # Chopping along 0X
    core_cone[0].chop(0, count = innerSquareTanCount )
    core_cone[1].chop(1, count = outerSquareTanCount )
    core_cone[4].chop(1, count = rad_count_core)
    core_cone[16].chop(0, count = innerSquareTanCount)
    core_cone[17].chop(1, count = outerSquareTanCount )

    # Chopping along 0Y
    core_cone[0].chop(1, count = innerSquareTanCount)
    core_cone[2].chop(0, count = outerSquareTanCount )
    core_cone[16].chop(1, count = innerSquareTanCount)
    core_cone[18].chop(0, count = outerSquareTanCount)

    # Chopping along 0Z
    core_cone[0].chop(2, count = z_count, c2c_expansion = z_c2c_cone)
    if H_dip != 0: core_dipleg[-1].chop(2, start_size = coreCellLength*(z_c2c_cone)**z_count )

    # Chopping in bin:
    if H_bin != 0:
        # Chopping along 0X
        core_bin[0].chop(0, count = innerSquareTanCount/2 )
        core_bin[1].chop(1, count = outerSquareTanCount/2 )
        core_bin[4].chop(1, count = rad_count_core/2)
        core_bin[16].chop(0, count = innerSquareTanCount/2 )
        core_bin[17].chop(1, count = outerSquareTanCount/2 )

        # Chopping along 0Y
        core_bin[0].chop(1, count = innerSquareTanCount/2 )
        core_bin[2].chop(0, count = outerSquareTanCount/2 )
        core_bin[16].chop(1, count = innerSquareTanCount/2 )
        core_bin[18].chop(0, count = outerSquareTanCount/2 )

        # Chopping along 0Z
        core_bin[0].chop(2, start_size = coreCellLength, c2c_expansion = bin_z_c2c)

    for i in range(len(core_cone)):
        mesh.add(core_cone[i])
        if H_dip != 0:
            core_dipleg[i].set_patch("top","merge")
            mesh.add(core_dipleg[i])
        if H_bin != 0:
            core_bin[i].set_patch("bottom","merge_s")
            binmesh.add(core_bin[i])
########################################

# Extrude top faces upward for core in cylinder region ###
core_cyl = []
bottomFace_outlet = []
for j in [0,1,2,3]:
    for i in range(len(core_faces)):
        if j == 0: bottomFace = core_faces[i].copy()
        else: bottomFace = core_faces[i].copy().rotate((j*np.pi/2), [0,0,1], [0,0,0])
        topFace = bottomFace.copy().translate([0,0,H4-H3])
        bottomFace_outlet.append(topFace)
        core_cyl.append( Loft(topFace, bottomFace))

core_cyl[-1].chop(2, start_size = z_size)

for i in range(len(core_cyl)):
    mesh.add(core_cyl[i])

########################################

# Extra core region in case inlet is larger than VX ###
if H_b > H_vx:
    core_cyl2 = []
    bottomFace_outlet = []
    for j in [0,1,2,3]:
        for i in range(len(core_faces)):
            if j == 0: bottomFace = core_faces[i].copy().translate([0,0,H4-H3])
            else: bottomFace = core_faces[i].copy().rotate(j*np.pi/2, [0,0,1], [0,0,0]).translate([0,0,H4-H3])
            topFace = bottomFace.copy().translate([0,0,H5-H4])
            bottomFace_outlet.append(topFace)
            core_cyl2.append( Loft(topFace, bottomFace) )

    core_cyl2[-1].chop(2, start_size = z_size)

    for i in range(len(core_cyl2)):
        mesh.add(core_cyl2[i])
########################################

# Construct outlet ###
outlet = []
for i in range(len(bottomFace_outlet)):
    topFace_outlet = bottomFace_outlet[i].copy().translate([0,0,H_out])
    outlet.append( Loft(topFace_outlet, bottomFace_outlet[i]))

outlet[-1].chop(2, start_size = z_size)

for i in range(len(outlet)):
    outlet[i].set_patch("bottom","outlet")
    mesh.add(outlet[i])
########################################

# Construct revolution solids (dipleg, cone and cylinder) ###
# Revolution faces for dipleg
if H_dip > 0:
    revolutionFace_vx_dipleg = Face([ [R1b,0,H1], [R1b,0,H2], [R2b,0,H2], [R2b,0,H1]  ])
    revolutionFace_outer_dipleg = Face([ [R2b,0,H1], [R2b,0,H2], [R3b,0,H2], [R3b,0,H1]  ])
    if R_dip > R_b: revolutionFace_rim_dipleg = Face([ [R3b,0,H1], [R3b,0,H2], [R4b,0,H2], [R4b,0,H1]  ])
    vx_dipleg = []
    outer_dipleg = []
    rim_dipleg = []

# Revolution faces for cone
revolutionFace_vx_cone = Face([ [R1b,0,H2], [R1,0,H3], [R2,0,H3], [R2b,0,H2]  ])
revolutionFace_outer_cone = Face([ [R2b,0,H2], [R2,0,H3], [R3,0,H3], [R3b,0,H2]  ])
vx_cone = []
outer_cone = []

# Revolution faces for cylinder
revolutionFace_vx_cyl1 = Face([ [R1,0,H3], [R1,0,H4], [R2,0,H4], [R2,0,H3] ])
revolutionFace_outer_cyl1 = Face([ [R2,0,H3], [R2,0,H4], [R3,0,H4], [R3,0,H3] ])
vx_cyl1 = []
outer_cyl1 = []

revolutionFace_outer_cyl2 = Face([ [R2,0,H4], [R2,0,H5], [R3,0,H5], [R3,0,H4] ])
outer_cyl2 = []

# Revolution solids for bin:
if H_bin != 0:
    revolutionFace_vx_bin = Face([ [R1b,0,H0], [R1b,0,H1], [R2b,0,H1], [R2b,0,H0] ])
    revolutionFace_outer_bin = Face([ [R2b,0,H0], [R2b,0,H1], [R3b,0,H1], [R3b,0,H0] ])
    if R_dip > R_b:
        revolutionFace_rim_bin = Face([ [R3b,0,H0], [R3b,0,H1], [R4b,0,H1], [R4b,0,H0] ])
        revolutionFace_rim2_bin = Face([ [R4b,0,H0], [R4b,0,H1], [R5b,0,H1], [R5b,0,H0] ])
    else:
        revolutionFace_rim2_bin = Face([ [R3b,0,H0], [R3b,0,H1], [R5b,0,H1], [R5b,0,H0] ])
    vx_bin = []
    outer_bin = []
    rim_bin = []
    rim2_bin = []

if H_b != H_vx: # In case inlet is taller or shorter than VX, there is at least one more block
    revolutionFace_outer_cyl3 = Face([ [R2,0,H5], [R2,0,H6], [R3,0,H6], [R3,0,H5] ])
    outer_cyl3 = []

    if H_b > H_vx: # If the inlet is taller than VX, there is another revolution block in VX radius
        revolutionFace_vx_cyl2 = Face([ [R1,0,H4], [R1,0,H5], [R2,0,H5], [R2,0,H4] ])
        vx_cyl2 = []


for i in range(len(Angles) - 1):
    if H_dip > 0:
        # Revolve dipleg
        vx_dipleg.append( Revolve( revolutionFace_vx_dipleg.copy().rotate(Angles[i], [0,0,1], [0,0,0]) , (Angles[i+1] - Angles[i]), [0,0,1], [0,0,0] ) )
        vx_dipleg[i].set_patch("left","merge")
        mesh.add(vx_dipleg[i]) # chop radially

        outer_dipleg.append( Revolve( revolutionFace_outer_dipleg.copy().rotate(Angles[i], [0,0,1], [0,0,0]) , (Angles[i+1] - Angles[i]), [0,0,1], [0,0,0] ) )
        outer_dipleg[i].set_patch("left","merge")
        mesh.add(outer_dipleg[i]) # chop radially

        if R_dip > R_b:
            rim_dipleg.append(Revolve( revolutionFace_rim_dipleg.copy().rotate(Angles[i], [0,0,1], [0,0,0]) , (Angles[i+1] - Angles[i]), [0,0,1], [0,0,0] ))
            rim_dipleg[i].chop(1, count = rad_count_rim)#, c2c_expansion = rad_c2c)
            rim_dipleg[i].set_patch("right","wall")
            rim_dipleg[i].set_patch("back","wall")
            rim_dipleg[i].set_patch("left","merge")
            mesh.add(rim_dipleg[i])

    if H_bin != 0:
        # Revolve bin:
        vx_bin.append( Revolve( revolutionFace_vx_bin.copy().rotate(Angles[i], [0,0,1], [0,0,0]) , (Angles[i+1] - Angles[i]), [0,0,1], [0,0,0] ) )
        if rad_count_vx%2 == 0:
            vx_bin[i].chop(1, count = rad_count_vx/2)#, c2c_expansion = rad_c2c)
        else:
            vx_bin[i].chop(1, count = rad_count_vx)#, c2c_expansion = rad_c2c)
        vx_bin[i].set_patch("right","merge_s")
        binmesh.add(vx_bin[i])

        outer_bin.append( Revolve( revolutionFace_outer_bin.copy().rotate(Angles[i], [0,0,1], [0,0,0]) , (Angles[i+1] - Angles[i]), [0,0,1], [0,0,0] ) )
        outer_bin[i].chop(1, count = rad_count_outer/2)#/2, c2c_expansion = rad_c2c)
        outer_bin[i].set_patch("right","merge_s")
        binmesh.add(outer_bin[i])

        if R_dip > R_b:
            rim_bin.append( Revolve( revolutionFace_rim_bin.copy().rotate(Angles[i], [0,0,1], [0,0,0]) , (Angles[i+1] - Angles[i]), [0,0,1], [0,0,0] ) )
            rim_bin[i].chop(1, count = rad_count_rim/2)#, c2c_expansion = rad_c2c)
            rim_bin[i].set_patch("right","merge_s")
            binmesh.add(rim_bin[i])

        rim2_bin.append( Revolve( revolutionFace_rim2_bin.copy().rotate(Angles[i], [0,0,1], [0,0,0]) , (Angles[i+1] - Angles[i]), [0,0,1], [0,0,0] ) )
        rim2_bin[i].chop(1, start_size = rad_cellN_rim, c2c_expansion = bin_rad_c2c)
        binmesh.add(rim2_bin[i])

    # Revolve cone
    vx_cone.append( Revolve( revolutionFace_vx_cone.copy().rotate(Angles[i], [0,0,1], [0,0,0]) , (Angles[i+1] - Angles[i]), [0,0,1], [0,0,0] ) )
    vx_cone[i].chop(1, count = rad_count_vx)#, c2c_expansion = rad_c2c)
    mesh.add(vx_cone[i])

    outer_cone.append( Revolve( revolutionFace_outer_cone.copy().rotate(Angles[i], [0,0,1], [0,0,0]) , (Angles[i+1] - Angles[i]), [0,0,1], [0,0,0] ) )
    outer_cone[i].chop(1, count = rad_count_outer)#, c2c_expansion = rad_c2c)
    mesh.add(outer_cone[i])

    # Revolve cone to cylinder
    vx_cyl1.append( Revolve( revolutionFace_vx_cyl1.copy().rotate(Angles[i], [0,0,1], [0,0,0]) , (Angles[i+1] - Angles[i]), [0,0,1], [0,0,0] ) )
    if H_b <= H_vx:
        vx_cyl1[i].set_patch("right","wall")
    mesh.add(vx_cyl1[i])


    outer_cyl1.append( Revolve( revolutionFace_outer_cyl1.copy().rotate(Angles[i], [0,0,1], [0,0,0]) , (Angles[i+1] - Angles[i]), [0,0,1], [0,0,0] ) )
    mesh.add(outer_cyl1[i])

    outer_cyl2.append( Revolve( revolutionFace_outer_cyl2.copy().rotate(Angles[i], [0,0,1], [0,0,0]) , (Angles[i+1] - Angles[i]), [0,0,1], [0,0,0] )  )
    if H_b <= H_vx: outer_cyl2[i].chop(0, start_size = z_size)
    mesh.add(outer_cyl2[i])

    # Check for H_vx vs H_b
    if H_b != H_vx:
        outer_cyl3.append( Revolve( revolutionFace_outer_cyl3.copy().rotate(Angles[i], [0,0,1], [0,0,0]) , (Angles[i+1] - Angles[i]), [0,0,1], [0,0,0] ) )
        outer_cyl3[i].chop(0, start_size = z_size)
        mesh.add(outer_cyl3[i])

        if H_b > H_vx:
            vx_cyl2.append( Revolve( revolutionFace_vx_cyl2.copy().rotate(Angles[i], [0,0,1], [0,0,0]) , (Angles[i+1] - Angles[i]), [0,0,1], [0,0,0] ) )
            mesh.add(vx_cyl2[i])
########################################

### Construct inlet
# Inlet structure is affected by two checks: H_b vs H_vx, and Angle_inlet vs pi/4
inletFaces = []

if Angle_inlet == np.pi/4: angles_pos = [1]
elif Angle_inlet > np.pi/4: angles_pos = [3]
else: angles_pos = [1,2,3]

if H_b < H_vx:
    for i in angles_pos:
        inletFaces.append( outer_cyl3[i].get_face("back") )
elif H_b == H_vx:
    for i in angles_pos:
        inletFaces.append( outer_cyl2[i].get_face("back") )
else:
    for i in angles_pos:
        inletFaces.append( outer_cyl2[i].get_face("back") )
        inletFaces.append( outer_cyl3[i].get_face("back") )

extr_inletFaces = []
for face in inletFaces:
    extr_inletFaces.append(Face([
        [L_in, face.points[0].position[1], face.points[0].position[2] ],
        [L_in, face.points[1].position[1], face.points[1].position[2]],
        [L_in, face.points[2].position[1], face.points[2].position[2]],
        [L_in, face.points[3].position[1], face.points[3].position[2] ]]))

inlet = []
for i in range(len(inletFaces)):
    inlet.append( Loft( inletFaces[i], extr_inletFaces[i] ) )

    if len(inletFaces) == 1: # These walls are assigned a special wall patch, so that boundary layers are not added here later on
        inlet[i].set_patch("front","wall")
        inlet[i].set_patch("back","wall")
        inlet[i].set_patch("left","wall")
        inlet[i].set_patch("right","wall")

    elif len(inletFaces) == 2:
        inlet[i].set_patch("front","wall")
        inlet[i].set_patch("back","wall")
        if i == 0:
            inlet[i].set_patch("right","wall")
        else:
            inlet[i].set_patch("left","wall")

    elif len(inletFaces) == 3:
        inlet[i].set_patch("right","wall")
        inlet[i].set_patch("left","wall")
        if i == 0:
            inlet[i].set_patch("back","wall")
        elif i == 2:
            inlet[i].set_patch("front","wall")

    else: # six blocks
        if i in [1, 3, 5]:
            inlet[i].set_patch("right","wall")
        elif i in [0, 2, 4]:
            inlet[i].set_patch("left","wall")
        if i in [0,1]:
            inlet[i].set_patch("back","wall")
        if i in [4, 5]:
            inlet[i].set_patch("front","wall")

    inlet[i].set_patch("top","inlet")
    mesh.add(inlet[i])

inlet[i].chop(2, start_size = inlet_size)
########################################

# Prepare mesh and write blockMeshDict
mesh.set_default_patch("wallLayers", "wall")
mesh.write('bodyMesh/system/blockMeshDict', debug_path="./body.vtk")

if H_bin != 0:
    binmesh.set_default_patch("wall","wall")
    binmesh.write("binMesh/system/blockMeshDict", debug_path="./bin.vtk")
