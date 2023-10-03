# classyCyclone
This is a tool for the generation of cyclone separator meshes using the classy_blocks library.

# Prerequisites
- OpenFOAM (any version)
- classy_blocks (https://github.com/damogranlabs/classy_blocks)

# How to use
1. Edit the geometry.py file to input the desired geometry data. Cylindrical diplegs and bins can be included in the mesh if so desired.
2. Run Allrun.sh script

The (unscaled) body mesh can be found in bodyMesh/constant/polyMesh. If a bin is enabled, the body mesh is merged with the bin mesh (grid is coarse here, 2:1 ratio) - the (unscaled) bin mesh is found at binMesh/constant/polyMesh, while the merged mesh is found at mergedMesh/1/polyMesh.
