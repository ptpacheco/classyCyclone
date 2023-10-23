#!/bin/sh

# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

python3.8 classyCyclone.py
cd bodyMesh/
rm -r constant/polyMesh
rm -r 1
blockMesh

stitchMesh -perfect master_topStitch1 slave_topStitch1
rm -r constant/polyMesh
cp -r 1/polyMesh constant/
rm -r 1/

stitchMesh -perfect master_topStitch2 slave_topStitch2
rm -r constant/polyMesh
cp -r 1/polyMesh constant/
rm -r 1/

stitchMesh master_inletStitch slave_inletStitch
rm -r constant/polyMesh
cp -r 1/polyMesh constant/
rm -r 1/

topoSet
createPatch
rm -r constant/polyMesh
cp -r 1/polyMesh constant/
rm -r 1/

renumberMesh
rm -r constant/polyMesh
cp -r 1/polyMesh constant/
rm -r 1/
transformPoints -scale '(0.001 0.001 0.001)'
