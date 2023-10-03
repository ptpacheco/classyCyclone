#!/bin/sh
# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# EDIT THIS INCASE THERE IS NO BIN
has_bin=true

# Make main mesh
python3.8 make_mesh.py
cd ./bodyMesh
./Allclean.sh
runApplication blockMesh

# Make bin, if it exists
if [ "$has_bin" = true ]; then
  cd ../binMesh
  ./Allclean.sh
  runApplication blockMesh
  # Merge two meshes
  cd ../mergedMesh
  ./Allclean.sh
  cp -r ../bodyMesh/constant/polyMesh constant/
  runApplication mergeMeshes -overwrite . ../binMesh/
  runApplication stitchMesh -partial -toleranceDict toleranceDict merge_s merge
else
  transformPoints -scale '(0.001 0.001 0.001)'
fi
