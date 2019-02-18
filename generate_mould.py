import scipy.io as sio
from stl import mesh
from skimage import measure
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import subprocess
from scipy.spatial import ConvexHull
from solid import *
from solid.utils import *


inputFilePath = sys.argv[1]
outputFilePath = sys.argv[2]
matMatrixName = 'threshold_smooth_box'

debugMode = False

#Everything in mm
gapBetweenSlabs = 1.5
sectionPlaneDistances = 10
sizeInXDirection = 210 # printer limitation
sizeInYDirection = 210 # printer limitation
basePlateThickness = 10
slabHeight = 60
#We cut along x-direction!

tumor= sio.loadmat(inputFilePath)
tree = tumor[matMatrixName]
tree = (tree==1)

print("#######################################\n####### COOKIE CUTTER GENERATOR #######\n#######################################")
print("### Input file: "+inputFilePath)
print("##### Matrix name: "+matMatrixName)
print("### Output file: "+outputFilePath)

print("# Performing marching cubes algorithm on tumor volume")
verts, faces = measure.marching_cubes_classic(tree,)
tumorMesh = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        tumorMesh.vectors[i][j] = verts[f[j],:]
print("# done.")
print("# Save intermediate mesh")
tumorMesh.save("intermediate.stl")
print("# done.")

print("# MeshLab messsages:")
in_file = 'intermediate.stl'
out_file = 'intermediate_proc.stl'
filter_script_path = 'new_smooth.mlx'
 # Add input mesh
command = "meshlabserver -i '"+os.getcwd()+"/" + in_file+"'"
# Add the filter script
command += " -s '"+os.getcwd()+"/" + filter_script_path+"'"
# Add the output filename and output flags
command += " -o '"+os.getcwd()+"/" + out_file + "' -om vn fn"
# Execute command

#print ("Going to execute: " + command)
output = subprocess.check_output(command, shell=True)
#last_line = output.splitlines()[-1]
#print(last_line)
#print("Done:")
#print (in_file + " > " + out_file)
#volume = ConvexHull(tree).volume
print("# done.")

processed_tumour = mesh.Mesh.from_file('intermediate_proc.stl')

testArray = np.array([processed_tumour.vectors[:,0,0],processed_tumour.vectors[:,0,1]]).T
testArray.shape
tumorHull = ConvexHull(testArray,incremental=True)

if debugMode:
    plt.figure(figsize=(10,10))
    plt.ylim([0,210])
    plt.xlim([0,210])
    plt.plot(processed_tumour.vectors[:,0,0],processed_tumour.vectors[:,0,1])
    plt.plot(testArray[tumorHull.vertices,0], testArray[tumorHull.vertices,1], 'r--', lw=2)
    plt.plot(testArray[tumorHull .vertices[0],0], testArray[tumorHull .vertices[0],1], 'ro')
    plt.show()

convexHullExtrude = linear_extrude(slabHeight)(translate([10,10,0])(offset(r=5)(polygon([[testArray[a,0], testArray[a,1]] for a in tumorHull.vertices]))))
convexHullExtrude2 = linear_extrude(slabHeight)(translate([10,10,0])(offset(r=10)(polygon([[testArray[a,0], testArray[a,1]] for a in tumorHull.vertices]))))
tmpSlabs = []
for i in range(sizeInYDirection//sectionPlaneDistances):
    tmpSlabs += translate( [0,i*sectionPlaneDistances,basePlateThickness])(
    cube([sizeInXDirection,sectionPlaneDistances-(gapBetweenSlabs/2),slabHeight])
    )
d = intersection()(cube([sizeInXDirection,sizeInYDirection,basePlateThickness]),convexHullExtrude2)
d += intersection()(tmpSlabs,convexHullExtrude)
d -= translate( [10,10,basePlateThickness])(import_stl("intermediate_proc.stl"))

scad_render_to_file(d,'test.scad')

# Add input mesh
in_file = "test.scad"
out_file = outputFilePath
command = "openscad " + in_file+""
# Add the output filename and output flags
command += " -o " + out_file
# Execute command
# We do:
# Merge close vertices
# Remesh to 10k faces
# Invert face normals
print(command)
print ("Going to execute: " + command)
output = subprocess.check_output(command, shell=True)
print("Done:")
print (in_file + " > " + out_file)
#volume = ConvexHull(tree).volume

# tidy up
os.remove("test.scad")
os.remove("intermediate.stl")
os.remove("intermediate_proc.stl")
