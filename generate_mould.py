import argparse
import os
import subprocess
import sys
from ftplib import FTP

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import stl
from scipy import ndimage
from scipy.spatial import ConvexHull
from shapely.geometry import LineString, Polygon
from skimage import measure
from solid import *
from solid.utils import *
from stl import mesh
import time
import math

# Parse all arguments from command line
parser = argparse.ArgumentParser()
parser.add_argument("inputFile", type=str,
                    help="name of .mat input file containing 3D volume matrix")
parser.add_argument("patientIdentifier", type=str,
                    help="name of the patient to be placed in final file names")
parser.add_argument("--tumourVolume", type=str, default="threshold_smooth_box",
                    help="name of the tumour volume in the .mat input file")
parser.add_argument("--tumourPoints", type=str, default="labels_smooth_box",
                    help="name of the volume points in the .mat input file")
parser.add_argument('--gapBetweenSlabs', type=float, default=1.5)
parser.add_argument('--sectionPlaneDistances', type=float, default=10)
parser.add_argument('--basePlateThickness', type=float, default=4)
parser.add_argument('--printer', type=str, default='prusa-mk3s')
# parser.add_argument("show_convex_hull",
#                help="displays the convex hull of tumor 2D projection")


# Some housekeeping for functions we need later
def ftp_upload(ftp_obj, path, ftype='TXT'):
    """
    A function for uploading files to an FTP server
    @param ftp_obj: The file transfer protocol object
    @param path: The path to the file to upload
    """
    if ftype == 'TXT' or ftype == 'STL':
        with open(path, 'rb') as fobj:
            ftp.storlines('STOR ' + path, fobj)
    else:
        with open(path, 'rb') as fobj:
            ftp.storbinary('STOR ' + path, fobj, 1024)


def convert_scad_to_stl(scadFileName, stlFileName):
    print("# Final conversion to .stl file:")
    # Define in and out files
    in_file = scadFileName + ".scad"
    out_file = stlFileName + ".stl"
    # Add the output filename and output flags
    command = "openscad " + in_file + " -o " + out_file

    # Execute command
    print("Going to execute: " + command)
    output = subprocess.check_output(command, shell=True)
    print("Done.")
    os.remove(in_file)


# Assign all arguments to more readable variables
args = parser.parse_args()

# High-level variables (files etc.)
tumorVolumeNameInMat = args.tumourVolume
tumorPointsNameInMat = args.tumourPoints
inputFilePath = args.inputFile
patientIdentifier = args.patientIdentifier
printerType = args.printer
outputFilePath = "Patient_" + patientIdentifier + "_"

# All units are given in millimeter
gapBetweenSlabs = args.gapBetweenSlabs  # Adapt to blade thickness
sectionPlaneDistances = args.sectionPlaneDistances  # Default in our case
basePlateThickness = args.basePlateThickness  # Thickness of mould base plate

mouldSlabHeight = 47  # this
contactHoleSize = 10  # size of holes in print

sectionDirection = 0  # 0 = Along X?, 1 = Along Y?

tumourVolumeIndex = 1
kidneyPointIndex = 2
tumorPointIndex = 3
vesselPointIndices = [6]  # DODGY, REPLACE!!!
contactPointIndices = [7, 8]  # DODGY, REPLACE!!! Sometimes with 8
pointLabels = ['V', 'CK', 'CT']

# Specify printer limitations (i.e. Prusa MK3S)
if printerType == 'prusa-mk3s':
    sizeInXDirection = 210
    sizeInYDirection = 210
    sizeinZDirection = 250
else:
    raise ValueError('Only Prusa i3 MK3s are supported in this version. Please modify the code to accomodate your printer limitations.')


annotatedSegmentationVolume = sio.loadmat(inputFilePath)
tree = annotatedSegmentationVolume[tumorVolumeNameInMat]
tree = (tree == tumourVolumeIndex)

treePoints = annotatedSegmentationVolume[tumorPointsNameInMat]

# Configuration for cutter
barDepth = 10
columnHeight = mouldSlabHeight + 25


print("#######################################\n####### TISSUE CUTTER GENERATOR #######\n#######################################")
print("### Input file: " + inputFilePath)
print("##### Matrix name: " + tumorVolumeNameInMat)
print("### Output file prefix: " + outputFilePath + "*")

print("# Performing marching cubes algorithm on tumor volume...")
verts, faces = measure.marching_cubes_classic(tree,)
tumorMesh = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        tumorMesh.vectors[i][j] = verts[f[j], :]
print("# done.")
print("# Save intermediate mesh...")
tumorMesh.save(outputFilePath + "tumour_preproc.stl")
print("# done.")

print("# Process and re-save intermediate mesh...")
in_file = outputFilePath + "tumour_preproc.stl"
out_file = outputFilePath + "tumour_postproc.stl"
filter_script_path = 'filter.mlx'
# Add input mesh
command = "meshlabserver -i '" + os.getcwd() + "/" + in_file + "'"
# Add the filter script
command += " -s '" + os.getcwd() + "/" + filter_script_path + "'"
# Add the output filename and output flags
command += " -o '" + os.getcwd() + "/" + out_file + "' -om vn fn"
# Execute command

#print ("Going to execute: " + command)
output = subprocess.check_output(command, shell=True, stderr=open(outputFilePath +'meshlab.log', 'w'))
#last_line = output.splitlines()[-1]
# print(last_line)
# print("Done:")
#print (in_file + " > " + out_file)
#volume = ConvexHull(tree).volume
print("# done.")

print("# Convex hull generation for 2D projection...")
processed_tumour = mesh.Mesh.from_file(outputFilePath + "tumour_postproc.stl")
testArray = np.array([processed_tumour.vectors[:, 0, 0],
                      processed_tumour.vectors[:, 0, 1]]).T
testArray.shape
tumorHull = ConvexHull(testArray, incremental=True)

#if args.show_convex_hull:
    # plt.figure(figsize=(10, 10))
    # plt.ylim([0, 210])
    # plt.xlim([0, 210])
    # plt.plot(processed_tumour.vectors[:, 0, 0],
    #          processed_tumour.vectors[:, 0, 1])
    # plt.plot(testArray[tumorHull.vertices, 0],
    #          testArray[tumorHull.vertices, 1], 'r--', lw=2)
    # plt.plot(testArray[tumorHull .vertices[0], 0],
    #          testArray[tumorHull .vertices[0], 1], 'ro')
    # plt.show()


convexHullArray = np.array([[testArray[a, 0], testArray[a, 1]]
                            for a in tumorHull.vertices])

convexHullPolygon = Polygon(
    [[testArray[a, 0], testArray[a, 1]] for a in tumorHull.vertices])
mbr_points = list(
    zip(*convexHullPolygon.minimum_rotated_rectangle.exterior.coords.xy))
mbr_lengths = [LineString(
    (mbr_points[i], mbr_points[i + 1])).length for i in range(len(mbr_points) - 1)]
minor_axis = min(mbr_lengths)
major_axis = max(mbr_lengths)

print('# Minor axis: ' + str(round(minor_axis, 2)) + 'mm / Major axis: ' + str(round(major_axis, 2)) + 'mm')
print("# done.")

print("# Build mould structure:")
convexHullExtrude = linear_extrude(mouldSlabHeight)(translate([0, 0, 0])(offset(r=5)(
    polygon([[testArray[a, 0], testArray[a, 1]] for a in tumorHull.vertices]))))
convexHullExtrude2 = linear_extrude(mouldSlabHeight)(translate([0, 0, 0])(offset(r=10)(
    polygon([[testArray[a, 0], testArray[a, 1]] for a in tumorHull.vertices]))))
tmpSlabs = []
numberOfSections = sizeInYDirection // sectionPlaneDistances if sectionDirection == 1 else sizeInXDirection // sectionPlaneDistances
for i in range(numberOfSections):
    translateInX = i * sectionPlaneDistances if sectionDirection == 0 else 0
    translateInY = i * sectionPlaneDistances if sectionDirection == 1 else 0
    slabSize = [sizeInXDirection, sectionPlaneDistances - (gapBetweenSlabs / 2), mouldSlabHeight] if sectionDirection == 1 else [
        sectionPlaneDistances - (gapBetweenSlabs / 2), sizeInYDirection, mouldSlabHeight]
    tmpSlabs += translate([translateInX, translateInY, basePlateThickness])(
        cube(slabSize)
    )
d = intersection()(cube(
    [sizeInXDirection, sizeInYDirection, basePlateThickness]), convexHullExtrude2)
d += intersection()(tmpSlabs, convexHullExtrude)
d -= translate([0, 0, basePlateThickness])(import_stl(outputFilePath + "tumour_postproc.stl"))


quit()



# Create vessel and contact point holes
for pointIdx, contactHole in enumerate(vesselPointIndices + contactPointIndices):
    centerOfMassCH = ndimage.measurements.center_of_mass(
        treePoints == contactHole)
    d = difference()(d, translate([centerOfMassCH[0], centerOfMassCH[1], 0])(
        cylinder(h=sizeinZDirection, r=contactHoleSize)))
    d = difference()(d, translate([centerOfMassCH[0] + 1.2 * contactHoleSize, centerOfMassCH[1] + 1.2 * contactHoleSize, 0])(
        mirror([1, 0, 0])(linear_extrude(1)(text(pointLabels[pointIdx], font="Liberation Sans", size=5)))))


# Create tumor color part
centerOfMassTumor = ndimage.measurements.center_of_mass(
    treePoints == tumorPointIndex)
centerOfMassKidney = ndimage.measurements.center_of_mass(
    treePoints == kidneyPointIndex)

d_kidney = intersection()(translate([centerOfMassKidney[0], centerOfMassKidney[1],
                                     basePlateThickness - 1])(cylinder(h=sizeinZDirection, r=contactHoleSize * 1.5)), d)
d = difference()(d, translate([centerOfMassKidney[0], centerOfMassKidney[1],
                               basePlateThickness - 1])(cylinder(h=sizeinZDirection, r=contactHoleSize * 1.5)))
d_tumour = intersection()(translate([centerOfMassTumor[0], centerOfMassTumor[1],
                                     basePlateThickness - 1])(cylinder(h=sizeinZDirection, r=contactHoleSize * 1.5)), d)
d = difference()(d, translate([centerOfMassTumor[0], centerOfMassTumor[1],
                               basePlateThickness - 1])(cylinder(h=sizeinZDirection, r=contactHoleSize * 1.5)))


# Build a baseplate
# Determind minimum and maximum extend in alldirections
print(convexHullArray.shape)
minX = np.min(convexHullArray[:, 0])
maxX = np.max(convexHullArray[:, 0])
minY = np.min(convexHullArray[:, 1])
maxY = np.max(convexHullArray[:, 1])
print(minX, maxX, minY, maxY)
# print(sizeInXDirection/2-(maxX-minX)/2)
#d = translate([sizeInYDirection/2-(maxY-minY)/2,sizeInXDirection/2-(maxX-minX)/2,0])(d)
#d += translate([0,0,0])(cube([200,30,basePlateThickness]))

d = translate([barDepth, barDepth, 0])(d)
d_tumour = translate([barDepth, barDepth, 0])(d_tumour)
d_kidney = translate([barDepth, barDepth, 0])(d_kidney)

scad_render_to_file(d, 'test.scad')
scad_render_to_file(d_tumour, 'test_tumour.scad')
scad_render_to_file(d_kidney, 'test_kidney.scad')
print("# done.")


# find the max dimensions, so we can know the bounding box, getting the height,
# width, length (because these are the step size)...


def find_mins_maxs(obj):
    minx = maxx = miny = maxy = minz = maxz = None
    for p in obj.points:
        # p contains (x, y, z)
        if minx is None:
            minx = p[stl.Dimension.X]
            maxx = p[stl.Dimension.X]
            miny = p[stl.Dimension.Y]
            maxy = p[stl.Dimension.Y]
            minz = p[stl.Dimension.Z]
            maxz = p[stl.Dimension.Z]
        else:
            maxx = max(p[stl.Dimension.X], maxx)
            minx = min(p[stl.Dimension.X], minx)
            maxy = max(p[stl.Dimension.Y], maxy)
            miny = min(p[stl.Dimension.Y], miny)
            maxz = max(p[stl.Dimension.Z], maxz)
            minz = min(p[stl.Dimension.Z], minz)
    return minx, maxx, miny, maxy, minz, maxz


main_body = mesh.Mesh.from_file('Patient.stl')

minx, maxx, miny, maxy, minz, maxz = find_mins_maxs(main_body)

# the logic is easy from there
print("X:", maxx - minx)
print("Y:", maxy - miny)
print("Z:", maxz - minz)

# Build cut guide ushape
d_guide = translate([0, 0, 0])(
    cube([(maxx - minx) + 2 * barDepth, barDepth, basePlateThickness]))
d_guide += translate([0, barDepth, 0]
                     )(cube([barDepth, maxy - miny,  basePlateThickness]))
d_guide += translate([maxx - minx + barDepth, barDepth, 0]
                     )(cube([barDepth, maxy - miny,  basePlateThickness]))

numberOfSectionsGuide = math.ceil((
    maxy - miny) // sectionPlaneDistances if sectionDirection == 1 else (maxx - minx) // sectionPlaneDistances)
for i in range(numberOfSectionsGuide):
    translateInX = i * sectionPlaneDistances if sectionDirection == 0 else 0
    translateInY = i * sectionPlaneDistances if sectionDirection == 1 else 0
    slabSize = [barDepth, sectionPlaneDistances - (gapBetweenSlabs / 2), columnHeight] if sectionDirection == 1 else [
        sectionPlaneDistances - (gapBetweenSlabs / 2), barDepth, columnHeight]
    d_guide += translate([barDepth + translateInX, translateInY, basePlateThickness])(
        cube(slabSize))
    d_guide -= translate([barDepth + translateInX-4, translateInY+4, basePlateThickness + columnHeight-1])(
        linear_extrude(3)(rotate([0,0,180])(text(str(i), font="Liberation Sans", size=3.5))))
    d_guide -= translate([barDepth + translateInX, translateInY+6, basePlateThickness + columnHeight-1])(
        cube([2,2,2]))
    #d_guide += translate([maxx - minx + 1, barDepth + translateInY, basePlateThickness])(
    #    cube(slabSize))

scad_render_to_file(d_guide, 'test_guide.scad')

scad_render_to_file(d+d_guide, 'test_combined_guide.scad')


filesToConvert = ['test.scad', 'test_tumour.scad',
                  'test_kidney.scad', 'test_guide.scad','test_combined_guide.scad']
targetFileNames = [outputFilePath, outputFilePath +
                   "_tumour", outputFilePath + "_kidney", outputFilePath + "_guide", outputFilePath + "_combined_guide"]
for fileIdx, fileToConvert in enumerate(filesToConvert):
    print("# Final conversion to .stl file:")
    # Add input mesh
    in_file = fileToConvert
    out_file = targetFileNames[fileIdx] + ".stl"
    command = "openscad " + in_file + ""
    # Add the output filename and output flags
    command += " -o " + out_file
    # Execute command

    print(command)
    print("Going to execute: " + command)
    output = subprocess.check_output(command, shell=True)
    print("Done:")
    print(in_file + " > " + out_file)
    os.remove(fileToConvert)
#volume = ConvexHull(tree).volume
print("# Projected tumor size (max): \n## Major axis:" +
      str(major_axis) + " \n## Minor axis:" + str(minor_axis))
print("# done. Happy 3D printing!")

# print("# Upload to kidney.cyted.io")
ftp = FTP('cyted.io')
print ("Welcome: ", ftp.getwelcome())
ftp.login("kidney@cyted.io", "#")
ftp_upload(ftp, "Patient.stl", ftype="STL")
ftp_upload(ftp, "Patient_tumour.stl", ftype="STL")
ftp_upload(ftp, "Patient_kidney.stl", ftype="STL")
ftp_upload(ftp, "Patient_guide.stl", ftype="STL")
ftp.quit()

#output = subprocess.check_output("/home/gehrun01/Applications/Slic3rPE-1.41.3/Slic3rPE-1.41.3+linux64-full-201902121303.AppImage "+outputFilePath+".stl "+outputFilePath+"_tumour.stl "+outputFilePath+"_kidney.stl", shell=True)
