# 3D printed tumor moulds for guided tissue sampling
## Requirements
* Python (with modules as listed in requirements.txt)
* Meshlab (with meshlabserver)
* OpenSCAD
## Example
`python generate_mould.py example_data/aligned_tumour.mat`
will result in three .stl files:
1. Patient.stl (containing the mould with holes for projected centroids of tumour and kidney)
2. Patient_tumour.stl (containing the fill-in for the tumour centroid hole)
3. Patient_kidney.stl (containing the fill-in for the kidney centroid hole)
