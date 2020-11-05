import pydicom as dicom
import numpy as np
import argparse
from PIL import Image  
import glob
import os
import json
from matplotlib.path import Path
from matplotlib import pyplot as plt
import timeit

parser = argparse.ArgumentParser(description='Create DSO')
parser.add_argument('series_path', type=str, help='Path to series')

args = parser.parse_args()

series_path = args.series_path

start = timeit.default_timer()

studies = glob.glob('/home/series/PatientSeries/*')
paths = []
for study in studies:
    series_paths = glob.glob(study + '/*')
    for item in series_paths:
        paths.append(item)

abs_path = ''
for path in paths:
    if series_path in path:
        abs_path = path
        break

def bubble_sort(series, series2):
    swapped = True
    while swapped:
        swapped = False
        for i in range(len(series) - 1):
            if series[i].ImagePositionPatient[2] > series[i + 1].ImagePositionPatient[2]:
                series[i], series[i + 1] = series[i + 1], series[i]
                series2[i], series2[i + 1] = series2[i + 1], series2[i]
                swapped = True
    return [series, series2]

seriesDCMPaths = glob.glob(str(abs_path) + "/*.dcm")

shape = dicom.read_file(seriesDCMPaths[0]).pixel_array.shape[::-1]

seriesDCM = seriesDCMPaths.copy()
for i in range(len(seriesDCMPaths)):
  seriesDCM[i] = dicom.read_file(seriesDCMPaths[i])  
seriesDCM, seriesDCMPaths = bubble_sort(seriesDCM, seriesDCMPaths)

seriesDCM.reverse()
seriesDCMPaths.reverse()

if not os.path.isdir('output'):
  os.mkdir('output')

seg_mask = np.zeros((shape[0], shape[1], len(seriesDCM)))

with open(glob.glob('/home/series/files/*')[0]) as f:
  init_data = json.load(f)
data = init_data["ImageAnnotationCollection"]["imageAnnotations"]["ImageAnnotation"][0]["markupEntityCollection"]["MarkupEntity"]

name = init_data["ImageAnnotationCollection"]["imageAnnotations"]["ImageAnnotation"][0]["name"]["value"]

def find_z(path):
  for i in range(len(seriesDCMPaths)):
    if path in seriesDCMPaths[i]:
      return i

for img in data:
  coords = img["twoDimensionSpatialCoordinateCollection"]["TwoDimensionSpatialCoordinate"]
  path = img["imageReferenceUid"]["root"]
  z = find_z(path)
  tupVerts = []
  for coord in coords:
    x_val = round(coord["x"]["value"])
    y_val = round(coord["y"]["value"])
    tupVerts.append((x_val, y_val))
  
  x, y = np.meshgrid(np.arange(seg_mask.shape[0]), np.arange(seg_mask.shape[1])) # make a canvas with coordinates
  x, y = x.flatten(), y.flatten()
  points = np.vstack((x,y)).T 

  p = Path(tupVerts) # make a polygon
  grid = p.contains_points(points)
  mask = grid.reshape(seg_mask.shape[0],seg_mask.shape[1])

  for r in range(seg_mask.shape[0]):
    for c in range(seg_mask.shape[1]):
      seg_mask[r][c][z] = mask[r][c]

print("Seg_mask array created")
stop = timeit.default_timer()
print("Runtime (mask creation): " + str(stop-start) + "s")

start2 = timeit.default_timer()

import write_dso
import numpy as np
import matlab
maskarray = matlab.double(seg_mask.tolist())

stop2 = timeit.default_timer()
print("Runtime (MATLAB imports and matlab.double): " + str(stop2-start2) + "s")

start3 = timeit.default_timer()

lib = write_dso.initialize()
lib.write_DSO(abs_path, name, maskarray, "./output/")

stop3 = timeit.default_timer()
print("DSO object saved to the folder 'output'")
print("Runtime (write_DSO): " + str(stop3-start3) + "s")