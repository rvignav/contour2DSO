import pydicom as dicom
import numpy as np
import argparse
from PIL import Image  
import glob
import os
import json
from matplotlib.path import Path
from matplotlib import pyplot as plt

series_path = "Series"

studies = glob.glob('SamplePatient/*')
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

seriesDCMPaths = glob.glob(str(abs_path) + "/*.dcm")

shape = dicom.read_file(seriesDCMPaths[0]).pixel_array.shape[::-1]

# TODO: Order by ImagePatientPosition, not SliceLocation
def bubble_sort(series, series2):
    swapped = True
    while swapped:
        swapped = False
        for i in range(len(series) - 1):
            if series[i].SliceLocation > series[i + 1].SliceLocation:
                series[i], series[i + 1] = series[i + 1], series[i]
                series2[i], series2[i + 1] = series2[i + 1], series2[i]
                swapped = True
    return [series, series2]

seriesDCM = seriesDCMPaths.copy()
for i in range(len(seriesDCMPaths)):
  seriesDCM[i] = dicom.read_file(seriesDCMPaths[i])  
seriesDCM, seriesDCMPaths = bubble_sort(seriesDCM, seriesDCMPaths)

seg_mask = np.zeros((shape[0], shape[1], len(seriesDCM)))

with open(glob.glob('files/*')[0]) as f:
  init_data = json.load(f)
data = init_data["ImageAnnotationCollection"]["imageAnnotations"]["ImageAnnotation"][0]["markupEntityCollection"]["MarkupEntity"]

name = init_data["ImageAnnotationCollection"]["imageAnnotations"]["ImageAnnotation"][0]["name"]["value"]

def find_z(path):
  for i in range(len(seriesDCMPaths)):
    if path in seriesDCMPaths[i]:
      return i

def visualize(seg_mask, filename):
  vis = np.zeros((seg_mask.shape[0], seg_mask.shape[1]))
  for r in range(vis.shape[0]):
    for c in range(vis.shape[1]):
      vis[r][c] = seg_mask[r][c]

  binary = vis > 0
  plt.imshow(binary, cmap='gray')
  plt.gca().set_axis_off()
  plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
  plt.margins(0,0)
  plt.gca().xaxis.set_major_locator(plt.NullLocator())
  plt.gca().yaxis.set_major_locator(plt.NullLocator())
  plt.savefig(filename, bbox_inches = 'tight', pad_inches = 0)

os.mkdir("vis")

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

  visualize(mask, 'vis/img' + str(z) + '.png')

print("Seg_mask array created")