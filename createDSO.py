import pydicom as dicom
import numpy as np
import argparse
from PIL import Image  
import glob
import os
import progressbar

parser = argparse.ArgumentParser(description='Create DSO')
parser.add_argument('series_path', type=str, help='Path to series')

args = parser.parse_args()

series_path = args.series_path

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

# TODO: Order by ImagePatientPosition, not SliceLocation
def bubble_sort(series):
    swapped = True
    while swapped:
        swapped = False
        for i in range(len(series) - 1):
            if series[i].SliceLocation > series[i + 1].SliceLocation:
                series[i], series[i + 1] = series[i + 1], series[i]
                swapped = True
    return series

seriesDCM = glob.glob(str(abs_path) + "/*.dcm")

bar = progressbar.ProgressBar(maxval=len(series1DCM)/5 + 1, \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])

shape = dicom.read_file(seriesDCM[0]).pixel_array.shape[::-1]

for i in range(len(seriesDCM)):
  seriesDCM[i] = dicom.read_file(seriesDCM[i])  
seriesDCM = bubble_sort(seriesDCM)

if not os.path.isdir('output'):
  os.mkdir('output')

seg_mask = np.zeros((shape[0], shape[1], len(seriesDCM))

bar.start()
for i in range(len(seriesDCM)):
  # TODO: Create seg_mask array
  if i % 5 == 0:
    bar.update(i/5 + 1)
bar.finish()

import matlab.engine
eng = matlab.engine.start_matlab()
eng.write_DSO([series_path, 'DSO', seg_mask, 'output', False, False])
eng.quit()

print("DSO object saved to the folder 'output'")