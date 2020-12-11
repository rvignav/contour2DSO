"""createDSO
createDSO, and the accompanying methods bubble_sort and find_z, are used to 
automatically create a DICOM Segmentation Object (DSO) from an AIM file 
storing the contours of a DICOM image series. Integrated into Stanford 
University's ePAD Imaging Platform (https://epad.stanford.edu/) as a plugin.

MIT License

Copyright (c) 2020 Vignav Ramesh, ePAD Team (Stanford University)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

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
from pydicom.uid import ExplicitVRLittleEndian, generate_uid
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from datetime import datetime

parser = argparse.ArgumentParser(description='Create DSO')

args = parser.parse_args()
start = timeit.default_timer()

# Store all DICOM image file paths in 'paths'.
series_paths = glob.glob('/home/series/PatientSeries/*')
paths = []
for item in series_paths:
    paths.append(item)
    
aimfiles = glob.glob('/home/series/files/*')
aimfile = aimfiles[0]

series_path = aimfile[aimfile.rindex('/home/series/files/') + len('/home/series/files/'):aimfile.rindex('.json')]

abs_path = ''
for path in paths:
    if series_path in path:
        abs_path = path
        break

def bubble_sort(series, series2):
    """bubble_sort

    :param series: list of DICOM objects in input series
    :param series2: list of filenames corresponding to each
                    element in series
    
    Sort 'series' and 'series2' by the z-index of the
    ImagePositionPatient attribute of the elements in 
    'series'.
    """
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

# Store dimensions of individual DICOM image in 'shape'.
shape = dicom.read_file(seriesDCMPaths[0]).pixel_array.shape[::-1]

seriesDCM = seriesDCMPaths.copy()
for i in range(len(seriesDCMPaths)):
  seriesDCM[i] = dicom.read_file(seriesDCMPaths[i])  
seriesDCM, seriesDCMPaths = bubble_sort(seriesDCM, seriesDCMPaths)

seriesDCM.reverse()
seriesDCMPaths.reverse()

# Create 'output' directory inside DICOM container.
if not os.path.isdir('/output'):
  os.mkdir('/output')

seg_mask = np.zeros((shape[0], shape[1], len(seriesDCM)))

# Load data from AIM file in '/home/series/files' directory in 'data'.
with open(aimfile) as f:
  init_data = json.load(f)
data = init_data["ImageAnnotationCollection"]["imageAnnotations"]["ImageAnnotation"][0]["markupEntityCollection"]["MarkupEntity"]

# Store name of annotation in 'name'.
name = init_data["ImageAnnotationCollection"]["imageAnnotations"]["ImageAnnotation"][0]["name"]["value"]

def find_z(path):
  """find_z

  :param path: path of DICOM file
  
  Return the z-index of the DICOM file specified by 'path'.
  """
  for i in range(len(seriesDCMPaths)):
    if path in seriesDCMPaths[i]:
      return i

for img in data:
  """Store outline of annotation as a list of points with values specified in 
  ["twoDimensionSpatialCoordinateCollection"]["TwoDimensionSpatialCoordinate"]
  attribute. Generate a polygon with all points contained by the outline and 
  store in 'seg_mask' at the appropriate indices.
  """

  coords = img["twoDimensionSpatialCoordinateCollection"]["TwoDimensionSpatialCoordinate"]
  path = img["imageReferenceUid"]["root"]
  z = find_z(path)
  tupVerts = []
  for coord in coords:
    x_val = round(coord["x"]["value"])
    y_val = round(coord["y"]["value"])
    tupVerts.append((x_val, y_val))
  
  x, y = np.meshgrid(np.arange(seg_mask.shape[0]), np.arange(seg_mask.shape[1]))
  x, y = x.flatten(), y.flatten()
  points = np.vstack((x,y)).T 

  p = Path(tupVerts)
  grid = p.contains_points(points)
  mask = grid.reshape(seg_mask.shape[0],seg_mask.shape[1])

  for r in range(seg_mask.shape[0]):
    for c in range(seg_mask.shape[1]):
      seg_mask[r][c][z] = mask[r][c]

print("Seg_mask array created")

stop = timeit.default_timer()
print("Runtime (mask creation): " + str(stop-start) + "s")

start = timeit.default_timer()

startIndex = 0
endIndex = 0

# Optimize seg_mask by removing all empty elements at the beginning and end of the array.
for z in range(0, len(seriesDCM)):
  no1s = True
  for r in range(seg_mask.shape[0]):
    for c in range(seg_mask.shape[1]):
      if (seg_mask[r][c][z]):
        no1s = False
        break 
    else:
      continue
    break
  if no1s:
    continue
  else:
    startIndex = z
    break

for z in range(0, len(seriesDCM)):
  no1s = True
  for r in range(seg_mask.shape[0]):
    for c in range(seg_mask.shape[1]):
      if (seg_mask[r][c][len(seriesDCM) - z - 1]):
        no1s = False
        break 
    else:
      continue
    break
  if no1s:
    continue
  else:
    endIndex = len(seriesDCM) - z - 1
    break

# Order needs to be z x y.
smask = np.zeros((endIndex-startIndex+1, seg_mask.shape[0], seg_mask.shape[1]),dtype=np.uint8)
for z in range(startIndex, endIndex + 1):
  for r in range(seg_mask.shape[0]):
    for c in range(seg_mask.shape[1]):
      smask[z - startIndex][r][c] = seg_mask[r][c][z]

# Create FileMetaDataset storing meta information of the DSO 'info_mask'.
file_meta = FileMetaDataset()
file_meta.FileMetaInformationVersion = b'\x00\x01'
file_meta.MediaStorageSOPClassUID='1.2.840.10008.5.1.4.1.1.66.4'
dicomuid = generate_uid()
instanceuid=dicomuid
file_meta.MediaStorageSOPInstanceUID=instanceuid
file_meta.ImplementationClassUID='1.2.840.10008.5.1.4.1.1.66.4'
file_meta.ImplementationVersionName='ePAD_python_1.0'

#Set is_little_endian and is_implicit_VR attributes of info_mask.
suffix = '.dcm'
info = seriesDCM[0]
info_mask = FileDataset('/output/' + str(name) + str(suffix), {},
                  file_meta=file_meta, preamble=b"\0" * 128)
if 'StudyDescription' in info:
  info_mask.StudyDescription=info.StudyDescription
info_mask.file_meta.TransferSyntaxUID=ExplicitVRLittleEndian
info_mask.is_little_endian = True
info_mask.is_implicit_VR = False

if 'ImageOrientationPatient' in info:
  # Pixel Measures Sequence	(0028,9110)
  pms = Dataset()
  pms.SliceThickness = info.SliceThickness
  pms.PixelSpacing = info.PixelSpacing

  # Plane Orientation Sequence (0020,9116)
  pos = Dataset()
  pos.ImageOrientationPatient = info.ImageOrientationPatient
  
  # Segment Identification Sequence (0062,000A)
  sis = Dataset()
  sis.ReferencedSegmentNumber = 1
  
  # Shared Functional Groups Sequence (5200,9229)
  sfgs = Dataset()
  sfgs.SegmentIdentificationSequence = [sis]
  sfgs.PlaneOrientationSequence = [pos]
  sfgs.PixelMeasuresSequence = [pms]
  info_mask.SharedFunctionalGroupsSequence = [sfgs]
else:
    if 'ImageOrientationPatient' in info.PerFrameFunctionalGroupsSequence[0].PlaneOrientationSequence[0]:
      # Plane Orientation Sequence (0020,9116)
      pos = Dataset()
      pos.ImageOrientationPatient = info.PerFrameFunctionalGroupsSequence[0].PlaneOrientationSequence[0].ImageOrientationPatient
      
      # Pixel Measures Sequence	(0028,9110)
      pms = Dataset()
      pms.SliceThickness=info.PerFrameFunctionalGroupsSequence[0].PixelMeasuresSequence[0].SliceThickness
      pms.PixelSpacing=info.PerFrameFunctionalGroupsSequence[0].PixelMeasuresSequence[0].PixelSpacing
            
      # Segment Identification Sequence (0062,000A)
      sis = Dataset()
      sis.ReferencedSegmentNumber = 1

      # Shared Functional Groups Sequence (5200,9229)
      sfgs = Dataset()
      sfgs.SegmentIdentificationSequence = [sis]
      sfgs.PlaneOrientationSequence = [pos]
      sfgs.PixelMeasuresSequence = [pms]
      info_mask.SharedFunctionalGroupsSequence = [sfgs]

if len(seriesDCM) > 1:
    # Referenced Instance Sequence (0008,114A)
    ris = []
    # Per-frame Functional Groups Sequence (5200,9230)
    pffgs = []
    
    for i in range(startIndex, endIndex+1):
        slice_info=seriesDCM[i]
        ib1=i-startIndex+1

        ris_item = Dataset()
        ris_item.ReferencedSOPClassUID=slice_info.SOPClassUID
        ris_item.ReferencedSOPInstanceUID=slice_info.SOPInstanceUID
        ris.append(ris_item)

        # Purpose Of Reference Code Sequence (0040,A170)
        porcs = Dataset()
        porcs.CodeValue='121322'
        porcs.CodingSchemeDesignator='DCM'
        porcs.CodeMeaning='Source image for image processing operation'

        # Source Image Sequence(0008,2112)
        sis = Dataset()
        sis.ReferencedSOPClassUID=slice_info.SOPClassUID
        sis.ReferencedSOPInstanceUID=slice_info.SOPInstanceUID
        sis.PurposeOfReferenceCodeSequence = [porcs]
        
        # Derivation Code Sequence (0008,9215)
        dcs = Dataset()
        dcs.CodeValue='113076'
        dcs.CodingSchemeDesignator='DCM'
        dcs.CodeMeaning='Segmentation'

        # Derivation Image Sequence (0008,9124)
        dis = Dataset()
        dis.SourceImageSequence = [sis]
        dis.DerivationCodeSequence = [dcs]
        
        # Frame Content Sequence (0020,9111)
        fcs = Dataset()
        fcs.StackID='1'
        fcs.InStackPositionNumber=ib1
        fcs.DimensionIndexValues= [1,i-startIndex+1,1] # frames are one indexed

        pffgs_item = Dataset()
        pffgs_item.DerivationImageSequence = [dis]
        pffgs_item.FrameContentSequence = [fcs]
        
        # Plane Position Sequence (0020,9113)
        if 'ImagePositionPatient' in info:
          pps = Dataset()
          pps.ImagePositionPatient=slice_info.ImagePositionPatient
          pffgs_item.PlanePositionSequence = [pps]
        else:
            if 'ImagePositionPatient' in info.PerFrameFunctionalGroupsSequence[0].PlanePositionSequence[0]:
                print('shouldnt come here. why the information is in frames and it is not a multiframe')
                pps = Dataset()
                pps.ImagePositionPatient=slice_info.PerFrameFunctionalGroupsSequence[0].PlanePositionSequence[0].ImagePositionPatient
                pffgs_item.PlanePositionSequence = [pps]
            else:
                print('shouldnt happen. why the information is not there and it is not a multiframe')
        pffgs.append(pffgs_item)

    # Referenced Series Sequence (0008,1115)
    rss = Dataset()
    rss.ReferencedInstanceSequence = ris
    rss.SeriesInstanceUID=info.SeriesInstanceUID
    info_mask.ReferencedSeriesSequence = [rss]
    info_mask.PerFrameFunctionalGroupsSequence = pffgs
else:
    # num of frames is the number of nonempty slices
    numOfFrames=endIndex-startIndex+1
    if numOfFrames<info.NumberOfFrames:
        print('The input mask is smaller than the frames! Assuming they start from the beginning')

    # Referenced Instance Sequence (0008,114A)
    ris_data = Dataset()
    ris_data.ReferencedSOPClassUID=info.SOPClassUID
    ris_data.ReferencedSOPInstanceUID=info.SOPInstanceUID
    ris = [ris_data]
    
    # Referenced Series Sequence (0008,1115)
    rss = Dataset()
    rss.ReferencedInstanceSequence = ris
    rss.SeriesInstanceUID=info.SeriesInstanceUID
    info_mask.ReferencedSeriesSequence = [rss]
    
    # Per-frame Functional Groups Sequence (5200,9230)
    pffgs = []
    for i in range(1, numOfFrames+1):
        # Purpose Of Reference Code Sequence (0040,A170)
        porcs = Dataset()
        porcs.CodeValue='121322'
        porcs.CodingSchemeDesignator='DCM'
        porcs.CodeMeaning='Source image for image processing operation'
        
        # Segment Identification Sequence (0062,000A)
        sis = Dataset()
        sis.ReferencedSOPClassUID=info.SOPClassUID
        sis.ReferencedSOPInstanceUID=info.SOPInstanceUID
        sis.PurposeOfReferenceCodeSequence = [porcs]
        
        # Derivation Code Sequence (0008,9215)
        dcs = Dataset()
        dcs.CodeValue='113076'
        dcs.CodingSchemeDesignator='DCM'
        dcs.CodeMeaning='Segmentation'

        # Derivation Image Sequence (0008,9124)
        dis = Dataset()
        dis.SourceImageSequence = [[sis]]
        dis.DerivationCodeSequence = [dcs]

        # Frame Content Sequence (0020,9111)
        fcs = Dataset()
        fcs.StackID='1'
        fcs.InStackPositionNumber=i
        fcs.DimensionIndexValues= [1,i,1]

        pffgs_item = Dataset()
        pffgs_item.DerivationImageSequence = [dis]
        pffgs_item.FrameContentSequence = [fcs]       
        
        if 'ImagePositionPatient' in info.PerFrameFunctionalGroupsSequence[0].PlanePositionSequence[0]:
          # Plane Position Sequence (0020,9113)
          pps = Dataset()
          pps.ImagePositionPatient=info.PerFrameFunctionalGroupsSequence[0].PlanePositionSequence[0].ImagePositionPatient
          pffgs_item.PlanePositionSequence = [pps]
        else:
            print('shouldnt happen. why the information is not there ')
        
        pffgs.append(pffgs_item)
    info_mask.PerFrameFunctionalGroupsSequence = pffgs

info_mask.ReferringPhysicianName=''
info_mask.PatientName=info.PatientName
info_mask.PatientID=info.PatientID
info_mask.PatientBirthDate= info.PatientBirthDate
info_mask.PatientSex= info.PatientSex
if 'PatientAge' in info:
    info_mask.PatientAge= info.PatientAge

if 'PatientWeight' in info:
    info_mask.PatientWeight= info.PatientWeight

info_mask.StudyID=info.StudyID

info_mask.ImageType='DERIVED\PRIMARY'
info_mask.SOPClassUID='1.2.840.10008.5.1.4.1.1.66.4'
info_mask.SOPInstanceUID= instanceuid

info_mask.AccessionNumber=info.AccessionNumber
info_mask.Modality='SEG'
info_mask.Manufacturer='Stanford University'

info_mask.ManufacturerModelName= 'ePAD Matlab'
info_mask.DeviceSerialNumber='SN123456'
info_mask.SoftwareVersions='1.0'
info_mask.StudyInstanceUID=info.StudyInstanceUID
info_mask.SeriesInstanceUID= dicomuid
info_mask.SeriesNumber= 1000
info_mask.ContentDate=datetime.today().strftime('%Y%m%d')
info_mask.StudyDate=info.StudyDate
info_mask.SeriesDate=datetime.today().strftime('%Y%m%d')
info_mask.AcquisitionDate=datetime.today().strftime('%Y%m%d')
currentTime=datetime.today().strftime('%H%M%S.')
f = datetime.today().strftime('%f')
currentTime += str(f[0:3])
info_mask.ContentTime=currentTime
info_mask.StudyTime=info.StudyTime
info_mask.SeriesTime=currentTime
info_mask.AcquisitionTime=currentTime
info_mask.InstanceNumber= 1
info_mask.FrameOfReferenceUID= info.FrameOfReferenceUID
info_mask.PositionReferenceIndicator= ''

# Dimension Index Sequence (0020,9222)
dis1 = Dataset()
dis1.DimensionIndexPointer = 0x00209056 # converted numbers to hex and put together
dis1.FunctionalGroupPointer = 0x00209111
dis1.DimensionDescriptionLabel='Stack ID'

dis2 = Dataset()
dis2.DimensionIndexPointer = 0x00209057
dis2.FunctionalGroupPointer = 0x00209111
dis2.DimensionDescriptionLabel='In-Stack Position Number'

dis3 = Dataset()
dis3.DimensionIndexPointer = 0x0062000B
dis3.FunctionalGroupPointer = 0x0062000A
dis3.DimensionDescriptionLabel='Referenced Segment Number'

dis4 = Dataset()
dis4.DimensionOrganizationUID= dicomuid

info_mask.DimensionIndexSequence = [dis1,dis2,dis3]
info_mask.DimensionOrganizationSequence = [dis4]

info_mask.SamplesPerPixel= 1
info_mask.PhotometricInterpretation= 'MONOCHROME2'
info_mask.NumberOfFrames=endIndex-startIndex+1 
info_mask.Rows= seg_mask.shape[0]
info_mask.Columns= seg_mask.shape[1]

# for binary masks it is 1 bit
info_mask.BitsAllocated= 1
info_mask.BitsStored= 1
info_mask.HighBit= 0

info_mask.PixelRepresentation= 0
info_mask.LossyImageCompression='00'
info_mask.SegmentationType='BINARY'

# Anatomic Region Sequence (0008,2218)
ars = Dataset()
ars.CodeValue = 'T-D0050'
ars.CodingSchemeDesignator='SRT'
ars.CodeMeaning='Tissue'

# Segmented Property Category Code Sequence (0062,0003)
spccs = Dataset()
spccs.CodeValue = 'T-D0050'
spccs.CodingSchemeDesignator='SRT'
spccs.CodeMeaning='Tissue'

# Segmented Property Type Code Sequence (0062,000F)
sptcs = Dataset()
sptcs.CodeValue = 'T-D0050'
sptcs.CodingSchemeDesignator='SRT'
sptcs.CodeMeaning='Tissue'

# Segment Sequence (0062,0002)
ss = Dataset()
ss.AnatomicRegionSequence = [ars]
ss.SegmentedPropertyCategoryCodeSequence = [spccs]
ss.SegmentNumber = 1
ss.SegmentLabel='Segmentation'
ss.SegmentAlgorithmType='SEMIAUTOMATIC'
ss.SegmentAlgorithmName='ePAD'
ss.SegmentedPropertyTypeCodeSequence = [sptcs]
info_mask.SegmentSequence = [ss]

info_mask.ContentCreatorName='ePAD^python'
info_mask.ContentLabel= 'ROI'

info_mask.ContentDescription=str(name) + 'segmentation'
info_mask.SeriesDescription=str(name) + ' segmentation'

file_name='/output/' + str(name) + '.dcm'
if (os.path.exists(file_name)):
  file_name += currentTime

# binary data needs to be packed for segmentations
info_mask.PixelData = np.packbits(smask,bitorder='little')

# check the dataset for DICOM standard by adding write_like_original=False
info_mask.save_as(file_name, write_like_original=False)
print("info_mask saved to " + file_name)

stop = timeit.default_timer()
print("Runtime (DSO creation): " + str(stop-start) + "s")