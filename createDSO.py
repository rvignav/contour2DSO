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
parser.add_argument('series_path', type=str, help='Path to series')

args = parser.parse_args()

series_path = args.series_path

start = timeit.default_timer()

studies = glob.glob('SamplePatient/*')
# studies = glob.glob('/home/series/PatientSeries/*')
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

with open(glob.glob('files/*')[0]) as f: #/home/series/files/*
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
  
  x, y = np.meshgrid(np.arange(seg_mask.shape[0]), np.arange(seg_mask.shape[1]))
  x, y = x.flatten(), y.flatten()
  points = np.vstack((x,y)).T 

  p = Path(tupVerts) # make a polygon
  grid = p.contains_points(points)
  mask = grid.reshape(seg_mask.shape[0],seg_mask.shape[1])

  for r in range(seg_mask.shape[0]):
    for c in range(seg_mask.shape[1]):
      seg_mask[r][c][z] = mask[r][c]

startIndex = 0
endIndex = 0

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

print("Seg_mask array created")

stop = timeit.default_timer()
print("Runtime (mask creation): " + str(stop-start) + "s")

file_meta = FileMetaDataset()

file_meta.MediaStorageSOPClassUID='1.2.840.10008.5.1.4.1.1.66.4'
dicomuid = generate_uid()
instanceuid=dicomuid
file_meta.MediaStorageSOPInstanceUID=instanceuid
file_meta.TransferSyntaxUID=ExplicitVRLittleEndian
file_meta.ImplementationClassUID='1.2.840.10008.5.1.4.1.1.66.4'
file_meta.ImplementationVersionName='ePAD_matlab_1.0'

suffix = '.dcm'
info = seriesDCM[0]
info_mask = FileDataset('output/' + str(name) + str(suffix), {},
                  file_meta=file_meta, preamble=b"\0" * 128)
info_mask.StudyDescription=info.StudyDescription

# if 'ImageOrientationPatient' in info:
#     info_mask.SharedFunctionalGroupsSequence[0].PlaneOrientationSequence[0].ImageOrientationPatient=info.ImageOrientationPatient
#     info_mask.SharedFunctionalGroupsSequence[0].PixelMeasuresSequence[0].SliceThickness=info.SliceThickness
#     info_mask.SharedFunctionalGroupsSequence[0].PixelMeasuresSequence[0].PixelSpacing=info.PixelSpacing
# else:
#     if 'ImageOrientationPatient' in info.PerFrameFunctionalGroupsSequence[0].PlaneOrientationSequence[0]:
#         info_mask.SharedFunctionalGroupsSequence[0].PlaneOrientationSequence[0].ImageOrientationPatient=info.PerFrameFunctionalGroupsSequence[0].PlaneOrientationSequence[0].ImageOrientationPatient
#         info_mask.SharedFunctionalGroupsSequence[0].PixelMeasuresSequence[0].SliceThickness=info.PerFrameFunctionalGroupsSequence[0].PixelMeasuresSequence[0].SliceThickness
#         info_mask.SharedFunctionalGroupsSequence[0].PixelMeasuresSequence[0].PixelSpacing=info.PerFrameFunctionalGroupsSequence[0].PixelMeasuresSequence[0].PixelSpacing
# info_mask.SharedFunctionalGroupsSequence[0].SegmentIdentificationSequence[0].ReferencedSegmentNumber=1

# if len(seriesDCM) > 1:
#     for i in range(startIndex, endIndex+1):
#         slice_info=seriesDCM[i]
#         ib1=i-startIndex+1
#         info_mask.ReferencedSeriesSequence[0].ReferencedInstanceSequence[ib1].ReferencedSOPClassUID=slice_info.SOPClassUID
#         info_mask.ReferencedSeriesSequence[0].ReferencedInstanceSequence[ib1].ReferencedSOPInstanceUID=slice_info.SOPInstanceUID
        
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].DerivationImageSequence[0].SourceImageSequence[0].ReferencedSOPClassUID=slice_info.SOPClassUID
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].DerivationImageSequence[0].SourceImageSequence[0].ReferencedSOPInstanceUID=slice_info.SOPInstanceUID
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].DerivationImageSequence[0].SourceImageSequence[0].PurposeOfReferenceCodeSequence[0].CodeValue='121322'
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].DerivationImageSequence[0].SourceImageSequence[0].PurposeOfReferenceCodeSequence[0].CodingSchemeDesignator='DCM'
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].DerivationImageSequence[0].SourceImageSequence[0].PurposeOfReferenceCodeSequence[0].CodeMeaning='Source image for image processing operation'
        
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].DerivationImageSequence[0].DerivationCodeSequence[0].CodeValue='113076'
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].DerivationImageSequence[0].DerivationCodeSequence[0].CodingSchemeDesignator='DCM'
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].DerivationImageSequence[0].DerivationCodeSequence[0].CodeMeaning='Segmentation'
        
        
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].FrameContentSequence[0].StackID='1'
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].FrameContentSequence[0].InStackPositionNumber=ib1
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].FrameContentSequence[0].DimensionIndexValues= [[1],[ib1],[1]]
#         if 'ImagePositionPatient' in info:
#             info_mask.PerFrameFunctionalGroupsSequence[ib1].PlanePositionSequence[0].ImagePositionPatient=slice_info.ImagePositionPatient
#         else:
#             if 'ImagePositionPatient' in info.PerFrameFunctionalGroupsSequence[0].PlanePositionSequence[0]:
#                 print('shouldnt come here. why the information is in frames and it is not a multiframe')
#                 info_mask.PerFrameFunctionalGroupsSequence[ib1].PlanePositionSequence[0].ImagePositionPatient=slice_info.PerFrameFunctionalGroupsSequence[0].PlanePositionSequence[0].ImagePositionPatient
#             else:
#                 print('shouldnt happen. why the information is not there and it is not a multiframe')
# else:
#     numOfFrames=len(seriesDCM)
#     if numOfFrames<info.NumberOfFrames:
#         print('The input mask is smaller than the frames! Assuming they start from the beginning')

#     info_mask.ReferencedSeriesSequence[0].ReferencedInstanceSequence[0].ReferencedSOPClassUID=info.SOPClassUID
#     info_mask.ReferencedSeriesSequence[0].ReferencedInstanceSequence[0].ReferencedSOPInstanceUID=info.SOPInstanceUID
    
    
#     for i in range(1, numOfFrames+1):
#         ib1 = i-1
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].DerivationImageSequence[0].SourceImageSequence[0].ReferencedSOPClassUID=info.SOPClassUID
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].DerivationImageSequence[0].SourceImageSequence[0].ReferencedSOPInstanceUID=info.SOPInstanceUID
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].DerivationImageSequence[0].SourceImageSequence[0].PurposeOfReferenceCodeSequence[0].CodeValue='121322'
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].DerivationImageSequence[0].SourceImageSequence[0].PurposeOfReferenceCodeSequence[0].CodingSchemeDesignator='DCM'
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].DerivationImageSequence[0].SourceImageSequence[0].PurposeOfReferenceCodeSequence[0].CodeMeaning='Source image for image processing operation'
        
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].DerivationImageSequence[0].DerivationCodeSequence[0].CodeValue='113076'
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].DerivationImageSequence[0].DerivationCodeSequence[0].CodingSchemeDesignator='DCM'
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].DerivationImageSequence[0].DerivationCodeSequence[0].CodeMeaning='Segmentation'
        
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].FrameContentSequence[0].StackID='1'
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].FrameContentSequence[0].InStackPositionNumber=i
#         info_mask.PerFrameFunctionalGroupsSequence[ib1].FrameContentSequence[0].DimensionIndexValues= [[1],[i],[1]]
#         if 'ImagePositionPatient' in info.PerFrameFunctionalGroupsSequence[0].PlanePositionSequence[0]:
#             info_mask.PerFrameFunctionalGroupsSequence[ib1].PlanePositionSequence[0].ImagePositionPatient=info.PerFrameFunctionalGroupsSequence[0].PlanePositionSequence[0].ImagePositionPatient
#         else:
#             print('shouldnt happen. why the information is not there ')
# info_mask.ReferencedSeriesSequence[0].SeriesInstanceUID=info.SeriesInstanceUID


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
# info_mask.InstanceCreatorUID='1.2.276.0.7230010.3'
info_mask.SOPClassUID='1.2.840.10008.5.1.4.1.1.66.4'
info_mask.SOPInstanceUID= instanceuid

info_mask.AccessionNumber=info.AccessionNumber
info_mask.Modality='SEG'
info_mask.Manufacturer='Stanford University'

info_mask.ManufacturerModelName= 'ePAD Matlab'
info_mask.DeviceSerialNumber='SN123456'
info_mask.SoftwareVersion='1.0'
info_mask.StudyInstanceUID=info.StudyInstanceUID
info_mask.SeriesInstanceUID= dicomuid
info_mask.SeriesNumber= 1000
info_mask.ContentDate=datetime.today().strftime('%Y%m%d')
info_mask.StudyDate=info.StudyDate
info_mask.SeriesDate=datetime.today().strftime('%Y%m%d')
info_mask.AcquisitionDate=datetime.today().strftime('%Y%m%d')
currentTime=datetime.today().strftime('%H%m%s.%f')
info_mask.ContentTime=currentTime
info_mask.StudyTime=info.StudyTime
info_mask.SeriesTime=currentTime
info_mask.AcquisitionTime=currentTime
info_mask.InstanceNumber= 1
info_mask.FrameOfReferenceUID= info.FrameOfReferenceUID
info_mask.PositionReferenceIndicator= ''

# info_mask.DimensionOrganizationSequence[0].DimensionOrganizationUID= dicomuid
# info_mask.DimensionIndexSequence[0].DimensionIndexPointer=[32, 36950]
# info_mask.DimensionIndexSequence[0].FunctionalGroupPointer=[32, 37137]
# info_mask.DimensionIndexSequence[0].DimensionDescriptionLabel='Stack ID'
# info_mask.DimensionIndexSequence[2].DimensionIndexPointer=[32, 36951]
# info_mask.DimensionIndexSequence[2].FunctionalGroupPointer=[32, 37137]
# info_mask.DimensionIndexSequence[2].DimensionDescriptionLabel='In-Stack Position Number'
# info_mask.DimensionIndexSequence[3].DimensionIndexPointer=[98, 11]
# info_mask.DimensionIndexSequence[3].FunctionalGroupPointer=[98,10]
# info_mask.DimensionIndexSequence[3].DimensionDescriptionLabel='Referenced Segment Number'
info_mask.SamplesPerPixel= 1
info_mask.PhotometricInterpretation= 'MONOCHROME2'
info_mask.NumberOfFrames= len(seriesDCM)
info_mask.Rows= seg_mask.shape[0]
info_mask.Columns= seg_mask.shape[1]

info_mask.BitsAllocated= 8
info_mask.BitsStored= 8
info_mask.HighBit= 7

info_mask.PixelRepresentation= 0
info_mask.LossyImageCompression='00'
info_mask.SegmentationType='BINARY'

# info_mask.SegmentSequence[0].AnatomicRegionSequence[0].CodeValue='T-D0050'
# info_mask.SegmentSequence[0].AnatomicRegionSequence[0].CodingSchemeDesignator='SRT'
# info_mask.SegmentSequence[0].AnatomicRegionSequence[0].CodeMeaning='Tissue'

# info_mask.SegmentSequence[0].SegmentedPropertyCategoryCodeSequence[0].CodeValue='T-D0050'
# info_mask.SegmentSequence[0].SegmentedPropertyCategoryCodeSequence[0].CodingSchemeDesignator='SRT'
# info_mask.SegmentSequence[0].SegmentedPropertyCategoryCodeSequence[0].CodeMeaning='Tissue'
# info_mask.SegmentSequence[0].SegmentNumber=1
# info_mask.SegmentSequence[0].SegmentLabel='Segmentation'
# info_mask.SegmentSequence[0].SegmentAlgorithmType='SEMIAUTOMATIC'
# info_mask.SegmentSequence[0].SegmentAlgorithmName='ePAD'

# info_mask.SegmentSequence[0].SegmentedPropertyTypeCodeSequence[0].CodeValue='T-D0050'
# info_mask.SegmentSequence[0].SegmentedPropertyTypeCodeSequence[0].CodingSchemeDesignator='SRT'
# info_mask.SegmentSequence[0].SegmentedPropertyTypeCodeSequence[0].CodeMeaning='Tissue'

info_mask.ContentCreatorsName='ePAD^matlab'
info_mask.ContentLabel= 'ROI'

info_mask.ContentDescription=str(name) + 'segmentation'
info_mask.SeriesDescription=str(name) + ' segmentation'

file_name='output/' + str(name) + '.dcm' # file_name='/home/output/' + str(name) + '.dcm'
if (os.path.exists(file_name)):
  file_name += currentTime

np_frame = np.array(seg_mask,dtype=np.uint8)
info_mask.PixelData = np_frame.tobytes()

info_mask.save_as(file_name)
print("info_mask saved to " + file_name)