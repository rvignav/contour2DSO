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

with open(glob.glob('files/*')[0]) as f: 
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

  p = Path(tupVerts)
  grid = p.contains_points(points)
  mask = grid.reshape(seg_mask.shape[0],seg_mask.shape[1])

  for r in range(seg_mask.shape[0]):
    for c in range(seg_mask.shape[1]):
      seg_mask[r][c][z] = mask[r][c]
  visualize(mask, 'vis/img' + str(z) + '.png')

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

if 'ImageOrientationPatient' in info:
  ds = Dataset()
  ds2 = Dataset()
  ds2.ImageOrientationPatient = info.ImageOrientationPatient
  ds.SliceThickness = info.SliceThickness
  ds.PixelSpacing = info.PixelSpacing
  pos = [ds2]
  pms = [ds]
  ds3 = Dataset()
  ds3.ReferencedSegmentNumber = 1
  sis = [ds3]
  fin = Dataset()
  fin.SegmentIdentificationSequence = sis
  fin.PlaneOrientationSequence = pos
  fin.PixelMeasuresSequence = pms
  info_mask.SharedFunctionalGroupsSequence = [fin]
else:
    if 'ImageOrientationPatient' in info.PerFrameFunctionalGroupsSequence[0].PlaneOrientationSequence[0]:
      ds = Dataset()
      ds.ImageOrientationPatient = info.PerFrameFunctionalGroupsSequence[0].PlaneOrientationSequence[0].ImageOrientationPatient
      ds2 = Dataset()
      ds2.SliceThickness=info.PerFrameFunctionalGroupsSequence[0].PixelMeasuresSequence[0].SliceThickness
      ds2.PixelSpacing=info.PerFrameFunctionalGroupsSequence[0].PixelMeasuresSequence[0].PixelSpacing
      pos = [ds]
      pms = [ds2]
      ds3 = Dataset()
      ds3.ReferencedSegmentNumber = 1
      sis = [ds3]
      fin = Dataset()
      fin.SegmentIdentificationSequence = sis
      fin.PlaneOrientationSequence = pos
      fin.PixelMeasuresSequence = pms
      info_mask.SharedFunctionalGroupsSequence = [fin]

if len(seriesDCM) > 1:
    ris = []
    pffgs = []
    for i in range(startIndex, endIndex+1):
        slice_info=seriesDCM[i]
        ib1=i-startIndex+1

        ds = Dataset()
        ds.ReferencedSOPClassUID=slice_info.SOPClassUID
        ds.ReferencedSOPInstanceUID=slice_info.SOPInstanceUID
        ris.append(ds)

        ds2 = Dataset()
        ds2.ReferencedSOPClassUID=slice_info.SOPClassUID
        ds2.ReferencedSOPInstanceUID=slice_info.SOPInstanceUID
        ds3 = Dataset()
        ds3.CodeValue='121322'
        ds3.CodingSchemeDesignator='DCM'
        ds3.CodeMeaning='Source image for image processing operation'
        ds2.PurposeOfReferenceCodeSequence = [ds3]
        da = Dataset()
        da.SourceImageSequence = [ds2]
        subdcs = Dataset()
        subdcs.CodeValue='113076'
        subdcs.CodingSchemeDesignator='DCM'
        subdcs.CodeMeaning='Segmentation'
        da.DerivationCodeSequence = [subdcs]
        di = Dataset()
        di.DerivationImageSequence = [da]
        fcs = Dataset()
        fcs.StackID='1'
        fcs.InStackPositionNumber=ib1
        # fcs.DimensionIndexValues= [[1],[ib1],[1]]
        di.FrameContentSequence = [fcs]
        
        if 'ImagePositionPatient' in info:
          pps = Dataset()
          pps.ImagePositionPatient=slice_info.ImagePositionPatient
          di.PlanePositionSequence = [pps]
        else:
            if 'ImagePositionPatient' in info.PerFrameFunctionalGroupsSequence[0].PlanePositionSequence[0]:
                print('shouldnt come here. why the information is in frames and it is not a multiframe')
                pps = Dataset()
                pps.ImagePositionPatient=slice_info.PerFrameFunctionalGroupsSequence[0].PlanePositionSequence[0].ImagePositionPatient
                di.PlanePositionSequence = [pps]
            else:
                print('shouldnt happen. why the information is not there and it is not a multiframe')
        pffgs.append(di)

    ds = Dataset()
    ds.ReferencedInstanceSequence = ris
    ds.SeriesInstanceUID=info.SeriesInstanceUID
    info_mask.ReferencedSeriesSequence = [ds]
    info_mask.PerFrameFunctionalGroupsSequence = pffgs
else:
    numOfFrames=len(seriesDCM)
    if numOfFrames<info.NumberOfFrames:
        print('The input mask is smaller than the frames! Assuming they start from the beginning')

    da = Dataset()
    da.ReferencedSOPClassUID=info.SOPClassUID
    da.ReferencedSOPInstanceUID=info.SOPInstanceUID
    ris = [da]
    rss = Dataset()
    rss.ReferencedInstanceSequence = ris
    rss.SeriesInstanceUID=info.SeriesInstanceUID
    info_mask.ReferencedSeriesSequence = [rss]
    
    pffgs = []
    for i in range(1, numOfFrames+1):
        ds = Dataset()
        ds.ReferencedSOPClassUID=info.SOPClassUID
        ds.ReferencedSOPInstanceUID=info.SOPInstanceUID
        ds2 = Dataset()
        ds2.CodeValue='121322'
        ds2.CodingSchemeDesignator='DCM'
        ds2.CodeMeaning='Source image for image processing operation'
        ds.PurposeOfReferenceCodeSequence = [ds2]
        sis = [ds]
        dis = Dataset()
        dis.SourceImageSequence = [sis]

        dcs = Dataset()
        dcs.CodeValue='113076'
        dcs.CodingSchemeDesignator='DCM'
        dcs.CodeMeaning='Segmentation'
        dis.DerivationCodeSequence = [dcs]

        di = Dataset()
        di.DerivationImageSequence = [dis]
        fcs = Dataset()
        fcs.StackID='1'
        fcs.InStackPositionNumber=i
        # fcs.DimensionIndexValues= [[1],[i],[1]]
        di.FrameContentSequence = [fcs]       
        
        if 'ImagePositionPatient' in info.PerFrameFunctionalGroupsSequence[0].PlanePositionSequence[0]:
          ipp = Dataset()
          ipp.ImagePositionPatient=info.PerFrameFunctionalGroupsSequence[0].PlanePositionSequence[0].ImagePositionPatient
          di.PlanePositionSequence = [ipp]
        else:
            print('shouldnt happen. why the information is not there ')
        
        pffgs.append(di)
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
info_mask.InstanceCreatorUID='1.2.276.0.7230010.3'
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

da1 = Dataset()
da2 = Dataset()
da3 = Dataset()
da4 = Dataset()
da4.DimensionOrganizationUID= dicomuid
da1.DimensionIndexPointer=[32, 36950]
da1.FunctionalGroupPointer=[32, 37137]
da1.DimensionDescriptionLabel='Stack ID'
da2.DimensionIndexPointer=[32, 36951]
da2.FunctionalGroupPointer=[32, 37137]
da2.DimensionDescriptionLabel='In-Stack Position Number'
da3.DimensionIndexPointer=[98, 11]
da3.FunctionalGroupPointer=[98,10]
da3.DimensionDescriptionLabel='Referenced Segment Number'

info_mask.DimensionIndexSequence = [da1,da2,da3]
info_mask.DimensionOrganizationSequence = [da4]

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

ds = Dataset()
ds.CodeValue = 'T-D0050'
ds.CodingSchemeDesignator='SRT'
ds.CodeMeaning='Tissue'
ana = [ds]

ds2 = Dataset()
ds2.CodeValue = 'T-D0050'
ds2.CodingSchemeDesignator='SRT'
ds2.CodeMeaning='Tissue'
segseq = [ds2]

ds3 = Dataset()
ds3.CodeValue = 'T-D0050'
ds3.CodingSchemeDesignator='SRT'
ds3.CodeMeaning='Tissue'
segseq2 = [ds3]

dataset = Dataset()
dataset.AnatomicRegionSequence = ana
dataset.SegmentedPropertyCategoryCodeSequence = segseq
dataset.SegmentNumber = 1
dataset.SegmentLabel='Segmentation'
dataset.SegmentAlgorithmType='SEMIAUTOMATIC'
dataset.SegmentAlgorithmName='ePAD'
dataset.SegmentedPropertyTypeCodeSequence=segseq2
seg = [dataset]
info_mask.SegmentSequence = seg

info_mask.ContentCreatorsName='ePAD^matlab'
info_mask.ContentLabel= 'ROI'

info_mask.ContentDescription=str(name) + 'segmentation'
info_mask.SeriesDescription=str(name) + ' segmentation'

file_name='output/' + str(name) + '.dcm'
if (os.path.exists(file_name)):
  file_name += currentTime

np_frame = np.array(seg_mask,dtype=np.uint8)
info_mask.PixelData = np_frame.tobytes()

info_mask.save_as(file_name)
print("info_mask saved to " + file_name)