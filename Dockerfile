# Dockerfile
# Contains commands used to create an automated build for createDSO.py.
# Running the commands
# 
#      docker build -t run .
#      docker run -it -v "local/aims:/home/series/files" -v "local/patient/series:/home/series/PatientSeries" run "name of AIM file" "filename of series"
# 
# will create a new Docker container and run createDSO.py within that container,
# resulting in a DICOM segmentation object stored in the 'output' directory.
# 
# MIT License
#
# Copyright (c) 2020 Vignav Ramesh, ePAD Team (Stanford University)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

FROM python:3.5

RUN apt-get update -y
RUN apt-get install -y libgl1-mesa-glx
RUN pip install matplotlib opencv-python pydicom numpy argparse Pillow glob2

ADD createDSO.py /

ENTRYPOINT [ "python", "./createDSO.py" ]