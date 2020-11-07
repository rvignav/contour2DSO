FROM python:3.5

RUN apt-get update -y
RUN apt-get install -y libgl1-mesa-glx
RUN apt install octave
RUN pip install matplotlib opencv-python pydicom numpy argparse Pillow glob2

ADD createDSO.py /

ENTRYPOINT [ "python", "./createDSO.py" ]
