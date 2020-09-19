FROM python:3.5

ADD dso_matlab /

RUN mkdir /mcr-install && \
    mkdir /opt/mcr && \
    cd /mcr-install && \
    wget -q http://ssd.mathworks.com/supportfiles/downloads/R2016b/deployment_files/R2016b/installers/glnxa64/MCR_R2016b_glnxa64_installer.zip && \
    unzip -q MCR_R2016b_glnxa64_installer.zip && \
    rm -f MCR_R2016b_glnxa64_installer.zip && \
    ./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    cd / && \
    rm -rf mcr-install

ENV LD_LIBRARY_PATH /opt/mcr/v91/runtime/glnxa64:/opt/mcr/v91/bin/glnxa64:/opt/mcr/v91/sys/os/glnxa64
ENV XAPPLRESDIR /opt/mcr/v91/X11/app-defaults

RUN cd dso_matlab && \
    python setup.py install

ADD createDSO.py /
RUN apt-get update -y
RUN apt-get install -y libgl1-mesa-glx
RUN pip install opencv-python pydicom numpy argparse Pillow glob2 oct2py
ENTRYPOINT [ "python", "./createDSO.py" ]