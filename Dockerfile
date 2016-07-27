FROM ubuntu:14.04

MAINTAINER Apaar Shanker <apaar92@gmail.com>

USER root

RUN apt-get update
RUN apt-get install -y build-essential
#RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test
RUN apt-get install -y libstdc++6
# Install python
RUN apt-get install -y python
RUN apt-get install -y build-essential python-dev
RUN apt-get install -y python-pip
#RUN apt-get install -y octave

#Install Python Dependencies
ADD requirements.txt requirements.txt
RUN pip install -r requirements.txt

RUN conda install libgcc
ADD Asap-3.8.4 Asap-3.8.4
RUN cd Asap-3.8.4 \
                  && make depend \
                  && make serial \
                  && make install

#CMD octave -q compiling.m
#CMD python Asap-3.8.4/setup.py install

RUN ls
RUN pwd
  
ENV export PATH=/home/main/notebooks/Asap-3.8.4/x86_64:$PATH
ENV export PYTHONPATH=/home/main/notebooks/Asap-3.8.4/Python:/home/main/notebooks/Asap-3.8.4/x86_64:$PYTHONPATH

#RUN export PATH=/home/main/notebooks/Asap-3.8.4/Python:$PATH
#RUN export PYTHONPATH=/home/main/notebooks/Asap-3.8.4/Python:$PYTHONPATH
