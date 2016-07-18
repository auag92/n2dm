FROM ouruser/ubuntu:py1

MAINTAINER Apaar Shanker <apaar92@gmail.com>

USER root

RUN apt-get update
# Install python
RUN apt-get install -y python
RUN apt-get install -y build-essential python-dev
RUN apt-get install -y python-pip

#Install Python Dependencies
ADD requirements.txt requirements.txt
RUN pip install -r requirements.txt

ADD Asap-3.8.4 Asap-3.8.4

RUN python Asap-3.8.4/setup.py install
RUN export PATH=/Asap-3.8.4/x86_64:$PATH
RUN export PYTHONPATH=/Asap-3.8.4/Python:$HOME/Asap/x86_64:$PYTHONPATH
