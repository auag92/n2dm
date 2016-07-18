FROM ubuntu:14.04

MAINTAINER Apaar Shanker <apaar92@gmail.com>

USER root

# Install python dependencies
RUN apt-get update
RUN apt-get install -y python
RUN apt-get install -y build-essential python-dev
RUN apt-get install -y python-pip
#RUN apt-get install -y python-numpy
#RUN apt-get install -y python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose

#USER main
# Install requirements for Python 2

ADD requirements.txt requirements.txt
RUN pip install -r requirements.txt
