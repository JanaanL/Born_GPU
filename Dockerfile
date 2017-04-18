From centos:7
MAINTAINER Bob Clapp <bob@sep.stanford.edu>
RUN yum -y install gcc-gfortran gcc-c++ cmake make git
RUN git clone http://zapad.Stanford.EDU/bob/genericIO.git /opt/genericIO
RUN mkdir -p /opt/genericIO/build
RUN cd /opt/genericIO/build; cmake -DCMAKE_INSTALL_PREFIX=/opt/genericIO  ..; make install;
