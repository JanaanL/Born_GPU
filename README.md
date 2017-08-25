
Pre-reqs 

You need  a c++ compiler with  c++11 support and tbb.

You need to install two small libraries before installing Born.

Install instructions

RUN git clone http://zapad.Stanford.EDU/bob/hypercube.git   /opt/hypercube/src && \
    mkdir /opt/hypercube/build && \
    cd /opt/hypercube/build &&\
    cmake -DCMAKE_INSTALL_PREFIX=/opt/hypercube ../src &&\
    make install &&\
    rm -rf /opt/hypercube/build
RUN git clone http://zapad.Stanford.EDU/bob/genericIO.git /opt/genericIO/src && \
    mkdir /opt/genericIO/build &&\
    cd /opt/genericIO/build &&\
    cmake  -Dhypercube_DIR=/opt/hypercube/lib  -DCMAKE_INSTALL_PREFIX=/opt/genericIO ../src &&\
    make install && \
    rm -rf /opt/genericIO/build
RUN git clone http://zapad.Stanford.EDU/SEP-external/Born.git /opt/born/src && \
  mkdir /opt/born/build && \
  cd /opt/born/build &&  \
    cmake  -Dhypercube_DIR=/opt/hypercube/lib -DgenericIO_DIR=/opt/genericIO/lib  -DCMAKE_INSTALL_PREFIX=/opt/genericIO ../src &&\
    make install &&\
    rm -rf /opt/born/build



