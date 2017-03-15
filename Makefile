# Cuda SDK
CUDASDKDIR = /opt/nvidia/cuda4/NVIDIA_GPU_Computing_SDK/C
CUDASDKINC = $(CUDASDKDIR)/common/inc
CUDALIB = /opt/nvidia/cuda4/lib64
CUDASDKLIB = $(CUDASDKDIR)/lib

# CUDA-specific libraries from runtime and SDK
CUINCFLAGS =  -I$(CUDASDKINC)
CULIBFLAGS =  -L/usr/local/cuda-6.5/lib64 -L$(CUDASDKLIB) -L${CUDASDKDIR}/common/lib/linux -lcuda -L${CUDALIB} -lcuda -lcudadevrt -lcuda -lcudart

CUFLAGS = 
CCFLAGS = -O3  -openmp -std=c99

SEP=/data/bob/SEP
SEPINC=${SEP}/include
SEPLIB=${SEP}/lib
SEPFLAG=-I${SEPINC} -L${SEPLIB} -L$(SEPLIB) -lsepauxf90 -lsepgeef90 -lsepf90 -lsep2df90 -lsep3df90 -lsep3d -lsepf90 -lsep

BOBCPP= /home/clapp/base
BOBFLAG= -L${BOBCPP} -lbob_cpp
CXX  = /lib/cpp

# Common Libraries
CMNINCFLAGS = -I. -I$(SEPINC) -I$(MPIINC) -I$(FFTWINC) -I$(S)
CMNLIBFLAGS = $(SEPLIBFLAGS) $(MPILIBFLAGS) $(FFTWLIBFLAGS) $(GOMPLIBFLAGS) -lm -lpthread -lgomp

# Fortran Includes
F90INCFLAGS  = $(CMNINCFLAGS) $(CUINCFLAGS)
F90LIBFLAGS  = $(CMNLIBFLAGS) $(CULIBFLAGS)
F90FLAGS     = -O3 -openmp -FR -g -traceback -module $(M) -static-intel -check bounds
F90  = ifort

NVCC = /usr/local/cuda-6.5/bin/nvcc   $(CUFLAGS) $(CUINCFLAGS) $(CULIBFLAGS)  -arch=sm_20 -Xptxas -v,-abi=no,-dlcm=cg -ftz=true -maxrregcount=32

RESDIR       = ./Fig
BINDIR       = ./Bin
OBJDIR       = ./Obj
SRCDIR       = ./Src
R            = ${RESDIR}
B            = ${BINDIR}
O            = ${OBJDIR}
P            = ${PARDIR}
D            = ${DATDIR}
S            = ${SRCDIR}

CU_OBJS_3D= gpu_funcs_3d.o

CPP_OBJS_3D= vel_fd_3d.o  data_rtm_3d.o  image_rtm_3d.o  wavefield_insert_3d.o  source_func_3d.o  rtm_zero_op_3d.o  deriv_3d.o laplac_3d.o fd_prop_3d.o  map_data_3d.o  hypercube_float.o  my_operator.o  hypercube.o   tmute.o cgstep.o lin_solver.o i_op.o super_vector.o combo_oper.o float_3d.o  oc_float.o  sreg.o  sregf.o  axis.o  sinc_bob.o  sep_params.o  param_func.o



%.o: %.cu
	${NVCC} -g $*.cu -c -o $@ $(CUFLAGS) $(CUINCPATH) $(CULIBPATH) $(CULIBS)  ${SEPFLAG}

cgstep.o: cgstep.cpp
	g++ -Wall -o $*.o  -c cgstep.cpp  ${SEPFLAG} -I. -I.  -O2 -g 

%.o: %.cpp
	g++ -Wall -o $*.o  -c $*.cpp ${SEPFLAG} -I. -I.  -O2 -g 


%.o: %.f90
	# Preprocessing FORTRAN $^
	$(CXX) $(F90INCFLAGS) -DSOURCE="'$^'" < $*.f90  > $*.fix.f90
	# Compiling FORTRAN $@
	$(F90) $(F90INCFLAGS) $(F90FLAGS)  -c -o $*.o $*.fix.f90

FIG          = >/dev/null out=

%.x: %.o
	#Linking FORTRAN $@
	${F90} $(F90FLAGS) -o $@ $^       ${SEPFLAG}

RTM2_3d.x: ${CU_OBJS_3D} ${CPP_OBJS_3D}   RTM2_3d.o
	g++  ${CULIBFLAGS} ${CULIBS} RTM2_3d.o  ${CPP_OBJS_3D} ${CU_OBJS_3D}  -o $@ ${CULIBFLAGS} ${CULIBS}  ${SEPFLAG}

RTM_fwd.x: ${CU_OBJS_3D} ${CPP_OBJS_3D}   RTM_fwd.o
	nvcc  ${CULIBFLAGS} ${CULIBS} RTM_fwd.o  ${CPP_OBJS_3D} ${CU_OBJS_3D}  -o $@ ${CULIBFLAGS} ${CULIBS}  ${SEPFLAG}

RTM_inv.x: ${CU_OBJS_3D} ${CPP_OBJS_3D}   RTM_inv.o
	nvcc   ${CULIBFLAGS} ${CULIBS}  ${CPP_OBJS_3D} ${CU_OBJS_3D}  -o $@ ${CULIBFLAGS} ${CULIBS}  ${SEPFLAG} RTM_inv.o 

wavelet.H: wave.P
	Wavelet par=wave.P |Window j1=4 >$@

data_single.H: 
	Spike n1=1800 o1=0 d1=.004 o2=0 d2=50 o3=0 d3=50 n2=208 n3=208 n4=9 n5=7 o5=3650 o4=3650 d4=400 d5=400 >$@

data2.H: 
	Spike n1=2500 o1=0 d1=.004 o2=0 d2=50 o3=0 d3=50 n2=208 n3=208 n4=9 n5=7 o5=3650 o4=3650 d4=400 d5=400 >$@

vel_t.H: mandy_go.H
	Transp maxsize=2000 < mandy_go.H plane=12 |Transp plane=23 >$@ maxsize=2000
	
refl_t.H: refl_go.H
	Transp maxsize=2000 < refl_go.H plane=12 |Transp plane=23 >$@ maxsize=2000
vel2_t.H: mandy_mo.H
	Transp maxsize=2000 < mandy_mo.H plane=12 |Transp plane=23 >$@ maxsize=2000
	
refl2_t.H: refl_mo.H
	Transp maxsize=2000 < refl_mo.H plane=12 |Transp plane=23 >$@ maxsize=2000


test2.H: RTM_fwd.x refl2_t.H vel2_t.H wavelet.H data2.H
	./RTM_fwd.x wavelet=wavelet.H image=refl2_t.H velocity=vel2_t.H src_depth=0 data=data2.H dataout=$@ n_gpus=4 rec_depth=900

test.H: RTM_fwd.x refl_t.H vel_t.H wavelet.H data_single.H
	./RTM_fwd.x wavelet=wavelet.H image=refl_t.H velocity=vel_t.H src_depth=900 data=data_single.H dataout=$@ n_gpus=4


wave.seam3d.H:
	${SEP}/bin/Wavelet n1=8192 o1=0 d1=0.004 wavelet=ricker2 phase=0 fund=20 tdelay=0.075 | ${SEP}/bin/Window3d n1=128 > tmp.H
	Math file1=tmp.H exp='file1*1000' > $@
	echo "label1='Time(s)'"  >> $@

VEL_DATA = /data/seam/earthmodel/Vp_xyz_10m.bin

seam_full_vel.H: ${VEL_DATA}
	echo 	o1=0 o2=0 o3=0  >>$@
	echo	n1=3501 n2=4001 n3=1501  >>$@
	echo	d1=10 d2=10 d3=10  >>$@
	echo	data_format=""xdr_float""  >>$@
	echo	label1="x(m)" label2="y(m)" label3="z(m)" >>$@
	echo	esize= 4  >>$@
	echo in="$<" >> $@

seam_vel.H: seam_full_vel.H
#	Window3d f1=1000 n1=1650 f2=2220 n2=1650 n3=725 < $<  > $*3.H
	Window3d f1=300 n1=1650 f2=2220 n2=1650 n3=725 < $<  > $*3.H
	Transp plane=12 reshape=2,3,4 < $*3.H > $*4.H
	cp $*4.H $@
	swapbytes < /scr1/chrisl/$*4.H@ > /scr1/chrisl/$@@
	echo in="/scr1/chrisl/$@@" >> $@
	Rm $*4.H
	Rm $*3.H

seam_small.H: seam_vel.H
#	Window3d n1=500 n2=800 f1=0 f2=200 < $< > $@ 
	Window3d n1=724 n3=300 < $< > tmps.H
	Window3d j1=2 j2=2 < tmps.H > $@
	Rm tmps.H
	echo "o2=0 o3=0 d1=10 d2=10" >> $@

%.borndata.H: RTM_fwd.x wave.seam3d.H window.x
	Window3d n2=800 < $*.H | Transp plane=12 reshape=1,3 > tmp_vel.H
	Pad beg3=20 extend=1 < tmp_vel.H > tmp_vel2.H
	Gpow gpow=-1 < tmp_vel2.H | Smooth rect1=10 rect2=10 rect3=10 | Gpow gpow=-1 > $*.smooth.H
	Math file1=$*.smooth.H file2=tmp_vel2.H exp='file2-file1' > $*.refl.H #calculate reflectivity
	window.x < tmp_vel2.H cut=1501 | Smooth rect1=8 rect2=8 rect3=8 > $*.mute.H #mute for water column
	Transp plane=12 reshape=2,3 < $*.refl.H > $*.dims.H #right size for dimensions (geometry) file
	echo "n1=1250 d1=0.004 o1=0. n4=20 o4=2000 d4=250 n5=5 o5=500 d5=500" >> $*.dims.H #put in geometry here
	RTM_fwd.x n_gpus=8 wavelet=wave.seam3d.H image=$*.refl.H velocity=$*.smooth.H data=$*.dims.H dataout=$@ bc_a=40 bc_b=0.0005 source_fields=srcfld.H src_depth=200. mute=$*.mute.H

