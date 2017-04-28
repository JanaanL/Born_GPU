#include <vel_fd_3d.h>
#include <data_rtm_3d.h>
#include <image_rtm_3d.h>
#include <ioModes.h>
#include "source_func_3d.h"
#include "rtm_zero_op_3d.h"
#include "cpu_prop.h"

main(int argc, char **argv){


	SEP::ioModes modes(argc,argv);



	std::shared_ptr<SEP::genericIO>  io=modes.getDefaultIO();
	std::shared_ptr<SEP::paramObj> pars=io->getParamObj();
	fprintf(stderr,"becore data\n");
	std::shared_ptr<data_rtm_3d> data(new data_rtm_3d("data",io));
	fprintf(stderr,"before waelet \n");
	std::shared_ptr<wavelet_source_func> wavelet(new wavelet_source_func(io,"wavelet"));
	fprintf(stderr,"after wavelet \n");
	float src_depth=pars->getFloat("src_depth",0.);
	fprintf(stderr,"before velicty \n");
	std::shared_ptr<vel_fd_3d> vel(new vel_fd_3d(io,"velocity"));
	fprintf(stderr,"after velicity \n");
	std::shared_ptr<image_rtm_3d> image(new image_rtm_3d(io,"image",vel));
	fprintf(stderr,"after image \n");
	SEP::axis asx=data->getAxis(4);
	SEP::axis asy=data->getAxis(5);


	wavelet->set_sources_axes(src_depth,asx,asy);


	image->zero();
	//image->set_source_file(image);

	//std::shared_ptr<oc_float> d2(new oc_float(io,"dtest.H",data));
	//data->set_source_file(d2);


	//int n_gpus=pars.get_int("n_gpus",1);
	//setup_cuda(n_gpus,argc,argv);

	float aper=pars->getFloat("aper",8.);

	std::vector<int> rand_vec(1.0);

	bool encode=false;
	std::shared_ptr<cpuProp> prop(new cpuProp(io));


	std::shared_ptr<rtm_zero_op> op(new
		rtm_zero_op(pars,prop,vel,wavelet,data,image,aper,true,encode,rand_vec,true));



	op->adjoint(false,image,data,1);

	return 0;
}
