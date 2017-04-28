#include <ioModes.h>
using namespace SEP;

int main(int argc, char **argv){

	ioModes modes(argc,argv);



	std::shared_ptr<genericIO>  io=modes.getDefaultIO();
	std::shared_ptr<paramObj> par=io->getParamObj();





	std::string out=par->getString(std::string("out"));

	std::vector<SEP::axis> axes;
	int nt=par->getInt("nt",750); float dt=par->getFloat("dt",.008); axes.push_back(SEP::axis(nt,0.,dt));
	int nx=par->getInt("nx",300);
	float dx=par->getFloat("dx",10.);
	float ox=par->getFloat("ox",500.);
	axes.push_back(SEP::axis(nx,ox,dx));
	int ny=par->getInt("ny",300);
	float dy=par->getFloat("dy",10.);
	float oy=par->getFloat("oy",500.);
	axes.push_back(SEP::axis(ny,oy,dy));
        axes.push_back(SEP::axis(1,ox+dx*nx/2,dx));
        axes.push_back(SEP::axis(1,oy+dy*ny/2,dy));
 

	std::shared_ptr<hypercube> hyp(new hypercube(axes));

	std::shared_ptr<genericRegFile> outp=io->getRegFile(out,usageOut);
	outp->setHyper(hyp);
	outp->writeDescription();

	std::vector<float> val(hyp->getN123(),0.);

	val[nt/2+nx*nt/2+nx*nt*ny/2]=1.;



	outp->writeFloatStream(val.data(),hyp->getN123());


}
