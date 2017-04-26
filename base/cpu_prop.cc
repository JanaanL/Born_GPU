#include "cpu_prop.h"
#define C0  0
#define CZ1 1
#define CX1 2
#define CY1 3
#define CZ2 4
#define CX2 5
#define CY2 6
#define CZ3 7
#define CX3 8
#define CY3 9
#define CZ4 10
#define CX4 11
#define CY4 12

cpuProp::cpuProp(std::shared_ptr<SEP::genericIO> io){
	storeIO(io);

}

void cpuProp::rtmForward(int n1, int n2, int n3, int jt, float *img,
	float *rec, int npts, int nt, int nt_big, int rec_nx, int rec_ny){

//	rtm_forward(n1,n2,n3,jt,img,rec,npts,nt,nt_big,rec_nx,rec_ny);

}
void cpuProp::rtmAdjoint(int n1, int n2, int n3, int jtd, float *src_p0, float *src_p1,
	float *img, int npts_s, int nt){
//    rtm_adjoint(ad1.n,ad2.n,ad3.n,jtd,src_p0->vals,src_p1->vals,img->vals,npts_s,nt/*,src,recx*/);
	std::vector<float> rec_p0(_n123,0.),rec_p1(_n123,0.);
	float *r_p0=rec_p0.data(), *r_p1=rec_p1.data();
	_dir=-1;
	for(int it=nt-1; it >=0; it--) {
		int id_s=(it+1)/_jtsS;
		int i_s=it+1-id_s*_jtsS;
		int id=it/_jtdD;
		int ii=it-id*_jtdD;
		if(it>0) prop(rec_p0.data(),rec_p1.data(),_vel2);
		if(it< nt-1) {
			prop(src_p0,src_p1,_vel1);
			injectSource(id_s,i_s,src_p1);                                                                                                                                                                                                                                                                                                                                                                                                                               //I think this should be p0;
		}
		damp(r_p0,r_p1);
		injectReceivers(id,ii,r_p0);
		imageCondition(r_p0,src_p0,img);
		float *pt=src_p0;
		src_p0=src_p1;
		src_p1=pt;
		pt=r_p0;
		r_p0=r_p1;
		r_p1=pt;
	}

}
void cpuProp::imageCondition(float *rec, float *src, float *img){


	for(long long i=0; i < _n123; i++) {
		img[i]+=src[i]*rec[i];
	}



}
void cpuProp::sourceProp(int nx, int ny, int nz, bool damp, bool getLast,
	float *p0, float *p1, int jts, int npts, int nt){


	int n12=_nx*_ny;
	_dir=1;
	for(int it=0; it<=nt; it++) {
		int id=it/_jtsS;
		int ii=it-id*_jtsS;
		prop(p0,p1,_vel1);
		injectSource(id,ii,p0);
		float *pt=p1; p1=p0; p1=pt;
	}

}

void cpuProp::damp(float *p0,float *p1){

	for(int i3=4; i3 <_nz-4; i3++) {
		int edge=std::min(i3-4,_nz-4-i3);
		for(int i2=4; i2 < _ny-4; i2++) {
			edge=std::min(edge,std::min(i2-4,_ny-4-i2));
			int ii=i2*_nx+4+_n12*i3;
			for(int i1=4; i1 < _nx-4; i1++,ii++) {
				edge=std::min(edge,std::min(i1-4,_nx-4-i1));
				if(edge>=0) {
					p0[ii]*=_bound[edge];
					p1[ii]*=_bound[edge];
				}
			}
		}
	}

}
void cpuProp::injectSource(int id, int ii, float *p){


	for(int i=0; i < _nptsS; i++) {
		p[_locsS[i]]+=_dir*(
			_tableS[ii][0]*_sourceV[_ntBig*i+id]+
			_tableS[ii][1]*_sourceV[_ntBig*i+id+1]+
			_tableS[ii][2]*_sourceV[_ntBig*i+id+2]+
			_tableS[ii][3]*_sourceV[_ntBig*i+id+3]+
			_tableS[ii][4]*_sourceV[_ntBig*i+id+4]+
			_tableS[ii][5]*_sourceV[_ntBig*i+id+5]+
			_tableS[ii][6]*_sourceV[_ntBig*i+id+6]+
			_tableS[ii][7]*_sourceV[_ntBig*i+id+7]);
	}
}
void cpuProp::injectReceivers(int id, int ii, float *p){


	for(int i=0; i < _nRecs; i++) {
		p[_locsR[i]]+=_dir*(
			_tableS[ii][0]*_rec[_ntBig*i+id]+
			_tableS[ii][1]*_rec[_ntBig*i+id+1]+
			_tableS[ii][2]*_rec[_ntBig*i+id+2]+
			_tableS[ii][3]*_rec[_ntBig*i+id+3]+
			_tableS[ii][4]*_rec[_ntBig*i+id+4]+
			_tableS[ii][5]*_rec[_ntBig*i+id+5]+
			_tableS[ii][6]*_rec[_ntBig*i+id+6]+
			_tableS[ii][7]*_rec[_ntBig*i+id+7]);
	}
}
void cpuProp::prop(float *p0, float *p1, float *vel){
	for(int i3=4; i3 < _nz-4; i3++) {
		for(int i2=4; i2 < _ny-4; i2++) {
			int ii=i2*_nx+4+_n12*i3;
			for(int i1=4; i1 < _nx-4; i1++,ii++) {
				p0[ii]=vel[ii]*
					      (
					coeffs[C0]*p1[ii]
					+coeffs[CX1]*(p1[ii-1]+p1[ii+1])+
					+coeffs[CX2]*(p1[ii-2]+p1[ii+2])+
					+coeffs[CX3]*(p1[ii-3]+p1[ii+3])+
					+coeffs[CX4]*(p1[ii-4]+p1[ii+4])+
					+coeffs[CY1]*(p1[ii-_nx]+p1[ii+_nx])+
					+coeffs[CY2]*(p1[ii-2*_nx]+p1[ii+2*_nx])+
					+coeffs[CY3]*(p1[ii-3*_nx]*p1[ii+3*_nx])+
					+coeffs[CY4]*(p1[ii-4*_nx]+p1[ii+4*_nx])+
					+coeffs[CZ1]*(p1[ii-1*_n12]+p1[ii+1*_n12])+
					+coeffs[CZ2]*(p1[ii-2*_n12]+p1[ii+2*_n12])+
					+coeffs[CZ3]*(p1[ii-3*_n12]+p1[ii+3*_n12])+
					+coeffs[CZ4]*(p1[ii-4*_n12]+p1[ii+4*_n12])
				        )
				        +p1[ii]+p1[ii]-p0[ii];
			}
		}
	}
}
void cpuProp::transferSincTableD(int nsinc, int jtd, float **table){
// transfer_sinc_table_d(nsinc,jtd,myr.table);
	_nsincD=nsinc;
	_jtdD=jtd;
	_tableD=table;
}
void cpuProp::transferSourceFunc(int npts,int nt_big,int *locs, float *vals){
	_nptsS=npts; _ntBig=nt_big; _locsS=locs; _sourceV=vals;
}
void cpuProp::transferVelFunc1(int nx, int ny, int nz, float *vloc){
	_nx=nx; _ny=ny; _nz=nz; _vel1=vloc;
}
void cpuProp::transferVelFunc2(int nx, int ny, int nz, float *vloc){
	_nx=nx; _ny=ny; _nz=nz; _vel2=vloc;

}

void cpuProp::transferReceiverFunc(int nx, int ny, int nt, int *locs,
	float *rec){
	_nRecs=nx*ny; _nt=nt; _locsR=locs; _rec=rec;
}
void cpuProp::transferSincTableS(int nsinc, int jts, float **table){
	_nsincS=nsinc;
	_jtsS=jts;
	_tableS=table;
}
#define C_C00(d) (8.0/(5.0*(d)*(d)))

void cpuProp::createSpace(float d1, float d2, float d3,float bc_a, float bc_b, float bc_y,
	int nx, int ny, int nz){
	_bcA=bc_a; _bcB=bc_b; bc_y=_bcY;
	_nx=nx; _ny=ny; _nz=nz;
	_n12=_nx*_ny;
	coeffs.resize(13);
	coeffs[0]=-1025.0/576.0*(C_C00(d1)+C_C00(d2)+C_C00(d3));
	coeffs[1]=C_C00(d1);
	coeffs[2]=C_C00(d2);
	coeffs[3]=C_C00(d3);
	coeffs[4]=-C_C00(d1)/8.0;
	coeffs[5]=-C_C00(d2)/8.0;
	coeffs[6]=-C_C00(d3)/8.0;
	coeffs[7]=C_C00(d1)/63.0;
	coeffs[8]=C_C00(d2)/63.0;
	coeffs[9]=C_C00(d3)/63.0;
	coeffs[10]=-C_C00(d1)/896.0;
	coeffs[11]=-C_C00(d2)/896.0;
	coeffs[12]=-C_C00(d3)/896.0;
	_n123=_n12*nz;
	_bcB=.0005;
	_bcA=40;
	_bound.resize(40);
	for(int i=0; i < 40; i++) {
		_bound[i]=expf(-_bcB*(_bcA-i));
	}


	//create_gpu_space(d1,d2,d3,bc_a,bc_b,bc_b_y,nx,ny,nz);
}
