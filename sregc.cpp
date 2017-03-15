#include<sregc.h>
using std::string;
using std::vector;

void sregc::read_all(string tag,hypercube_complex *cube){

   if(this->get_n123() !=cube->get_n123()) 
   seperr("Dataset(%d) and sepcube(%d) not the same size\n", this->get_n123(),cube->get_n123());
   if(8*this->get_n123()!=sreed(tag.c_str(),cube->vals,this->get_n123()*8)) 
     seperr("Trouble reading data from tag %s\n",tag.c_str());
}
void sregc::read_next(string tag, hypercube_complex *cube){
  if(this->get_n123() < cube->get_n123()) seperr("Datasize smaller than hypercube\n");
  if(8*cube->get_n123()!=sreed(tag.c_str(),cube->vals,cube->get_n123()*8)) 
    seperr("Trouble reading data from tag %s\n",tag.c_str());
}
void sregc::write_all(string tag,hypercube_complex *cube){
   if(this->get_n123() !=cube->get_n123()) 
     seperr("Dataset (%d) and sepcube (%d) not the same size\n",
     this->get_n123(),cube->get_n123());
   if(8*cube->get_n123()!=srite(tag.c_str(),cube->vals,cube->get_n123()*8)) 
     seperr("Trouble reading data from tag %s\n",tag.c_str());
}
void sregc::write_next(string tag, hypercube_complex *cube){
  if(8*cube->get_n123()!=srite(tag.c_str(),cube->vals,cube->get_n123()*8)) 
    seperr("Trouble reading data from tag %s\n",tag.c_str());
}
void sregc::read_slice(string tag, vector<int> nwind,hypercube_complex *cube){
  vector<int> fwind(nwind),jwind(nwind);
  fwind.assign(jwind.size(),0);
  jwind.assign(jwind.size(),1);
  this->read_subset(tag,nwind,fwind,jwind,cube);
}
void sregc::read_slice(string tag, vector<int> nwind, vector<int> fwind, hypercube_complex *cube){
  vector<int> jwind(nwind);
  jwind.assign(jwind.size(),1);
  this->read_subset(tag,nwind,fwind,jwind,cube);
}
void sregc::read_slice(string tag, vector<int> nwind, vector<int>fwind, vector<int> jwind, hypercube_complex *cube){

  this->read_subset(tag,nwind,fwind,jwind,cube);
}
void sregc::write_slice(string tag, vector<int> nwind, hypercube_complex *cube){
  vector<int> fwind(nwind),jwind(nwind);
  fwind.assign(jwind.size(),0);
  jwind.assign(jwind.size(),1);
  this->write_subset(tag,nwind,fwind,jwind,cube);
}
void sregc::write_slice(string tag, vector<int> nwind, vector<int> fwind, hypercube_complex *cube){
  vector<int> jwind(nwind);
  jwind.assign(jwind.size(),1);
  this->write_subset(tag,nwind,fwind,jwind,cube);
}
void sregc::write_slice(string tag, vector<int> nwind, vector<int> fwind, vector<int> jwind, hypercube_complex *cube){

  this->write_subset(tag,nwind,fwind,jwind,cube);
}
void sregc::read_subset(string tag, vector<int> nwind, vector<int> fwind, vector<int> jwind, hypercube_complex *cube){
  if(nwind.size() !=jwind.size()
  || nwind.size()!=jwind.size() 
  || nwind.size()!=this->get_ndim())
    seperr("Size of windowing params and data dimensions not all equiv\n");
  int idim,size=8,ndims=this->get_ndim();
  int ng[ndims],nw[ndims],fw[ndims],jw[ndims];


  axis ax;
  for(int i=0 ;i < ndims;i++){ 
    ax=this->get_axis(i+1); 
    ng[i]=ax.n;
    nw[i]=nwind[i]; jw[i]=jwind[i]; fw[i]=fwind[i];
  }
  const char *c_str = tag.c_str ( );
 // fprintf(stderr,"VV %d %d %d - %d %d %d -%d %d %d    \n",
  // fw[0],fw[2],fw[8],nw[0],nw[2],nw[8],ng[0],ng[2],ng[8]);
  int ierr=sreed_window(c_str,&ndims,ng,nw,fw,jw,size,cube->vals);
 
  if(ierr!=0){
    for(int i=0; i < nwind.size(); i++){
      fprintf(stderr,"Window parameters (tag=%s) iax=%d ng=%d nw=%d fw=%d jw=%d \n",
        tag.c_str(),i+1,this->get_axis(i+1).n,nwind[i],fwind[i],jwind[i]);
    }
    seperr("trouble reading in window\n");
   }
  
}
void sregc::write_subset(string tag, vector<int> nwind, vector<int> fwind, vector<int> jwind, hypercube_complex *cube){
  if(nwind.size() !=jwind.size() || nwind.size()!=jwind.size() || nwind.size()!=this->get_ndim())
    seperr("Size of windowing params and data dimensions not all equiv\n");
  int idim,size=8,ndims=this->get_ndim();
  int ng[ndims],nw[ndims],fw[ndims],jw[ndims];
  for(int i=0 ; i < ndims;i++){ 
    ng[i]=this->get_axis(i+1).n; nw[i]=nwind[i]; jw[i]=jwind[i]; fw[i]=fwind[i];
  }
  int ierr=srite_window(tag.c_str(),&ndims,ng,nw,fw,jw,size,cube->vals);
  
  if(ierr!=0){
    for(int i=0; i < nwind.size(); i++){
      fprintf(stderr,"Window parameters (tag=%s) iax=%d ng=%d nw=%d fw=%d jw=%d \n",
        tag.c_str(),i+1,this->get_axis(i+1).n,nwind[i],fwind[i],jwind[i]);
    }
    seperr("trouble writing out window\n");
  }
}
  

