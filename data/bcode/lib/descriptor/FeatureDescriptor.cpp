#include "feature.h"
#include "../matrix/matrix.h"
#define ORI_THRESHOLD 0.90

 
void Segment::copy(Segment * seg){
  x.insert(x.begin(),seg->x.begin(),seg->x.end());
  y.insert(y.begin(),seg->y.begin(),seg->y.end());
  angle.insert(angle.begin(),seg->angle.begin(),seg->angle.end());
  grad.insert(grad.begin(),seg->grad.begin(),seg->grad.end());  
  ang_mean=seg->ang_mean; 
  x_mean=seg->x_mean; 
  y_mean=seg->y_mean; 
  //cout <<"xsize "<< x.size()<< " " << seg->x.size()<< endl;
}

char * getName(unsigned long  type){ 
  char *name = new char[512];
 // cout <<  " tp:"<< type << " " << (type&DETECTOR)<< endl;
  sprintf(name,"FEATURES");
  if((type&DHESSIAN)==DHESSIAN){
    sprintf(name,"%s-DHESSIAN",name);
  }
  if((type&DHARHES)==DHARHES){
    sprintf(name,"%s-DHARHES",name);
  }
  if((type&DHARRIS)==DHARRIS){
    sprintf(name,"%s-DHARRIS",name);
  }
  if((type&DSEDGE)==DSEDGE){
    sprintf(name,"%s-DSEDGE",name);
  }
  if((type&DMSER)==DMSER){
    sprintf(name,"%s-DMSER",name);
  }
  if((type&DJLA)==DJLA){
    sprintf(name,"%s-DJLA",name);
  }
  if((type&DSHAPE)==DSHAPE){
    sprintf(name,"%s-DSHAPE",name);
  }
  if((type&DCC)==DCC){
    sprintf(name,"%s-DCC",name);
  }
  if((type&DCLBP)==DCLBP){
    sprintf(name,"%s-CLBP",name);
  }
  if((type&DMOM)==DMOM){
    cout << (type&DMOM) << " " <<DMOM<< endl;
    sprintf(name,"%s-DMOM",name);
  }
  if((type&DKOEN)==DKOEN){
    sprintf(name,"%s-DKOEN",name);
  }
  if((type&DCLBP)==DCLBP){
    sprintf(name,"%s-DCLBP",name);
  }
  if((type&DPCA)==DPCA){
    sprintf(name,"%s-DPCA",name);
  }
  if((type&DSIFT)==DSIFT){
    sprintf(name,"%s-DSIFT",name);
  }
  if((type&DCSIFT)==DCSIFT){
    sprintf(name,"%s-DCSIFT",name);
  }
  if((type&DCOLOR)==DCOLOR){
    sprintf(name,"%s-DCOLOR",name);
  }
  if((type&DGLOH)==DGLOH){
    sprintf(name,"%s-DGLOH",name);
  }

  return name;
}


void computeHistAngle(DARY *grad_im, DARY *angle_im,vector<float> &angles);

FeatureDescriptor::~FeatureDescriptor(){
  if(par[DSIZE]>0)delete[] vec;
  delete[] par;
  if(cluster_dist!=NULL)delete cluster_dist;
  if(imagename!=NULL)delete[] imagename; 
  for(uint i=0;i<occurences.size();i++)delete occurences[i];
  occurences.clear();
  occ.clear();
  neg_match.clear();
  pos_weight.clear();
}
 


FeatureDescriptor::FeatureDescriptor(float xin, float yin, float scale_in, float featureness_in){
  init();
  //x=xin;y=yin;c_scale=scale_in;featureness=featureness_in;
  par[XPOS]=xin;par[YPOS]=yin;par[C_SCALE]=scale_in;par[FEATURNESS]=featureness_in;par[ANGLE]=1000;
 
}

void deleteDescriptors(vector<FeatureDescriptor *> &desc){
  for(uint i=0;i<desc.size();i++){
    //cout << i << " " << desc.size()<< endl;
    if(desc[i]!=NULL){
      if(desc[i]->features.size()>0)deleteDescriptors(desc[i]->features);
      delete desc[i];
    }
  }
  desc.clear(); 
}
void FeatureDescriptor::Cout(int dim){
  cout << "psize " << par[PSIZE] << " id " << par[ID] << "  x " << par[XPOS] << " y " << par[YPOS] << " f " << par[FEATURNESS] << " s " << par[C_SCALE] << " r " << par[RADIUS]<< " intL " << par[INT_LEV]<< " nbf " <<  par[NBF]<< " var " << par[VAR]<< " angle " << par[ANGLE]<< " extr " << par[EXTR]<< " type " << par[TYPE]<< " m11 " << par[MI11]<< " m21 " << par[MI21]<< " m22 " << par[MI22]<< " type " << par[TYPE];
  if(par[DSIZE]>0){
    cout<< " vec "; 
    for(int i=0;i<dim;i++)cout<< vec[i]<< " ";
  }
  cout << endl;
}

void  FeatureDescriptor::findNearestNeighbor(vector<FD*> features, int dim, int &ind, float &sim){
  sim=MAX_FLOAT;
  float d;
  ind=-1;
  float var=par[VAR];
  for(uint i=0;i<features.size();i++){
    if(features[i]!=NULL){
      d=var+features[i]->par[VAR]+distEuc(vec,features[i]->getVec(),dim);
      if(d<sim){
	sim=d;
	ind=i;
      }      
    }
  }  
}

 
/**********************************************/
void FeatureDescriptor::copy(FeatureDescriptor* ds){

  memcpy(par, ds->getPar(), sizeof(float)*((int)(ds->par[PSIZE])));
  imagename=NULL;
  cluster_dist=NULL; 
  if(par[DSIZE]>0){
    allocVec((int)par[DSIZE]);
    float *vec2=ds->getVec();
    memcpy(vec, vec2, sizeof(float)*(int)par[DSIZE]);
  }
   //Cout(10);
 //for(int i=0;i<size;i++)vec[i]=ds->getV(i);
}

void FeatureDescriptor::init(){
  init(PARAMS);
}
 
void FeatureDescriptor::init(int psize){
  par = new float[psize+1];
  bzero(par,(psize+1)*sizeof(float));
  par[PSIZE]=psize+1;
  if(psize==PARAMS){
    par[FEATURNESS]=1;
    par[ANGLE]=10000;
    par[C_SCALE]=1;
    par[L1]=1;
    par[MI11]=1;
    par[MI22]=1;
    par[NBF]=1;
    par[INT_LEV]=-1;
  }
  vec=NULL;
  imagename=NULL;
  cluster_dist=NULL;
}

/**********************************************/
void FeatureDescriptor::allocVec(int size_in){
    if(size_in >0){
	if(vec!=NULL)delete [] vec;
	par[DSIZE]=size_in;
	vec = new float[size_in];
	bzero(vec,size_in*sizeof(float));
    }
}
  
/**********************************************/
void  FeatureDescriptor::setImageName(const char* name) { 
  if(imagename==NULL)imagename=new char[512];
  strcpy(imagename,name);
}




/**********************************************/
void FeatureDescriptor::changeBase(float *mat){
  for(int v=0;v<(int)par[DSIZE];v++){
    vec[v] = vec[v]*mat[v];   
  }
}

void FeatureDescriptor::normalizeVect()
{
   int i;
   float val, fac, sqlen = 0.0;
  int size=(int)par[DSIZE];

   for (i = 0; i < size; i++) {
     val = vec[i];
     sqlen += val * val;
   }
   if(sqlen<0.0001){
     sqlen+=1;
   }
   fac = 200.0 / sqrt(sqlen);
   float uval;
   for (i = 0; i < size; i++){
      vec[i]= vec[i]*fac;
//       uval = vec[i]+127;
//       uval = (255 < uval) ? 255 : uval;
//       uval = (0 > uval) ? 0 : uval;
//       vec[i]=(int)(0.5+uval);

   }
   
}

/**********************************************/
void FeatureDescriptor::pca(int dim, float *avg, float *base){
  float *outvec = new float[dim];
  bzero(outvec,dim*sizeof(float));
  //for(int v=0;v<dim;v++)outvec[v]=0;
  int size=(int)par[DSIZE];
  for(int v=0;v<size;v++)vec[v]-=avg[v];  
  uint cnt=0;
  for(int i=0;i<dim;i++){
    for(int v=0;v<size;v++,cnt++){
      outvec[i] += vec[v]*base[cnt];   
    }
    //cout << outvec[i]<< endl;
  }
  //for(int i=0;i<dim;i++)vec[i]=outvec[i];
   
  memcpy(vec,outvec,dim*sizeof(float));
  //getchar();
  delete []outvec;
  par[DSIZE]=dim;
  //vec = outvec;
}



/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
void FeatureDescriptor::readCommon( ifstream &input, int size_in){
  double a,b,c;
  Matrix U(2,2,0.0),D,Vi,V;
  input >> par[XPOS];
  input >> par[YPOS];
  input >> a;
  input >> b; 
  input >> c;
  //cout << x << " " << y << " " << a << " ac " << c <<" " << b <<  endl;
  U(1,1)=a;
  U(1,2)=b;
  U(2,1)=b;
  U(2,2)=c;
    //cout << x << " " << y << " " << a << " ac " << c << endl;
    U.svd(Vi,D,V);
    D(1,1)=(1.0/sqrt(D(1,1)));
    D(2,2)=(1.0/sqrt(D(2,2)));
    a=sqrt(D(2,2)*D(1,1));
    //cout << D(1,1)<< " " <<  D(2,2)<< "  " << tmp1<< endl;
    D.tabMat[2][2]/=a;
    D.tabMat[1][1]/=a;
    U=V*D*V.transpose();
    par[C_SCALE]=a;
    //cout << x << " " << y << " " << a << " "<< b << " "<< c << " " << size_in<<endl;getchar();
    par[MI11]=U(1,1);
    par[MI12]=U(1,2);
    par[MI21]=U(2,1);
    par[MI22]=U(2,2);
  
  if(size_in>0){
    allocVec(size_in);
    for(int j=0;j<size_in;j++){
      input >> vec[j];
    }
  } 
} 
/**********************************************/
void FeatureDescriptor::writeCommon(ofstream &output){

    Matrix U(2,2,0.0),D,Vi,V;
    U(1,1)=par[MI11];
    U(1,2)=par[MI12];
    U(2,1)=par[MI21];
    U(2,2)=par[MI22];
    U.svd(Vi,D,V);
    D=D*(par[C_SCALE]);
    D(1,1)=1.0/(D(1,1)*D(1,1));
    D(2,2)=1.0/(D(2,2)*D(2,2));
    U=V*D*V.transpose();
    output << par[XPOS] << " " << par[YPOS]  << " " << U(1,1)<< " "<< U(1,2)<< " "<< U(2,2);
    if(par[DSIZE] >0){
      for(int j=0; j<(int)par[DSIZE] ;j++){
	output  << " " << vec[j];
      }
    }
    output << endl;
} 

/**********************************************/
void FeatureDescriptor::read(ifstream &input, int psize_in, int dsize_in){
  
  for(int i=0;i<psize_in;i++){
    input >> par[i];
  }
  if(psize_in==5){
    par[C_SCALE]=par[4];
    par[FEATURNESS]=1.0;
    
  }
  par[PSIZE]=(psize_in<PARAMS)?PARAMS:psize_in;//some internal parameters at positions beyond psize wouldn't be copied if <PARAMS
 
  if(dsize_in>0){
    allocVec(dsize_in);
    for(int j=0;j<(int)par[DSIZE];j++){
      input >> vec[j];
    }
  } 
} 

/**********************************************/
void FeatureDescriptor::write(ofstream &output){
  for(int i=0;i<(int)par[PSIZE];i++){
    output << par[i] << " ";
  }
  if(par[DSIZE]>0){
    for(int j=0; j<(int)par[DSIZE];j++){
      output << vec[j] << " ";
    }
  }
  output << endl;
}
/**********************************************/
void FeatureDescriptor::writeBin(FILE *fid, float *buf){
  fwrite(par, sizeof(float),(int)par[PSIZE], fid);
  fwrite(vec, sizeof(float),(int)par[DSIZE], fid);
}

 /**********************************************/
void FeatureDescriptor::readBin(FILE *fid, int size, float *buf){
  allocVec(size);
  fread(par, sizeof(float), (int)par[PSIZE], fid);
  fread(vec, sizeof(float), (int)par[DSIZE] , fid);
}



 
/**********************************************/
void loadTextFeatures( const char* points1, vector<FeatureDescriptor*> &cor1, int format){
    //cout << "Loading features from " << points1 << " with format " << format << endl;
    FeatureDescriptor* cor;
    ifstream input1(points1);
    if(!input1)return;
    int cor_nb1,dsize,psize; 
     if(format==FTEXT){
      char *bigbuf=new char[512];
      input1.getline(bigbuf,512);
      delete []bigbuf;
      input1 >> cor_nb1;
      input1 >> psize;
      input1 >> dsize;      
      
    }else if(format==OXFORD_FTEXT){
      psize=PARAMS;
      float tmp;
      input1 >> tmp;
      input1 >> cor_nb1;
      if(tmp<=1.0)dsize=0;
      else dsize=(int)tmp;
    } else if(format==AUDIO_FTEXT){
      cor_nb1=1000000;
      dsize=39;
      psize=0;
    }
 
    if(cor_nb1==0)return;
    uint cn=0,flag,tm;
    while(input1.good() && cn<cor_nb1){
      cor = new FeatureDescriptor(psize);
      flag=1;
      if(format==OXFORD_FTEXT)
        cor->readCommon(input1,dsize);
      else if(format==FTEXT) {
        cor->read(input1,psize,dsize);
      }else if(format==AUDIO_FTEXT){
	cor->read(input1,psize,dsize);
	tm=(uint)(cn/100);//time/100nanosec=1milisec
	cor->setX_Y(tm,0);
	cor->setScale(1.0);
	cor->setType(DAUDIO);	
	if(fabs(cor->vec[38])>100.0 || fabs(cor->vec[0])==0.0  || fabs(cor->vec[12])>100.0 || fabs(cor->vec[25])>100.0)flag=0;
      }
      cn++;
      //cor->Cout(39);
      if(flag)cor1.push_back(cor); //cout << "read " << psize << "  " << cor1.size()<<endl;cor->Cout();
      else delete cor;
      //cor->Cout(0);cout << dsize << " DSIZE " << cor->par[DSIZE]<< endl;getchar();
    }
    if(cor_nb1!=(int)cor1.size()){
      cout << "warning:"<< endl<<"in header: "<<cor_nb1 << ", in file: " << points1 << " "<< cor1.size()<< endl; 
    }
   
    input1.close();
}


/**********************************************/

void loadBinFeatures(const char* points_out, vector<FeatureDescriptor*> &cor){
  //cout << "Loading features from " << points_out << "... "<<  flush;
    
  FILE *fid;
  fid=fopen(points_out, "rb");

  if(!fid){
	cout << "loadBinFeatures: error opening " << points_out<< endl;
 	return;
  }
   
  uchar format;  
  fread(&format, sizeof(uchar), 1, fid);
  uchar psi;  
  fread(&psi, sizeof(uchar), 1, fid);
  uint psize=(uint)psi;
  uint size;  
  fread(&size, sizeof(uint), 1, fid);
  uint nb;
  fread(&nb, sizeof(uint), 1, fid);
  //cout << "psize "<< psize << " size " << size << " nb " << nb<< endl;
  if(nb==0){ 
    cout << "No features in this file "<< points_out << endl;
    fclose(fid); 
    return;
  }

  
  if(format==FBIN){
  uint totsize=(size+psize)*nb;
  float *write_vec = new float[totsize];
  fread(write_vec, sizeof(float), totsize , fid);

  int pos=0;
  uint fstart=cor.size();
  nb+=fstart;
  for(uint i=fstart; i<nb;i++){
    cor.push_back(new FeatureDescriptor());
    cor[i]->allocVec(size);
    for(uint p=0;p<psize;p++)cor[i]->par[p]=write_vec[pos++];
    for(uint p=0;p<size;p++)cor[i]->vec[p]=write_vec[pos++];
    //cor[i]->Cout();getchar();
  }  
  delete []write_vec;
    
  }else if(format==COMPACT_FBIN){
   uint vecsize=(psize+size)*nb;
   uchar *write_vec_byte = new uchar[vecsize];
   fread(write_vec_byte, sizeof(uchar), vecsize , fid);


  uint pos=0;
  uint cn=0,tmp1,tmp2;
  uint fstart=cor.size();
  nb+=fstart;
  uchar byte1,byte2,byte3,byte4;
  for(uint i=fstart; i<nb;i++){
   cor.push_back(new FeatureDescriptor());
   cor[i]->allocVec(size);   
   byte1=write_vec_byte[pos++];byte2=write_vec_byte[pos++];
   cor[i]->setX_Y((float)byte1,(float)byte2);
   byte1=write_vec_byte[pos++];
   cor[i]->setScale((float)byte1);
   byte1=write_vec_byte[pos++];  
   cor[i]->setAngle((float)byte1);
   
    for(int p=0;p<size;p++)cor[i]->vec[p]=(float)write_vec_byte[pos++];
    //cor[i]->Cout(10);getchar();
  }  
  delete []write_vec_byte;
    
  }else {
        cout << "loadBinFeatures: wrong format read "<< (int)format << " in "<< points_out << endl;
        fclose(fid); 
	exit(0);
    
  }
  fclose(fid); 
} 


/**********************************************/
void loadAndProjectFeatures( const char* points1, vector<FeatureDescriptor*> &cor1, int format, float *spca_mean, float *spca_base, uint desc_dim){
	loadFeatures( points1, cor1,  format);
	pca(cor1, spca_mean, spca_base, desc_dim);
	
}

void loadAndProjectFeatures(vector<char *> filenames, vector<FD *> &features, ulong max_number_features, ulong max_memory_features, int format, float *spca_mean, float *spca_base, uint desc_dim) {
	
    char filename[512];
    int nbp;
    ulong current_mem=0;
    uint file_nb=0,fsize;
    while (features.size() < max_number_features && file_nb< filenames.size() && current_mem <  max_memory_features) {
        //cout << "loadFeatures " << filenames.size() << " " << file_nb<< "  " <<filenames[file_nb] << endl;
		fsize=features.size();
		loadAndProjectFeatures( filenames[file_nb], features, format, spca_mean, spca_base, desc_dim);
		for(uint i=fsize;i<features.size();i++)features[i]->par[IM_ID]=file_nb;
		current_mem=4*(features[0]->par[PSIZE]+features[0]->par[DSIZE])*features.size();
        if(!(file_nb%100))cout << "loadFeatures file " << file_nb  << " of  " << filenames.size()  << " memory: " << current_mem/1048576 << "kB <? "<< max_memory_features/1048576 << "kB  features: " << features.size() << " <? " << max_number_features << endl;
        file_nb++;
    }
	
}


/**********************************************/
void loadFeatures( const char* points1, vector<FeatureDescriptor*> &cor1, int format){
  if(format==FBIN || format==COMPACT_FBIN)loadBinFeatures(points1,cor1);
  else loadTextFeatures( points1, cor1,  format);

}

void loadFeatures(vector<char *> filenames, vector<FD *> &features, ulong max_number_features, ulong max_memory_features, int format) {

    char filename[512];
    int nbp;
    ulong current_mem=0;
    uint file_nb=0,fsize;
    while (features.size() < max_number_features && file_nb< filenames.size() && current_mem <  max_memory_features) {
        //cout << "loadFeatures " << filenames.size() << " " << file_nb<< "  " <<filenames[file_nb] << endl;
	fsize=features.size();
	loadFeatures( filenames[file_nb], features, format);
	for(uint i=fsize;i<features.size();i++)features[i]->par[IM_ID]=file_nb;
        current_mem=4*(features[0]->par[PSIZE]+features[0]->par[DSIZE])*features.size();
        if(!(file_nb%100))cout << "loadFeatures file " << file_nb  << " of  " << filenames.size()  << " memory: " << current_mem/1048576 << "kB <? "<< max_memory_features/1048576 << "kB  features: " << features.size() << " <? " << max_number_features << endl;
        file_nb++;
    }

}

/**********************************************/
void writeTextFeatures(vector<FeatureDescriptor*> cor, const char* points_out, int format){
    ofstream output(points_out);   
    if(!output)cout << "error opening " << points_out<< endl;
    int size=0,psize=0;;
    
    if (cor.size()==0){
      cout << " descriptors nb " << cor.size() << endl;
     
    }else {
      psize=cor[0]->par[PSIZE];
      size=cor[0]->par[DSIZE];
    }
    
    
    
    if(format==FTEXT){
      output << "#comments psize dsize x y cornerness scale angle ...." << endl; 
      output << cor.size()<< endl;
      output << psize << endl;
      output << size  << endl;
    }else{
      output << size<< endl;
      output << cor.size()<< endl;
    }
    for(uint i=0; i<cor.size();i++){
      if(format==OXFORD_FTEXT)
        cor[i]->writeCommon(output);
      else  cor[i]->write(output);
    }
    output.close();  
} 

/**********************************************/

void loadBOcc(const char* points_out, vector<BOccList> &wbocc){
  cout << "Loading occlists from " << points_out << "... "<<  flush;
    
  FILE *fid;
  fid=fopen(points_out, "rb");

  if(!fid){
	cout << "loadBOcc: error opening " << points_out<< endl;
 	return;
  }
   
  uint nb;
  fread(&nb, sizeof(uint), 1, fid);
  cout << "nb "<< nb << endl;
  if(nb==0 || wbocc.size()>0){ 
    cout << "No occ in this file or wbocc is not empty"<< points_out << endl;
    fclose(fid); 
    return;
  }
  uint *write_vec_uint = new uint[nb];
  fread(write_vec_uint, sizeof(uint), nb , fid);
  uint occnb=0;
  for(uint i=0; i<nb;i++){
    wbocc.push_back(BOccList(0));
    occnb+=write_vec_uint[i];
    //cout << i << " " <<  write_vec_uint[i] << endl;
  }
  //cout  << "occnb " << occnb << endl;getchar();
  uchar psize=17;//9 bytes for x,y,s,i:4B,a,w,df1:4B,df2:4B
   //cout << "psize "<< psize << " size " << size << " nb " << nb<< endl;
  uint vecsize=(psize)*occnb;
  uchar *write_vec_byte = new uchar[vecsize];
  fread(write_vec_byte, sizeof(uchar), vecsize , fid);

  uint pos=0;
  uint cn=0,tmp;
  uint byte1,byte2,byte3,byte4;
  BOcc *occ;
  for(uint i=0; i<nb;i++){
    for(uint j=0; j<write_vec_uint[i];j++){
    occ = new BOcc();
    byte1=write_vec_byte[pos++];byte2=write_vec_byte[pos++];byte3=write_vec_byte[pos++];byte4=write_vec_byte[pos++];
    byte42int(byte1, byte2, byte3, byte4, tmp );
    occ->o=tmp;
    byte1=write_vec_byte[pos++];byte2=write_vec_byte[pos++];byte3=write_vec_byte[pos++];byte4=write_vec_byte[pos++];
    byte42int(byte1, byte2, byte3, byte4, tmp );
    occ->df1=tmp;
    byte1=write_vec_byte[pos++];byte2=write_vec_byte[pos++];byte3=write_vec_byte[pos++];byte4=write_vec_byte[pos++];
    byte42int(byte1, byte2, byte3, byte4, tmp );
    occ->df2=tmp;
    occ->x=write_vec_byte[pos++];
    occ->y=write_vec_byte[pos++];
    occ->s=write_vec_byte[pos++];
    occ->a=write_vec_byte[pos++];
    occ->w=write_vec_byte[pos++];
    wbocc[i].push_back(occ);
    //cout << i << " "<< j << " load " << (int)occ->o << " " << (int)occ->x << " " << (int)occ->y << " " << (int)occ->s << " " << (int)occ->a << " " << (int)occ->w << endl;
    }
  }  
 
  delete []write_vec_byte;
  delete []write_vec_uint;
  fclose(fid); 
  cout  << occnb << " occnb."  << endl;

} 
 

/**********************************************/
void saveBOcc(vector<BOccList> wbocc, const char* points_out){
  //cout << "Saving " << wbocc.size() << " occlists  in "  << points_out <<  "  "<< flush;
  if (wbocc.size()==0){
    cout << "saveBOcc: saving  0 number of lists" << endl;
    //return; 
  }else{ 
    //cout << cor[0]->getSize() << " dim in " << points_out << "... "<< flush;
  }	

  FILE *fid;
  fid=fopen(points_out, "wb");

  if(!fid){
	cout << "error opening " << points_out<< endl;
        exit(1);
  }	

  uchar psize=17;//9 bytes for x,y,s,i:4B,a,w
  //if(format==COMPACT_FBIN)psize=10;//14 bytes for
  
  uint nb=wbocc.size();
  fwrite(&nb, sizeof(uint), 1, fid);
  if(nb==0){
     cout << "done with 0 occ "<< endl; 
     fclose(fid); 
     return;
  }

  uint *write_vec_uint = new uint[nb];

  uint occnb=0;
  for(uint i=0; i<nb;i++){
    occnb+=wbocc[i].size();
    write_vec_uint[i]=wbocc[i].size(); 
  }
  //cout  << "occnb " << occnb << endl;

  fwrite(write_vec_uint, sizeof(uint), nb, fid);	
    if(occnb>0){
      uint vecsize=(psize)*occnb;
      uchar *write_vec_byte = new uchar[vecsize],b1,b2,b3,b4;
      uint pos=0,tmp;
      for(uint i=0; i<nb;i++){
        for(uint j=0; j<write_vec_uint[i];j++){
	tmp=(uint)wbocc[i][j]->o;
	int24byte((uint)wbocc[i][j]->o, b1, b2, b3, b4);
	write_vec_byte[pos++]=b1;write_vec_byte[pos++]=b2;write_vec_byte[pos++]=b3;write_vec_byte[pos++]=b4;
	int24byte((uint)wbocc[i][j]->df1, b1, b2, b3, b4);
	write_vec_byte[pos++]=b1;write_vec_byte[pos++]=b2;write_vec_byte[pos++]=b3;write_vec_byte[pos++]=b4;
	int24byte((uint)wbocc[i][j]->df2, b1, b2, b3, b4);
	write_vec_byte[pos++]=b1;write_vec_byte[pos++]=b2;write_vec_byte[pos++]=b3;write_vec_byte[pos++]=b4;
	write_vec_byte[pos++]=wbocc[i][j]->x;
	write_vec_byte[pos++]=wbocc[i][j]->y;
	write_vec_byte[pos++]=wbocc[i][j]->s;
	write_vec_byte[pos++]=wbocc[i][j]->a;
	write_vec_byte[pos++]=wbocc[i][j]->w;
	//cout <<"save " <<  wbocc[i][j]->o << " " << (int)wbocc[i][j]->x << " " << (int)wbocc[i][j]->y << " " << (int)wbocc[i][j]->s << " " << (int)wbocc[i][j]->a << " " << (int)wbocc[i][j]->w << endl;

	}  
      }
      fwrite(write_vec_byte, sizeof(uchar), vecsize , fid);
      delete [] write_vec_byte;
    }
  delete [] write_vec_uint;
  fclose(fid); 
  cout << occnb << " done"<< endl;
} 
 


/**********************************************/
void saveBinFeatures(vector<FeatureDescriptor*> cor, const char* points_out, int format){
  cout << "Saving " << cor.size() << " features of format "<< format <<" in "  << points_out <<  "  "<< flush;
  if (cor.size()==0){
    cout << "saveBinFeatures: saving  0 number of features" << endl;
    //return; 
  }else{ 
    //cout << cor[0]->getSize() << " dim in " << points_out << "... "<< flush;
  }	

  FILE *fid;
  fid=fopen(points_out, "wb");

  if(!fid){
	cout << "error opening " << points_out<< endl;
        exit(1);
  }	

  uint size=0;
  uchar psize=1;
  if(cor.size()){
    size=(uint)cor[0]->getSize();
    psize=(uchar)cor[0]->par[PSIZE];
  }
  uint nb=cor.size();
  if(format==COMPACT_FBIN){
    psize=4;//14 bytes for
  }
  
  uchar form = (uchar)format;
  fwrite(&form, sizeof(uchar), 1, fid);
  fwrite(&psize, sizeof(uchar), 1, fid);
  fwrite(&size, sizeof(uint), 1, fid);
  fwrite(&nb, sizeof(uint), 1, fid);
  if(nb==0){
     cout << "done with 0 features "<< endl; 
     fclose(fid); 
     return;
  }


 if(format==FBIN){
  int pos=0;
  uint totsize=(size+psize)*nb;
  float *write_vec = new float[totsize];
  for(uint i=0; i<cor.size();i++){    
    for(int p=0;p<psize;p++)write_vec[pos++]=cor[i]->par[p];
    for(int p=0;p<size;p++)write_vec[pos++]=cor[i]->vec[p];
  }  
  fwrite(write_vec, sizeof(float),totsize , fid);
  delete [] write_vec;
      
 }else if(format==COMPACT_FBIN){
   
  uint vecsize=(psize+size)*nb;
  uchar *write_vec_byte = new uchar[vecsize];
  float uval;
  uint pos=0,tmp;
  uchar byte1, byte2, byte3, byte4;
  for(uint i=0; i<cor.size();i++){
      write_vec_byte[pos++]=(uchar)cor[i]->getX();
      write_vec_byte[pos++]=(uchar)cor[i]->getY();
      write_vec_byte[pos++]=(uchar)cor[i]->getScale();
      write_vec_byte[pos++]=(uchar)cor[i]->getAngle(); 
      for(int p=0;p<size;p++){
       if(cor[i]->vec[p]<0 || cor[i]->vec[p]>255){
	cout << "ERROR saveBinFeatures: saving values outside of byte range cor "<< i << " val " << p << endl;
	exit(1);
       }
        write_vec_byte[pos++]=(uchar)cor[i]->vec[p];
       }
  }  
  fwrite(write_vec_byte, sizeof(uchar), vecsize , fid);
  delete [] write_vec_byte;
   
 }

  fclose(fid); 
  cout << "done"<< endl;
} 
 
/**********************************************/
void saveFeatures( vector<FeatureDescriptor*> cor1, const char* points_out, int format){
  if(format==FBIN || format==COMPACT_FBIN)saveBinFeatures(cor1,points_out,format);
  else writeTextFeatures(  cor1,  points_out, format);
}





/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/




/**********************************************/
void smoothHistogram(float *hist, int bins)
{
	int i;
	float prev, temp;
	float first=hist[0];
	prev = hist[bins - 1];
	for (i = 0; i < bins; i++) {
		temp = hist[i];
		hist[i] = (prev + hist[i] + ((i + 1 == bins) ? first : hist[i + 1])) / 3.0;
		prev = temp;
	}
}

void findDominantAngles(float *hist, int OriBins, vector<float> &angles, float range_angle, float max_angle){
  for (int i = 0; i < 6; i++)
    smoothHistogram(hist, OriBins);
	
	
  /* Find maximum value in histogram. */
  int maxi=-1; 
  float maxval=0;
  int maxi2=-1;
  int vec_pos=-1;
  float maxval2=0;
  for (int i = 0; i < OriBins; i++){          
    if (hist[i] > maxval){
      maxval = hist[i];
      maxi=i;     
    } 
  }
    
  int prev = (maxi == 0 ? OriBins - 1 : maxi - 1);
  int next = (maxi == OriBins - 1 ? 0 : maxi + 1);
  float interp = interpPeak(hist[prev], hist[maxi], hist[next]);
  float langle=range_angle * (maxi + interp) / OriBins;
  while(langle>max_angle)langle-=range_angle;
  angles.push_back(langle); 
        //cout << "an=[ " ;
  for (int i = 0; i < OriBins; i++){
    prev = (i == 0 ? OriBins - 1 : i - 1);
    next = (i == OriBins - 1 ? 0 : i + 1);
	  //cout  << " "<< hist[i] << " ";
    if (hist[i] > hist[prev] && hist[i] > hist[next] && hist[i]/maxval>ORI_THRESHOLD && i!=maxi){
      if(hist[i]/maxval>0.95 && hist[i]>maxval2){
	maxval2=hist[i];
	maxi2=i;
	vec_pos=angles.size();
      }
      interp = interpPeak(hist[prev], hist[i], hist[next]);
	    //angles.push_back(M_2PI * (i + 0.5 + interp) / OriBins - M_PI);
      langle=range_angle * (i + interp) / OriBins;
      while(langle>max_angle)langle-=range_angle;
      angles.push_back(langle); 
	    //cout <<i <<" " <<  360*angles[angles.size()-1]/M_2PI;      
    }	   //cout << endl;
  }
  
  if(maxi2>0){
    if(((maxi-maxi2)>0 && (maxi-maxi2)<(OriBins>>1)) || ((maxi2-maxi)>0 && (maxi2-maxi)>(OriBins>>1))){
      float tmp=angles[0];
      angles[0]=angles[vec_pos];
      angles[vec_pos]=tmp;
    }
  }
  
  return;

  if(angles.size()>1){
    //measure the minimum distance between angles
    float da = fabs(angles[0]-angles[angles.size()-1]);
    int mini=0;
    if(da>M_PI){
    da=M_2PI-da;
    mini=angles.size()-1;
    }else{
      mini=0;
    }
    float dat;
    for(uint i=1;i<angles.size();i++){
      dat=fabs(angles[i]-angles[i-1]);
      if(dat>M_PI){
        dat=M_2PI-dat;
      }
      if(dat<da){
        if(fabs(angles[i]-angles[i-1])>M_PI){
          mini=i;
        }else{
          mini=i-1;
        }
      }
      
    }
    langle=angles[mini];
    angles.clear();
    angles.push_back(langle); 
  }
}

/******************Lowe method*********************/
void computeHistAngle(DARY *grad, DARY *ori, int x, int y, int size,float &angle){
	if(ori==NULL ||grad==NULL ){return;}
	int i, r, c, rows, cols, bin, prev, next,OriBins=36;
	float hist[OriBins],  gval, langle,dbin,fbin,  interp, maxval = 0.0;/*OriBins=36*/
   
	rows = grad->y();
	cols = grad->x();
   
	for (i = 0; i < OriBins; i++)
		hist[i] = 0.0;

	int xs=((x-size)<0)?0:x-size;
	int xe=((x+size)<(int)grad->x())?x+size:(grad->x()-1);
	int ys=((y-size)<0)?0:y-size;
	int ye=((y+size)<(int)grad->y())?y+size:(grad->y()-1);
	float non_zero_grad=0;
	//cout <<x << " y " << y << " size " << size << " xs "<< xs <<" xe " << xe <<" ys "<< ys << " ye " << ye << endl;
	for (r = ys; r <= ye; r++)
	  for (c = xs; c <= xe; c++){
	    gval = grad->fel[r][c];   
	    if (gval > 15.0) {
	      /* Ori is in range of -PI to PI. */
	      langle = ori->fel[r][c] + M_PI;
	      fbin =  (OriBins * langle / (M_2PI));
	      bin=(int)floor(fbin);
	      dbin=fbin-bin;
	      //assert(bin >= 0 && bin <= OriBins);
	      bin = (bin < OriBins)?bin:(0);
	      hist[bin] +=  (1-dbin)*gval;
	      bin = (bin+1 < OriBins)?bin+1:(0);
	      hist[bin] +=  (dbin) * gval;
	      non_zero_grad++;
	    }
	  }
	/* Apply smoothing 6 times for accurate Gaussian approximation. */
	for (i = 0; i < 6; i++)
	  smoothHistogram(hist, OriBins);
	
	
	/* Find maximum value in histogram. */
	int maxi=-1; 
	for (i = 0; i < OriBins; i++){
	  //cout << hist[i] << " ";
	  if (hist[i] > maxval){
	    maxval = hist[i];
	    maxi=i;
	  }
	}
	//cout << endl;

	non_zero_grad=(non_zero_grad/((2*size+1)*(2*size+1)));
	if(non_zero_grad < 0.03 || maxi==-1){
	  //cout << (non_zero_grad/(4*size*size))<< endl;
	  angle=NO_MOTION;
	  return;
	}


	prev = (maxi == 0 ? OriBins - 1 : maxi - 1);
	next = (maxi == OriBins - 1 ? 0 : maxi + 1);
	interp = interpPeak(hist[prev], hist[maxi], hist[next]);
	angle=(M_2PI * (maxi + interp) / OriBins); 
	
	while(angle>M_PI)angle-=M_2PI;

	/*	
	cout << maxval<< " " << 180*angle/M_PI<< " maxi " << maxi<< " " << non_zero_grad << endl;

    DARY*  grad2 = new DARY(grad);
    grad2->fel[y-1][x]=-0.01;grad2->fel[y][x+1]=-0.01;
    grad2->fel[y][x-1]=-0.01;grad2->fel[y+1][x]=-0.01;
    grad2->writePNG("grad.png");
    delete grad2;	
    DARY*  ori2 = new DARY(ori);
    ori2->writePNG("ori.png");
    delete ori2;	
	getchar();
	*/
}
/******************Lowe method*********************/
void computeHistAngle(int x, int y, int radius, DARY *grad, DARY *ori,vector<float> &angles){
	if(ori==NULL ||grad==NULL ){return;}
	int i, r, c, rows, cols, bin, prev, next, OriBins=18;
	float hist[OriBins],  gval, langle, angle, dbin, fbin,  interp, maxval = 0.0;/*OriBins=36*/
	for (i = 0; i < OriBins; i++)
		hist[i] = 0.0;
	float gauss=0;
        //patch_mask->normalize(0,1);patch_mask->writePNG("patch_mask.png");
        //cout << radius << endl;
	for (int j = -radius; j <= radius; j++){
	  for (int i = -radius; i <= radius; i++){
	   if(x+i>=0 && y+j>=0 && x+i<grad->x() && y+j<grad->y()){
	    gval = grad->fel[y+j][x+i];
           {
	      gauss=0;
              r=(int)(0.5+GAUSS_SIZE*(sqrt(i*i+j*j)/((float)radius)));
              if(r<GAUSS_SIZE){
                gauss  =   gauss_100[r]; 
              }
            }
	    if (gval > 1.0  &&  gauss>0) {
	      /* Ori is in range of -PI to PI. */
	      langle = ori->fel[y+j][x+i] + M_PI;
	      fbin =  (OriBins * langle / (M_2PI));
	      bin=(int)floor(fbin);
	      dbin=fbin-bin;
	      //assert(bin >= 0 && bin <= OriBins);
	      bin = (bin < OriBins)?bin:(0);
	      hist[bin] +=  (1-dbin)*gval * gauss;
	      bin = ((bin+1) < OriBins)?(bin+1):(0);
	      hist[bin] +=  (dbin) * gval * gauss;
	    }
           } 
	  }
	}
        
 	//for (i = 0; i < 3; i++)
	//  smoothHistogram(hist, OriBins);
	
       findDominantAngles(hist,  OriBins, angles,M_2PI,M_PI);
        return;
        //cout << "OK " << endl;
	/* Apply smoothing 6 times for accurate Gaussian approximation. */
	for (i = 0; i < 6; i++)
	  smoothHistogram(hist, OriBins);
	
	
	/* Find maximum value in histogram. */
	int maxi=-1; 
	for (i = 0; i < OriBins; i++){          
	  if (hist[i] > maxval){
	    maxval = hist[i];
	    maxi=i;
	  }
	}
        prev = (maxi == 0 ? OriBins - 1 : maxi - 1);
        next = (maxi == OriBins - 1 ? 0 : maxi + 1);
        interp = interpPeak(hist[prev], hist[maxi], hist[next]);
        langle=M_2PI * (maxi + interp) / OriBins;
        while(langle>M_PI)langle-=M_2PI;
        angles.push_back(langle); 
        //cout << "an=[ " ;
        for (i = 0; i < OriBins; i++){
	  prev = (i == 0 ? OriBins - 1 : i - 1);
	  next = (i == OriBins - 1 ? 0 : i + 1);
	  //cout  << " "<< hist[i] << " ";
	  if (hist[i] > hist[prev] && hist[i] > hist[next] && hist[i]/maxval>ORI_THRESHOLD && i!=maxi){
	    interp = interpPeak(hist[prev], hist[i], hist[next]);
	    //angles.push_back(M_2PI * (i + 0.5 + interp) / OriBins - M_PI);
	    langle=M_2PI * (i + interp) / OriBins;
	    while(langle>M_PI)langle-=M_2PI;
	    angles.push_back(langle); 
	    //cout <<i <<" " <<  360*angles[angles.size()-1]/M_2PI;      
	  }
	   //cout << endl;
	}   
        //cout << "]"<<endl<< maxi << endl;
}


/******************Lowe method*********************/
void computeHistAngle(DARY *grad, DARY *ori,vector<float> &angles){
	if(ori==NULL ||grad==NULL ){return;}
	int i, r, c, rows, cols, bin, prev, next,OriBins=36;
	float hist[OriBins],  gval, langle, angle, dbin, fbin,  interp, maxval = 0.0;/*OriBins=36*/
   
	rows = grad->y();
	cols = grad->x();
   
	for (i = 0; i < OriBins; i++)
		hist[i] = 0.0;

	for (r = 1; r < rows - 1; r++)
	  for (c = 1; c < cols - 1; c++){
	    gval = grad->fel[r][c];   
	    if (gval > 1.0  &&  patch_mask->fel[r][c]>0) {
	      /* Ori is in range of -PI to PI. */
	      langle = ori->fel[r][c] + M_PI;
	      fbin =  (OriBins * langle / (M_2PI));
	      bin=(int)floor(fbin);
	      dbin=fbin-bin;
	      //assert(bin >= 0 && bin <= OriBins);
	      bin = (bin < OriBins)?bin:(0);
	      hist[bin] +=  (1-dbin)*gval * patch_mask->fel[r][c];
	      bin = (bin+1 < OriBins)?bin+1:(0);
	      hist[bin] +=  (dbin) * gval * patch_mask->fel[r][c];
	    } 
	  }
          findDominantAngles(hist,  OriBins, angles,M_2PI,M_PI);
          return;
	/* Apply smoothing 6 times for accurate Gaussian approximation. */
	for (i = 0; i < 6; i++)
	  smoothHistogram(hist, OriBins);
	
	
	/* Find maximum value in histogram. */
	int maxi=-1; 
	for (i = 0; i < OriBins; i++){
	  if (hist[i] > maxval){
	    maxval = hist[i];
	    maxi=i;
	  }
	}
	for (i = 0; i < OriBins; i++){
	  prev = (i == 0 ? OriBins - 1 : i - 1);
	  next = (i == OriBins - 1 ? 0 : i + 1);
	  //cout << i  << " "<< hist[i] << " ";
	  if (hist[i] > hist[prev] && hist[i] > hist[next] && hist[i]/maxval>ORI_THRESHOLD){
	    interp = interpPeak(hist[prev], hist[i], hist[next]);
	    //angles.push_back(M_2PI * (i + 0.5 + interp) / OriBins - M_PI);
	    langle=M_2PI * (i + interp) / OriBins;
	    while(langle>M_PI)langle-=M_2PI;
	    angles.push_back(langle); 
	    //cout <<i <<" " <<  360*angles[angles.size()-1]/M_2PI;      
	  }
	   //cout << endl;
	}   
	//cout << endl;getchar();
}

void normalize(DARY * img, int x, int y, float radius){
	float sum=0;
	float gsum=0; 

	for(uint j=0;j<img->y();j++){ 
		for(uint i=0;i<img->x();i++){ 
			if(patch_mask->fel[j][i]>0){
				sum+=img->fel[j][i]; 
				gsum++;
			}
		} 
	}    
	sum=sum/gsum;
	float var=0;
	for(uint j=0;j<img->y();j++){ 
		for(uint i=0;i<img->x();i++){ 
			if(patch_mask->fel[j][i]>0){	
				var+=(sum-img->fel[j][i])*(sum-img->fel[j][i]);	
			}
		}
	}     
	var=sqrt(var/gsum);    

    //  cout << "mean "<<sum<< " " <<img->fel[y][x] << " var " << var << endl;
	float fac=50.0/var;
	float max=0,min=1000;
	for(uint j=0;j<img->y();j++){ 
		for(uint i=0;i<img->x();i++){ 
			img->fel[j][i]=128+fac*(img->fel[j][i]-sum);
			if(max<img->fel[j][i])max=img->fel[j][i];
			if(min>img->fel[j][i])min=img->fel[j][i];
			if(img->fel[j][i]>255)img->fel[j][i]=255;
			if(img->fel[j][i]<0)img->fel[j][i]=0;
		}
	}   
    // cout << "max " << max << " min "<< min <<endl;
}



/************NORMALIZATION PATCH****************/
//int PATCH_SIZE=33;
DARY *patch_mask;
//float PATCH_SUM;
void initPatchMask(int size){ 
        patch_mask = new DARY(size,size,FLOAT1,1.0);
	int center=size>>1;
	float radius = center*center;
	float sigma=8.0*radius;
	float disq;
	for(int i=0;i<size;i++)
	  for(int j=0;j<size;j++){
	    disq=(i-center)*(i-center)+(j-center)*(j-center);
	    if(disq < radius){
	      patch_mask->fel[j][i]= exp(- disq / sigma);
	      //mask->fel[j][i]= 255*exp(- disq / sigma);   
	      //cout << patch_mask->fel[j][i]<< endl; 
	      //PATCH_SUM+=patch_mask->fel[j][i];
	    }else { 
	      patch_mask->fel[j][i]=0;
	    }		
	  } 
	
	//patch_mask->normalize(0,1);patch_mask->write("mask.pgm");cout << "mask "<< endl;getchar();
} 


void  compDesc(DARY *patch, FD *ds, uint type){
  
  if((type&DJLA)==DJLA){
    //computeJLA(patch, ds);
  }else if((type&DSHAPE)==DSHAPE){
    computeShape(patch, ds);
  }else if((type&DCC)==DCC){
    //computeCC(patch, ds);
  }else if((type&DMOM)==DMOM){
    //computeMoments(patch, ds);
  }else if((type&DKOEN)==DKOEN){
    //computeKoen(patch, ds);
  }else if((type&DSPIN)==DSPIN){
    //computeSpin(patch, ds);
  }else if((type&DPCA)==DPCA){
    //computePcaSift(patch, ds);
  }else if((type&DSIFT)==DSIFT){
    //computeSift(patch, ds);
  }else {
    //cout <<"compDesc: no such descriptor "<< type<<endl;// exit(0);
  }

}


void  compDesc(DARY *dx, DARY *dy, DARY *grad, DARY *ori, FD *ds, uint type){
  

  if((type&DSIFT)==DSIFT){
    computeSift(dx, dy ,grad, ori , ds);
  }else {
    //cout <<"compDesc dx: no such descriptor "<< type<< " " << DSIFT << endl; //exit(0);
  }

}

void computeDescriptor(DARY *dx, DARY *dy, uint DESC,
		       FeatureDescriptor *cor, vector<FeatureDescriptor *> &desc, int noangle){

  vector<float> angles;
  DARY *gpatch = new DARY(dx->y(),dx->x(),FLOAT1);
  DARY *opatch = new DARY(dx->y(),dx->x(),FLOAT1);
  gradAngle(dx,dy,gpatch,opatch);  
  if(!noangle)
    computeHistAngle(gpatch,opatch,angles);
  else angles.push_back(0.0);

  for(uint i=0;i<angles.size();i++){
    FeatureDescriptor *ds= new FeatureDescriptor();
    ds->copy(cor);
    ds->setAngle(angles[i]);
    compDesc(dx,dy,gpatch,opatch,ds,DESC);
    desc.push_back(ds);
    //{gpatch->writePNG("norm.png");cout << angles[i]<< " " << i << " angle  of " << angles.size() << endl;getchar();}
  }      
  delete gpatch;delete opatch;
  angles.clear();
}

/*****************************************************************************************/

/**********************COMPUTE DESCRIPTORS FROM EXTERNAL DETECTORS************************/
void computeAffineDescriptor( DARY *imgbs, DARY *patch, float scal, int DESC,
			FeatureDescriptor *cor, vector<FeatureDescriptor *> &desc){

  if(cor->getAngle()<M_PI && cor->getAngle()>-M_PI){  
    normalize(patch,patch->x()>>1,patch->y()>>1,patch->x()>>1);
    FeatureDescriptor *ds= new FeatureDescriptor();
    ds->copy(cor);
    compDesc(patch,ds, DESC);
    desc.push_back(ds);
    return;
  }   

  DARY * grad = new DARY(patch->x(),patch->x(),FLOAT1);
  DARY * ori = new DARY(patch->x(),patch->x(),FLOAT1);
  vector<float> angles;
  gradAngle(patch, grad, ori);
  computeHistAngle(grad,ori,angles);
  if(angles.size()==0)angles.push_back(0);
  delete grad;delete ori;
  for(uint i=0;i<angles.size();i++){
    patch->interpolate(imgbs,imgbs->x()>>1,imgbs->y()>>1,1.0/scal,1.0/scal,-180*angles[i]/M_PI); // normalization done
    normalize(patch,patch->x()>>1,patch->y()>>1,patch->x()>>1);
    FeatureDescriptor *ds = new FeatureDescriptor();
    ds->copy(cor);
    ds->setAngle(angles[i]);
    compDesc(patch,ds, DESC);
    desc.push_back(ds);
    //if(scal>0){patch->writePNG("norm.png");cout << angles[i]<< " " << i << " angle  of " << angles.size() << endl;getchar();}
  }      
  angles.clear();
}

DARY* normalizeAffine(DARY *image, float x, float y, float c_scale, 
		      float angle, float mi11, float mi12,float mi21, float mi22, float &scal,
		      DARY *patch, float DESC_SCALE){
  
  int imb_size=2*(int)(1.414*DESC_SCALE*c_scale)+1;//9 pixels margin for blur and rotation
  scal=(2.0*DESC_SCALE*c_scale+1)/((float)patch->x());//scale for sampling without margin
  float lecos=1;
  float lesin=0; 

  if(angle<M_PI && angle>-M_PI){// if angle is already  estimated 
    lecos=cos(-angle);
    lesin=sin(-angle);
  } 
     
  float m11=(mi11*lecos-mi12*lesin);
  float m12=(mi11*lesin+mi12*lecos);
  float m21=(mi21*lecos-mi22*lesin);
  float m22=(mi21*lesin+mi22*lecos); 
  
  //  smooth before sampling
  //cout <<" x "<< x << " y "<< y << " cs "<< c_scale << " sc " << scal<< " imb " <<  imb_size<<" " <<  1.5*scal << " "<< PATCH_SIZE*scal<<  endl;
    //cout << angle << " m " <<mi11<< " " << mi12 << " " << mi21 << " " << mi22 << endl;  
  DARY * imgb=new DARY(imb_size,imb_size,FLOAT1);
  imgb->interpolate(image,x,y,m11,m12,m21,m22);//imgb->writePNG("bsnorm.png");//normalize affine with scale=1 
  if(scal>1.3){
    DARY * imgbs=new DARY(imb_size,imb_size,FLOAT1);
    smooth(imgb,imgbs,(scal)/1.3);//imgbs->writePNG("smnorm.png");//smooth 
    imgb->set(imgbs);
    //cout << c_scale << " " << 1.5*scal << endl; 
    delete imgbs;
  }
 
  patch->interpolate(imgb,imb_size>>1,imb_size>>1,1.0/scal,1.0/scal,0.0);//not rotated
  //patch->writePNG("patch.png"); getchar();
  return imgb;
}

 


/**********************COMPUTE DESCRIPTORS FROM EXTERNAL DETECTORS************************/
void computeAffineDescriptor( DARY *image, DARY *patch, int DESC, FeatureDescriptor *cor, vector<FeatureDescriptor *> &desc, float DESC_SCALE){
     
  float angle=cor->getAngle();
  float scal; 

  // DARY *imgbs= normalizeAffine(image, cor, patch, scal);
  //cor->Cout();
  DARY *imgbs = normalizeAffine(image, cor->getX(), cor->getY(), cor->getScale(), 
			       angle, cor->getMi11(), cor->getMi12(), cor->getMi21(),cor->getMi22(),
				scal,patch,DESC_SCALE);  
  computeAffineDescriptor(imgbs, patch, scal,DESC,cor,desc);
  delete imgbs;
  //estimate orientation angle
 
} 
 
/**********************COMPUTE DESCRIPTORS FROM EXTERNAL DETECTORS************************/
void computeAffineDescriptors(DARY *image,  vector<FeatureDescriptor *> &desc, int DESC, float DESC_SCALE){
  if(desc.size()==0)return;
  int size=desc[0]->getRadius();
    initPatchMask(size);
    DARY * patch = new DARY(size,size,FLOAT1);   
    vector<FeatureDescriptor *> tmpdesc;
    for(uint c=0;c<desc.size();c++){
      computeAffineDescriptor( image, patch, DESC, desc[c], tmpdesc,DESC_SCALE);
	if(!(c%100))cout << "\rdescriptor "<< c<< " of "<< desc.size() << "    " << flush;
    }
    delete patch;
    desc.clear();
    desc=tmpdesc;
    cout << desc.size() << endl;
    //tmpdesc.clear();
}
 

void pca(vector<FD*> &features, float *mvec, float *base, int newdim){
  float maxf=0,minf=0;
/*  float *var = new float[newdim];
  bzero(var,newdim*sizeof(float));*/
  FD *f;
  float *vec,uval;
  if(features.size() > 100000)cout << endl;
  if(newdim>0)
  for(uint i=0;i<features.size();i++){
    if(!(i%100000) && features.size() > 100000)cout << "\rdone " << i << " of " << features.size() << "      "<< flush;
    f=features[i];
    f->pca(newdim,mvec,base);
    //f->normalizeVect();
    vec=f->getVec();
    for(int j=0;j<newdim;j++){
         uval = vec[j]+127;
         uval = (255 < uval) ? 255 : uval;
         uval = (0 > uval) ? 0 : uval;
         vec[j]=(int)(uval);
    }
    //getchar();
  }   
/*  for(int j=0;j<newdim;j++){
    var[j]=(var[j]/features.size());cout << var[j]<< " ";
  }
  cout << endl << maxf <<" max min feat " <<minf<< endl;*/
    if(features.size() > 100000)cout << endl;


}


void meanVarBase(vector<FD*> &features, char * filename){
  cout << endl<< "Estimating audio mean variance basis and saving in " << filename<< endl;
  cout << "Estimating pca base for "<<  features.size() << " features of dim " <<features[0]->getSize()  <<  endl;
  if(features.size()<100){
    cout << "Too few features for pca"<< endl;
    return; 
  }
  int dim=features[0]->getSize();
 
  float *stdvar =new float[dim];
  float *mvec =new float[dim];
  bzero(mvec,dim*sizeof(float));
   
  float *vec;
  float fnb=0;
  for(uint i=0;i<features.size();i++){
    vec=features[i]->getVec();
    if(fabs(vec[0])!=0.0 ){// && fabs(vec[38])<100.0 &&  fabs(vec[12])<100.0 && fabs(vec[25])<100.0){
      for(int d=0;d<dim;d++){
        mvec[d]+=vec[d];
      }
      fnb++;
    }// else features[i]->Cout(39);
  }
  cout << "FNB " << fnb << endl;
    ofstream output(filename);
  if(!output){
    cout<< "Saving error "<< endl;
    exit(0);
  }
  output <<  dim << endl;       

  for(int d=0;d<dim;d++){
    mvec[d]=mvec[d]/fnb;
    output <<  mvec[d] <<  endl; 
    //  cout << mvec[d] << " ";
  }    
  //  cout << endl;getchar();

  for(uint i=0;i<features.size();i++){
    vec=features[i]->getVec();
    for(int d=0;d<dim;d++){
      vec[d]-=mvec[d];
    }    
  }

  Matrix cov(dim,dim,0.0);
  float *f1;
  for(uint i=0;i<features.size();i++){
    f1=features[i]->getVec();
     if( fabs(f1[0])!=0.0)//  &&  fabs(f1[38])<100.0 && fabs(f1[12])<100.0 && fabs(f1[25])<100.0)
    for(int d1=0;d1<dim;d1++){
      for(int d2=d1;d2<dim;d2++){
	cov.tabMat[d1+1][d2+1]+=((f1[d1])*(f1[d2]));
      }
    }
  }
  
  for(int d1=1;d1<=dim;d1++){
      stdvar[d1-1]=sqrt(cov.tabMat[d1][d1]/fnb);
      output << stdvar[d1-1] <<  endl; 
          //  cout << stdvar[d1-1] << " ";
    
  }
  output.close();
  //cout << cov<< endl;
  delete []stdvar;
  delete []mvec;
}

void meanBase(vector<FD*> features, float *mvec, int dim){
  bzero(mvec,dim*sizeof(float));
  float *vec;
  for(uint i=0;i<features.size();i++){
    if(!(i%1000))cout << "\rmean " << i << " of "<< features.size()<< "  "<< flush;
    vec=features[i]->getVec();
    for(int d=0;d<dim;d++){
      mvec[d]+=vec[d];
     }
  }
  
  float fnb=(float)features.size();
  for(int d=0;d<dim;d++){
    mvec[d]=mvec[d]/fnb;
    //  cout << mvec[d] << " ";
  }    

  cout << endl;
  
}

void varianceBase(vector<FD*> features, float *vvec, int dim){
  bzero(vvec,dim*sizeof(float));
  float *f1;
  for(uint i=0;i<features.size();i++){
    if(!(i%1000))cout << "\rvariance " << i << " of "<< features.size()<< "  "<< flush;
    f1=features[i]->getVec();
    for(int d1=0;d1<dim;d1++){
      vvec[d1]+=f1[d1]*f1[d1];
    }
  }
}

void pcaBase(vector<FD*> &features, float *eigen, float *mvec, float *base, int newdim){
  cout << "Estimating pca base for "<<  features.size() << " features of dim " <<features[0]->getSize()  <<  endl;
  if(features.size()<100){
    cout << "Too few features for pca"<< endl;
    return; 
  }
  int dim=features[0]->getSize();
  if(dim<newdim){
    cout << "pcaBase dim > newdim "<<dim<< " new  " << newdim <<  endl;exit(0);
  }
  // compute mean  

  //  for(uint d=0;d<dim;d++){
  //    mvec[d]=0;
  //  } 

  meanBase(features, mvec,  dim);
  
  float *vec;

  for(uint i=0;i<features.size();i++){
    vec=features[i]->getVec();
    for(int d=0;d<dim;d++){
      vec[d]-=mvec[d];
    }    
  }

  float *vvec=new float[dim];
  varianceBase(features, vvec,  dim);


  Matrix cov(dim,dim,0.0);
  float *f1;
  for(uint i=0;i<features.size();i++){
    if(!(i%1000))cout << "\rco-variance, " << i << " "<< features.size()<< "  "<< flush;
    f1=features[i]->getVec();
    for(int d1=0;d1<dim;d1++){
      for(int d2=d1+1;d2<dim;d2++){
	cov.tabMat[d1+1][d2+1]+=((f1[d1])*(f1[d2]));
      }
    }
  }
  
  for(int d1=1;d1<=dim;d1++){
    cov.tabMat[d1][d1]=vvec[d1-1];
    for(int d2=d1+1;d2<=dim;d2++){
      cov.tabMat[d2][d1]=cov.tabMat[d1][d2];
    }
  }
  //cout << cov<< endl;
  
  Matrix V(dim,dim);
  Vector d(dim);
  //cov.write("cov.mat",1);
  cov.jacobi(d, V);
  //d.write("d.mat",1);
  //V.write("V.mat",1);
  for(int d1=1;d1<=dim;d1++)eigen[d1-1]=d(d1)/d(1);
  int cnt=0;
  for(int d1=1;d1<=newdim;d1++){
    for(int d2=1;d2<=dim;d2++){
      base[cnt++]=V.tabMat[d2][d1];
    }
  }

  float *outvec = new float[newdim];
  for(uint f=0;f<features.size();f++){
    if(!(f%1000))cout << "\rprojecting features, " << f << " of "<< features.size()<< "  "<< flush;
    vec=features[f]->getVec();
    //for(int v=0;v<newdim;v++)outvec[v]=0;
    bzero(outvec,newdim*sizeof(float));
    uint cnt=0;
    for(int i=0;i<newdim;i++){
      for(int v=0;v<dim;v++,cnt++){
	outvec[i] += vec[v]*base[cnt];   
      }    
    }
    //for(int i=0;i<newdim;i++)vec[i]=outvec[i];     
    memcpy(vec, outvec, newdim*sizeof(float));
    features[f]->setSize(newdim);
    //features[f]->normalizeVect();
  }
  delete []outvec;
  
/*  float *mmvec = new float[newdim];
  meanBase(features, mmvec,  newdim);*/
  varianceBase(features, vvec,  newdim);

  //float mvar=0;
  for(int i=0;i<newdim;i++){
    //mvar+=vvec[i];
    vvec[i]=50.0/sqrt(vvec[i]/((float)features.size()));//50 to have the values distributed within a byte
  }
  //mvar/=((float)newdim*features.size());

/*  for(int i=0;i<newdim;i++){
    cout << i << "  "  << mmvec[i] << " " << (vvec[i]/((float)features.size()))<< endl;
  }
  delete []mmvec;*/
  cnt=0;
  for(int i=0;i<newdim;i++){
    for(int v=0;v<dim;v++,cnt++){
        base[cnt]= vvec[i]*base[cnt];   
    }    
  }
  delete []vvec;
  /*
  writeFeatures(features, "features.max.pca",1);

  Matrix cov2(newdim,newdim,0.0);
  for(uint i=0;i<features.size();i++){
    f1=features[i]->getVec();
    for(uint d1=0;d1<newdim;d1++){
      for(uint d2=d1;d2<newdim;d2++){
	cov2.tabMat[d1+1][d2+1]+=((f1[d1])*(f1[d2]));
      }
    }
  }
  
  for(uint d1=1;d1<=newdim;d1++){
    for(uint d2=d1+1;d2<=newdim;d2++){
      cov2.tabMat[d2][d1]=cov2.tabMat[d1][d2];
    }
  }

  cov2.write("cov2.mat",1);
  */


  //
}

int loadPCAVectors(const char *filename, float *&eigenv,  float *&vmean, uint &mean_dim, float *&vbase, uint &base_dim, uint &desc_dim){
  cout << "Loading projection vectors from " <<filename << flush;
  ifstream input(filename);
  if(!input){
    cout << "No pca vectors in " <<filename << endl;
    return 0;
  }
  input >> desc_dim;  
  input >> mean_dim;  

  if(desc_dim==0){
    input.close();
    return 1;
  }

  eigenv=new float[mean_dim];

  vmean=new float[mean_dim];
  float val;
  for(uint i=0;i<mean_dim;i++){
    input >> val;
    eigenv[i]=val;
  }

  for(uint i=0;i<mean_dim;i++){
    input >> val;
    vmean[i]=val;
  }
  input >> base_dim;  
  vbase=new float[base_dim];  
  for(uint i=0;i<base_dim;i++){
    input >> val;
    vbase[i]=val;
  }
  
  input.close();
  cout << " from dim "<< mean_dim << " to " << desc_dim << endl;
  return 1;
}

void savePCAVectors(const char *filename, float *eigenv,  float *vmean, uint size1 ,float *vbase, uint size2, uint desc_dim){
  cout << "Saving pca vectors in " <<filename << endl;

  ofstream output(filename);
  if(!output){
    cout<< "Saving error "<< endl;
    exit(0);
  }
  output <<  desc_dim << endl;    
  output <<  size1 << endl;    

  for(uint i=0;i<size1;i++){
    output  <<  eigenv[i] << endl;
  }

  for(uint i=0;i<size1;i++){
    output  <<  vmean[i] << endl;
  }

  output << size2 << endl;   
  for(uint i=0;i<size2;i++){
    output << vbase[i] << endl; 
  }

  output.close();
}


void computeColorDesc(DARY *cimg, FD *f){

  float *vec=f->getVec();
  uint size=f->getSize();
  bzero(vec,size*sizeof(float));
  uint sz=(uint)(0.5+f->getScale());
  //sz=(sz>100)?100:sz;
  f->setType(((DETECTOR & f->getType())| DCOLOR));

  int x=(int)(0.5+f->getX());
  int y=(int)(0.5+f->getY());
  
  int xs=x-sz;
  xs=(xs<0)?0:xs;
  int xe=x+sz;
  xe=(xe<(int)cimg->x())?xe:cimg->x();
  int ys=y-sz;
  ys=(ys<0)?0:ys;
  int ye=y+sz;
  ye=(ye<(int)cimg->y())?ye:cimg->y();
  
  float rf,gf,bf,dr,dg,db,rweight,gweight,bweight;
  int ri,gi,bi,rindex,gindex,bindex;
  int cell=5,cell2=25;
  int step=1;
  //if(size>20)step=2;
  for(int i=xs;i<xe;i+=step){
    for(int j=ys;j<ye;j+=step){      
      rf=cimg->belr[j][i]/64.0;
      gf=cimg->belg[j][i]/64.0;
      bf=cimg->belb[j][i]/64.0;     
      
      ri=(int)rf;
      gi=(int)gf;
      bi=(int)bf;

      //vec[ri*cell2+gi*cell+bi]+=1;
      
      dr=rf-ri;
      dg=gf-gi;
      db=bf-bi;
    
      for (int r = 0; r < 2; r++) {
	rindex = (ri + r)*cell2;
	rweight = ((r == 0) ? 1.0 - dr : dr);
	
	for (int g = 0; g < 2; g++) {
	  gindex = rindex+((gi + g)*cell);
	  gweight = rweight*((g == 0) ? 1.0 - dg : dg);
	  for (int b = 0; b < 2; b++) {
	    bindex = gindex+bi + b;
	    bweight = gweight*((b == 0) ? 1.0 - db : db);
	    //if(bindex>128)cout << bindex << " " << gindex<< " " << rindex<< " "<< ri <<" " << gi << " "<< bi <<   endl;
	    vec[bindex]+=bweight; 
	  }
	}
      }
      
    }
  } 
  
  normalizeFeature(vec,  size);
  //f->Cout(128);
}



void AddColorSample(float *index, float rgbf, float rf, float gf, float rx, float cx,
		    float rpos, float cpos, int colorSize, int locSize)
{
  int r, c, gc, rc, rgbi, lc,  ri, ci,  rindex, cindex, gindex, rgindex,rfi,gfi;
  float  rfrac, cfrac, rffrac,  gffrac,rgbffrac,rweight, cweight, gweight,rgweight,rgbweight;
   
   
   ri = (int)((rx >= 0.0) ? rx : rx - 1.0);  /* Round down to next integer. */
   ci = (int)((cx >= 0.0) ? cx : cx - 1.0);
   rfi = (int)rf;//red
   gfi = (int)gf;//green
   rgbi=(int)rgbf;//lumin
   rfrac = rx - ri;         /* Fractional part of location. */
   cfrac = cx - ci;
   rffrac = rf - rfi;         /* Fractional part of color. */
   gffrac = gf - gfi;
   rgbffrac = rgbf - rgbi;
   
   //    cout << "in r " << " rf "<<  rf << " gf " << gf <<" rx " << rx << " cx " << cx  << " rfi "<<  rfi << " gfi " << gfi <<" ri " << ri << " ci " << ci <<  " rfrac " << rfrac   << " cfrac " << cfrac <<   " rffrac " << rffrac << " gffrac " << gffrac <<   endl;

   /* Put appropriate fraction in each of 8 buckets around this point
      in the (row,col,ori) dimensions.  This loop is written for
      efficiency, as it is the inner loop of key sampling. */
   for (r = 0; r < 2; r++) {
     rindex = (ri + r);
      if (rindex >=0 && rindex < locSize) {
         rweight =  ((r == 0) ? 1.0 - rfrac : rfrac);
         rindex = rindex*locSize;

         for (c = 0; c < 2; c++) {
            cindex = ci + c;
            if (cindex >=0 && cindex < locSize) {
	      cindex=(rindex+cindex)*colorSize;
               cweight = rweight * ((c == 0) ? 1.0 - cfrac : cfrac);

               for (gc = 0; gc < 2; gc++) {
                  gindex = gfi + gc;
                  if (gindex < colorSize && gindex>=0){		    
		    gweight = cweight * ((gc == 0) ? 1.0 - gffrac : gffrac);
		    gindex=(cindex+gindex)*colorSize;

		    for (rc = 0; rc < 2; rc++) {
		      rgindex = rfi + rc;
		      if (rgindex < colorSize && rgindex >=0){
			rgweight=gweight*((rc == 0) ? 1.0 - rffrac : rffrac);
			rgindex=(gindex+rgindex)*colorSize;

			for (lc = 0; lc < 2; lc++) {			  
			  if (rgbi+lc < colorSize && rgbi+lc>=0){
			    rgbweight = rgweight * ((lc == 0) ? 1.0 - rgbffrac : rgbffrac);		      
			    index[rgindex+rgbi+lc]+=rgbweight;
			  }
			}
		      }
		    }
		  }
               }
            }  
         }
      }
   } 
} 


void rgb_to_hsv(float r, float g, float b, float &h, float &s, float &v){

  float rgbmax = ( r> g)?r:g;
  rgbmax = ( b> rgbmax)?b:rgbmax;
  float rgbmin=( r< g)?r:g;
  rgbmin = ( b< rgbmin)?b:rgbmin;
  v = rgbmax/255.0;
  s = ( rgbmax - rgbmin ) / rgbmax;

  float  rc = ( rgbmax - r ) / ( rgbmax - rgbmin );
  float  gc = ( rgbmax - g ) / ( rgbmax - rgbmin );
  float  bc = ( rgbmax - b ) / ( rgbmax - rgbmin );

  if ( r == rgbmax )  
    h = bc - gc;
  else if ( g == rgbmax ) 
    h = rc - bc;
  else
    h =  gc - rc;
} 
void rgb_to_yuv(float r, float g, float b, float &y, float &u, float &v){
  y = 0.299*r +0.587*g +0.114*b;
  u = 0.5 - 0.1687*r -0.3313*g +0.5*b;
  v = 0.5 + 0.5*r-0.4187*g -0.0813*b;
  
}

void computeColorSift(DARY *cimg, FD *f){

  int   locSize= 2;
  int   colorSize= 3; 
  f->allocVec(locSize*locSize*colorSize*colorSize*colorSize);
  float *vec=f->getVec();
  uint size=f->getSize();
  int sz=(int)(0.5+f->getScale());
  //sz=(sz>100)?100:sz;
  f->setType(((DETECTOR & f->getType())| DCOLOR));

  float angle=f->getAngle();
  int x=(int)(0.5+f->getX());
  int y=(int)(0.5+f->getY());
  
  int xs=(x<sz)?x:sz;
  int xe=(x+sz<(int)cimg->x())?sz:((int)cimg->x()-x-1);
  int ys=(y<sz)?y:sz;
  int ye=(y+sz<(int)cimg->y())?sz:((int)cimg->y()-y-1);
  
  float rf,gf,rgb,Y,U,V;
  int ri,gi,bi;
  int step=1;
  

  float spacing, rpos, cpos, rx, cx;
  spacing = (locSize + 1) / (2.0*sz);
  float sine = sin(-angle);
  float cosine = cos(-angle);
  
  //if(size>20)step=2;
  for(int i=-xs;i<=xe;i+=step){
    for(int j=-ys;j<=ye;j+=step){      
      
       rpos = (cosine*i + sine*j) * spacing;
       cpos = (-sine*i + cosine*j) * spacing;
       rx = rpos + (locSize - 1) / 2.0;
       cx = cpos + (locSize - 1) / 2.0;
       if (rx > -1.0 && rx < (float) locSize  &&
	   cx > -1.0 && cx < (float) locSize){
	 ri=cimg->belr[y + j][x + i];
	 gi=cimg->belg[y + j][x + i];
	 bi=cimg->belb[y + j][x + i];
	 //ri=(255.0*rand())/((float)RAND_MAX);
	 //gi=(255.0*rand())/((float)RAND_MAX);
	 //bi=(255.0*rand())/((float)RAND_MAX);

	 
	 rgb=ri+gi+bi;	 
	 rf=colorSize*ri/rgb;
	 gf=colorSize*gi/rgb;
	 rgb=colorSize*rgb/765.0;

	 rgb_to_yuv( ri/255.0,  gi/255.0,  bi/255.0, Y, U ,V);
	 rgb=colorSize*Y;
	 rf=colorSize*U;
	 gf=colorSize*V;

	 AddColorSample(vec,rgb , rf, gf, rx, cx,  rpos, cpos, colorSize, locSize);
       }       
    }
  }   
  
  

   normalizeFeature(vec,  size);
  //f->Cout(128);
}



void initSegmenSift(vector<segment *> epoints, vector<FD *> &features){
  float mx=0,my=0,sx=0,sy=0,sxy=0,l1=1,el=1,l2=1,ea=1,ssx,ssy,maxpf,qssx,qssy;
  int maxp;
  for(uint i=0;i<epoints.size();i++){
    //cout << "psize  "<< epoints[i]->points.size() << endl;
      features.push_back(new FD());
      mx=0;my=0;sx=0;sy=0;sxy=0;
      for(uint p=0;p<epoints[i]->points.size();p++){
        mx+=epoints[i]->points[p].a;
        my+=epoints[i]->points[p].b; 
      }
      mx/=epoints[i]->points.size();
      my/=epoints[i]->points.size();
      maxp=0;maxpf=0;
      for(uint p=0;p<epoints[i]->points.size();p++){
        ssx=(epoints[i]->points[p].a-mx);
        ssy=(epoints[i]->points[p].b-my);
        //epoints[i]->points[p].a=ssx;
        //epoints[i]->points[p].b=ssy;
        qssx=square(ssx);
        qssy=square(ssy);
        sx+=(qssx);
        sy+=(qssy);
        sxy+=((ssy)*(ssx));
        if(maxpf<qssx){
          maxpf=qssx;
          //maxp=p;
        }
        if(maxpf<qssy){
          maxpf=qssy;
          //maxp=p;
        }
        
   //   cout << mx << " my " << my << " a " <<  epoints[i]->points[p].a<< " a-mx " <<   (epoints[i]->points[p].a-mx) << " b " <<  (epoints[i]->points[p].b) <<" b-my " << (epoints[i]->points[p].b-my) <<  endl; 
        //if(!(p%100))getchar();
      }
      sx/=epoints[i]->points.size();
      sy/=epoints[i]->points.size();
      sxy/=epoints[i]->points.size();
      
      Matrix U(2,2,0.0),D,Vi,V;
      U(1,1)=(sx);
      U(1,2)=sxy;
      U(2,1)=sxy;
      U(2,2)=(sy);
      U.svd(Vi,D,V);
      
      //D(1,1)=(1.0/(D(1,1)));
     // D(2,2)=(1.0/(D(2,2)));
      float a=(sqrt(D(2,2)*D(1,1)));
    //cout << D(1,1)<< " " <<  D(2,2)<< "  " << tmp1<< endl;
      D.tabMat[2][2]/=a;
      D.tabMat[1][1]/=a;
      U=V*D*V.transpose();
      el=D.tabMat[2][2]/D.tabMat[1][1];
      features[i]->setScale(1.5*sqrt(maxpf));
//      features[i-1]->setScale(sqrt(a));
  //    fetaures[i-1]->setScale(a);
      //invSqrRoot(sx,sxy,sy,l2,l1,ea);		
      el=l2/l1;
      features[i]->setX_Y((int)mx,(int)my);
 //     fetaures[i-1]->setMi(sx,sxy,sxy,sy,l1,el,ea);
      features[i]->setMi(U(1,1),U(1,2),U(2,1),U(2,2),D(1,1),el,ea);
 
    
     // cout << i << " " << mx << " x " << my << "  sx "<<sx << " sy " << sy <<" sxy " << sxy <<  " l1 " << l1<< " l2 " << l2 << " scale " <<   fetaures[i-1]->getScale()<< " " << epoints[i]->points.size() << endl;
     // fetaures[i-1]->Cout();getchar();
  }
}

void histAngle(vector<edge> points, float **grad, float **ori, float max_angle, vector<float> &angle){
  int oriBins=18;
  float *hist=new float[oriBins];
  for (int i = 0; i < oriBins; i++)
    hist[i] = 0.0;
 
  float gval=0,oval=0,fbin,dbin;
  int bin;
  for(uint i=0;i<points.size();i++){
    gval = grad[points[i].b][points[i].a];  
    oval = M_PI/2+ori[points[i].b][points[i].a];  
    oval=(oval<0)?0:oval;
    fbin =  (oriBins * oval / max_angle);
    bin=(int)floor(fbin);
    dbin=fbin-bin;
    //assert(bin >= 0 && bin <= OriBins);
    //cout << bin << " " << fbin <<  " " << oval << " " << oriBins << " " << hist[bin] << " " << points[i].b<< " " << points[i].a  << endl;
    bin = (bin < oriBins)?bin:(0);
    hist[bin] +=  (1-dbin)*gval;
    bin = (bin+1 < oriBins)?bin+1:(0);
    hist[bin] +=  (dbin) * gval;
  }
  findDominantAngles(hist, oriBins, angle, max_angle, max_angle);
  
  delete []hist;
}

void extractSegmentSift(DARY *sr,DARY *sg, DARY *sb, vector<FD *> &features,  Params *params){
 
  int width=sr->x();
  int height=sr->y(); 
  float segmenation_k=params->getValue("segmentation_col_threshold.int");
  int min_size=params->getValue("segmentation_min_size.int");
  float grad_thres=params->getValue("segmentation_min_size.int");
  float bound_thres=params->getValue("segmentation_min_size.int");

  DARY *igrad= new DARY(height,width,FLOAT1,0.0);
  DARY *iori= new DARY(height,width,FLOAT1,0.0);
  DARY *labout= new DARY(height,width,FLOAT1,0.0);
  
  int num_ccs=0;
  vector<segment *> epoints;
  segment_image(sr, sg, sb, segmenation_k, grad_thres, bound_thres, min_size, num_ccs, igrad, iori, labout, epoints);
  //labout->writePNG("labout.png");
  vector<FD *> desc;
  initSegmenSift(epoints, desc);
  int lab=0;
  float **labs=labout->fel;
  int oriSize=8;
  int locSize=4;
  float max_angle=M_PI;
  vector<float> angle;
  vector<FD *> angdesc;
  int dsize=oriSize*locSize*locSize;
  for(uint i=0;i<epoints.size();i++){
    histAngle(epoints[i]->points, igrad->fel, iori->fel,  max_angle, angle);
    //cout << i << " " << angle.size()<< " " << epoints[i]->points.size()<< endl;
    for(uint a=0;a<angle.size();a++){
    //cout << i << " " << a << " " <<  angle[a]<< endl;
      FD *fd = new FD();
      fd->copy(desc[i]);
      fd->setAngle(angle[a]-M_PI/2);
      fd->allocVec(dsize);
      fd->setType(DCLBP);
      KeySample(fd, epoints[i]->points, igrad->fel, iori->fel,  oriSize,  locSize,  max_angle);
      angdesc.push_back(fd);
    }
    angle.clear();
  }
  //getchar();
 // iori->writePNG("ori.png");
  features.insert(features.begin(),angdesc.begin(),angdesc.end());
  for(uint i=0;i<epoints.size();i++){
    epoints[i]->points.clear();
  }
  epoints.clear();
  angdesc.clear();
  deleteDescriptors(desc);
  delete igrad;
  delete iori;delete labout;
}

void extractSegmentSift(DARY *image, vector<FD *> &features,  Params *params){
  unsigned long type=(unsigned long)params->getValue("feature_type.int");
  //cout << " " << DCLBP << " " << type << "  "<< ( type & DCLBP) << "  "<< ( type & DSIFT) << "  "  <<(DSEDGE|DSIFT|DCLBP)<< endl;
  if(( type & DCLBP)==0)
    return;
//   if(image->getType()==UCHAR)image->toCOLOR();
//   else if(image->getType()==CUCHAR);
//   else {
//     cout <<  "extractColorSift image error " <<endl;
//     exit(0); 
//   }

  int width=image->x();
  int height=image->y(); 
  DARY *r= new DARY(height,width,FLOAT1,0.0);
  DARY *g= new DARY(height,width,FLOAT1,0.0);
  DARY *b= new DARY(height,width,FLOAT1,0.0);
  DARY *sr= new DARY(height,width,FLOAT1,0.0);
  DARY *sg= new DARY(height,width,FLOAT1,0.0);
  DARY *sb= new DARY(height,width,FLOAT1,0.0);
  image->get(r,g,b);
  smooth5(r,sr);
  smooth5(g,sg);
  smooth5(b,sb);
  int num_ccs=0;
  vector<FD *> feats;

  extractSegmentSift(sr,sg, sb, feats,  params);
  features.insert(features.begin(),feats.begin(),feats.end());
  cout << "segment feats "<<feats.size()<< endl;
  feats.clear();  
  
   height=height>>1;
   width=width>>1;
   DARY *srs= new DARY(height,width,FLOAT1,0.0);
   DARY *sgs= new DARY(height,width,FLOAT1,0.0);
   DARY *sbs= new DARY(height,width,FLOAT1,0.0);
   DARY *ssrs= new DARY(height,width,FLOAT1,0.0);
   DARY *ssgs= new DARY(height,width,FLOAT1,0.0);
   DARY *ssbs= new DARY(height,width,FLOAT1,0.0);
  
  float scale=2.0;
  srs->scale(sr,scale,scale);
  sgs->scale(sg,scale,scale);
  sbs->scale(sb,scale,scale);
  

  smooth5(srs,ssrs);
  smooth5(sgs,ssgs);
  smooth5(sbs,ssbs);
  extractSegmentSift(ssrs,ssgs, ssbs, feats,  params);
  
  for(uint i=0;i<feats.size();i++){
    feats[i]->setScale(scale*feats[i]->getScale());
    feats[i]->setX_Y(scale*feats[i]->getX(),scale*feats[i]->getY());
  }
  cout << "segment feats "<<feats.size()<< endl;
  features.insert(features.begin(),feats.begin(),feats.end());
  feats.clear();

  cout << "segment features "<<features.size()<< endl;

 feats.clear();
  delete srs;
  delete sgs;
  delete sbs;
  delete ssrs;
  delete ssgs;
  delete ssbs;
  delete r;
  delete g;
  delete b;
  delete sr;
  delete sg;
  delete sb;

  
  
}


void extractColorLBP(DARY *img, DARY *gray,vector<FD *> &features,  Params *params){
  
  unsigned long type=(unsigned long)params->getValue("feature_type.int");
  //cout << " " << DCLBP << " " << type << "  "<< ( type & DCLBP) << "  "<< ( type & DSIFT) << "  "  <<(DSEDGE|DSIFT|DCLBP)<< endl;
  if(( type & DCLBP)==0)
    return;

  if(img->getType()==UCHAR1)img->convert(UCHAR3);
  else if(img->getType()==UCHAR3);
  else {
    cout <<  "extractColorLBP image error " <<endl;
    exit(0); 
  }
    
  float max_clbp_scale=params->getValue("max_clbp_scale.float");
  cout << "computing color lbp... " << flush;
  DARY *iclbp  = new DARY(img->y(),img->x(),UCHAR1);
  iclbp->set(255);
  DARY *iclbp0 = new DARY(img->y(),img->x(),UCHAR1);
  DARY *iclbp1 = new DARY(img->y(),img->x(),UCHAR1);
  DARY *iclbp2 = new DARY(img->y(),img->x(),UCHAR1);
  DARY *iclbp3 = new DARY(img->y(),img->x(),UCHAR1);
  FD *fc;
  int lbp,lbp0,lbp1,lbp2,lbp3,x=0,y=0,rad,xs,xe,ys,ye;
  int c,cc,ccr,ccb,ccg,rc,gc,bc,rl,rr,rt,rb,gl,gr,gt,gb,bl,br,bt,bb,ggl,ggr,ggtl,ggtr,ggbl,ggbr; 
  uchar **ir=img->belr;
  uchar **ig=img->belg;
  uchar **ib=img->belb;
  uchar **mlbp=iclbp->bel;
  uchar **mlbp0=iclbp0->bel;
  uchar **mlbp1=iclbp1->bel;
  uchar **mlbp2=iclbp2->bel;
  uchar **mlbp3=iclbp3->bel;
  float **igray=gray->fel;
  float *vec;
  int dsize=128,max=0;
  uint fsize=features.size();
  for(uint i=0;i<fsize;i++){
    if(features[i]->getScale()<max_clbp_scale && !(((int)features[i]->getX())==x && ((int)features[i]->getY())==y)){
      fc=new FD();
      fc->copy(features[i]);
      x=fc->getX(); 
      y=fc->getY();
      fc->allocVec(dsize);
      rad=(int)(fc->getScale()*10);
      xs=((x-rad)<2)?2:(x-rad); 
      ys=((y-rad)<2)?2:(y-rad);
      xe= ((x+rad)<img->x()-2)?(x+rad):(img->x()-3); 
      ye= ((y+rad)<img->y()-2)?(y+rad):(img->y()-3); 
      vec=fc->getVec();
      int nb=0;
      for(int jj=ys;jj<ye;jj++){
        for(int ii=xs;ii<xe;ii++){
          if(mlbp[jj][ii]==255){
            
            c=igray[jj][ii];
            cc=(c<128)?0:1;
            rc=(ir[jj][ii]<c)?0:1;
            gc=(ig[jj][ii]<c)?0:1;
            bc=(ib[jj][ii]<c)?0:1; 
                     
            ccr=(ir[jj][ii])>>6;//64
            ccg=(ig[jj][ii])>>6;//64
            ccb=(ib[jj][ii])>>6;//64
            //if(ccr>=4 || ccg>=4 || ccb>=4){cout << "COLOR ERROR ccr " << ccr << " ccg " << ccg << " ccb " << ccb << " " << ir[jj][ii] << " " <<ig[jj][ii] << ib[jj][ii] <<  endl;getchar();}
           
            ggl=(igray[jj][ii-2]<c)?0:1;
            ggr=(igray[jj][ii+2]<c)?0:1;
            ggtl=(igray[jj-2][ii-1]<c)?0:1;
            ggtr=(igray[jj-2][ii+1]<c)?0:1;
            ggbl=(igray[jj+2][ii-1]<c)?0:1;
            ggbr=(igray[jj+2][ii+1]<c)?0:1;
            
            
            rl=(ir[jj][ii-2]<ir[jj][ii])?0:1;
            rr=(ir[jj][ii+2]<ir[jj][ii])?0:1;
            rt=(ir[jj-2][ii]<ir[jj][ii])?0:1;
            rb=(ir[jj+2][ii]<ir[jj][ii])?0:1;
           
            gl=(ig[jj][ii-2]<ig[jj][ii])?0:1;
            gr=(ig[jj][ii+2]<ig[jj][ii])?0:1;
            gt=(ig[jj-2][ii]<ig[jj][ii])?0:1;
            gb=(ig[jj+2][ii]<ig[jj][ii])?0:1;
           
            bl=(ib[jj][ii-2]<ib[jj][ii])?0:1;
            br=(ib[jj][ii+2]<ib[jj][ii])?0:1;
            bt=(ib[jj-2][ii]<ib[jj][ii])?0:1;
            bb=(ib[jj+2][ii]<ib[jj][ii])?0:1;
           
            lbp=(ggl<<5)|(ggr<<4)|(ggtl<<3)|(ggtr<<2)|(ggbl<<1)|(ggbr);
//            lbp=(cc<<3)|(rc<<2)|(gc<<1)|(bc);
/*            lbp1=16+((rl<<3)|(rr<<2)|(rt<<1)|(rb));
            lbp2=32+((gl<<3)|(gr<<2)|(gt<<1)|(gb));
            lbp3=48+((bl<<3)|(br<<2)|(bt<<1)|(bb));*/
            lbp0=64+(ccr<<4)|(ccg<<2)|(ccb);
            mlbp[jj][ii]=lbp;
            mlbp0[jj][ii]=lbp0;
            //mlbp1[jj][ii]=lbp1;
            //mlbp2[jj][ii]=lbp2;
            //mlbp3[jj][ii]=lbp3;
            //if(lbp0>=64){cout << "COLOR ERROR ccr " << ccr << " ccg " << ccg << " ccb " << ccb << " " << ir[jj][ii] << " " <<ig[jj][ii] << ib[jj][ii] <<  endl;getchar();}
           
          }else{
            lbp =mlbp[jj][ii];
            lbp0=mlbp0[jj][ii];
            //lbp1=mlbp1[jj][ii];
            //lbp2=mlbp2[jj][ii];
            //lbp3=mlbp3[jj][ii];
          }
          vec[lbp]++;
          vec[lbp0]++;
/*          vec[lbp1]++;
          vec[lbp2]++;
          vec[lbp3]++;*/
          //if(fc->getExtr()<0){cout << "error "<< endl;getchar();}
          nb++;
          //cout << lbp<<  endl;getchar();
        }
      }
      normalizeFeature(vec, dsize);
 //     fc->setType((ADETECTOR & fc->getType())|DCLBP);
     // unsigned int t=fc->getType();
      fc->setType((ADETECTOR&fc->getType())|DCLBP);

      //cout << fc->getType() << "  "<< nb << "  ";fc->Cout(dsize);getchar();
      features.push_back(fc);
    }
  }
  cout  << (features.size()-fsize) << endl;
  //generate features 
  
  delete iclbp;
  delete iclbp0;
  delete iclbp1;
  delete iclbp2;
  delete iclbp3;
}
