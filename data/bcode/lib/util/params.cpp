#include "util.h"

void Params::save(const char *filename){
  cout <<"Saving parameters in "<< filename<< " ..."<< flush;
  ofstream output(filename);
  if(!output)cout << "error opening " <<filename << endl;
   output.setf ( ios::fixed, ios::floatfield );       // set hex as the basefield
   //cout.setf ( ios::showbase ); 
  for(uint i=0;i<names.size();i++){
    output << names[i] << endl;
    char *type = strchr(names[i],'.');    
    if(!strcmp(type,".int") || !strcmp(type,".float")){
      output << fvalues[i] << endl;    
      //      cout <<"........ TEST " <<names[i]<< " " <<   fvalues[i] << endl;    
      
    }else if(!strcmp(type,".char"))
      output << cvalues[i] << endl; 
  }
  output.close();
  cout << "done"<< endl;
}

void Params::load(const char *filename){
  cout <<"Loading parameters from "<< filename<< " ..."<< flush;
  ifstream input(filename);
  if(!input){
    cout << "no file " <<filename << endl;
    return;
  }
  char *bigbuf=NULL;
  double value; 
  char *bigbuf2 = new char [512];

  while(!input.eof()){    
    bigbuf = new char [512];
    input.getline(bigbuf,512);    
    char *type = strchr(bigbuf,'.');    
    if(type!=NULL){
      if(!strcmp(type,".char")){
	char *bigbuf3 = new char [512];
	input.getline(bigbuf3,512);
	put(bigbuf, bigbuf3);
      } else {
	input >> value;
	put(bigbuf, value);
	input.getline(bigbuf2,512);//read end of line
      }
      //cout << type<< " " << names[names.size()-1] << " " << fvalues[names.size()-1] <<  " " << cvalues[names.size()-1] << endl; 
    }
  }
  if(bigbuf!=NULL)delete []bigbuf2;
  delete []bigbuf;
  input.close();
  cout << "done"<< endl;

}


void Params::put(const char *name, double value){
   
  int  found=0;
  for(uint i=0;i<names.size();i++){
    if(!strcmp(names[i],name)){
      found=1;
      fvalues[i]=value;
    }
  }
  //cout << name << " " << (unsigned long)value << endl;

  if(!found){
    char *buf= new char[512];
    strcpy(buf,name);
    names.push_back(buf);
    fvalues.push_back(value);
    cvalues.push_back("null");
  }  
}
void Params::put(const char *name, const char *value){
   
  int  found=0;
  for(uint i=0;i<names.size();i++){
    if(!strcmp(names[i],name)){
      found=1;
      strcpy(cvalues[i],value);
    }
  }
  if(!found){
    char *buf= new char[512];
    strcpy(buf,name);
    names.push_back(buf);
    fvalues.push_back(0.0);
    char *buff= new char[512];
    strcpy(buff,value);
    cvalues.push_back(buff);
  }  
}

double Params::getValue(const char *name){
  uint i=0;
  int found=0;
  double value=0;
  while(i<names.size() && !found){
    if(!strcmp(names[i],name)){
      value=fvalues[i];
      found=1;
    }
    i++;
  }
  if(!found)cout << name << " not found in params "<< endl;
  return value;
}
char * Params::getString(const char *name){
  uint i=0;
  int found=0; 
  char *value=NULL;
  while(i<names.size() && !found){
    if(!strcmp(names[i],name)){
      value=cvalues[i];
      found=1;
    }
    i++;
  }
  if(!found)cout << name << " not found in params "<< endl;
  return value;
}



void loadFileNames(const char *filein, vector<char *> &filenames){

  cout << "Loading filenames from "<< filein << "... "<<flush;
  ifstream input(filein);   
  if(!input){
   cout << " no such file " << endl;
   return;
  }
  char *bigbuf;
  do{
    bigbuf = new char [512];
    input.getline(bigbuf,512);    
    if( bigbuf[ 0 ] )
       filenames.push_back(bigbuf);
  }while(!input.eof());
//  filenames.erase(filenames.end());
  input.close();
  cout << filenames.size()<< " files done."<< endl;
}

void saveFileNames(vector<char *> filenames, const char *filein ){

  cout << "Saving filenames in "<< filein << "... "<<flush;
  ofstream output(filein);
  
  for(uint i=0;i<filenames.size();i++){
    output << filenames[i]<< endl;
    
  }
  output.close();
  cout << filenames.size()<< " files done."<< endl;
}








