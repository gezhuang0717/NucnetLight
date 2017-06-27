void makeP0n(TString name,TString infile,TString outfile, Int_t isUpdate){
  ifstream f(infile);
  ofstream of;
  if (isUpdate==1) of.open(outfile,ofstream::out | ofstream::app);
  else of.open(outfile);
  Int_t num;
  f>>num;
  TString temp;
  for (Int_t i=0;i<num;i++){
    of<<"single_rate\n"<<name<<"\n1\n";
    f>>temp;
    of<<temp;
    of<<"\n3\n";
    for (Int_t j=0;j<4;j++){
      f>>temp;
      of<<temp<<endl;
    }
  }
}
void makeP1n(TString name,TString infile,TString outfile, Int_t isUpdate){
  ifstream f(infile);
  ofstream of;
  if (isUpdate==1) of.open(outfile,ofstream::out | ofstream::app);
  else of.open(outfile);
  Int_t num;
  f>>num;
  TString temp;
  for (Int_t i=0;i<num;i++){
    of<<"single_rate\n"<<name<<"\n1\n";
    f>>temp;
    of<<temp;
    of<<"\n4\n";
    for (Int_t j=0;j<5;j++){
      f>>temp;
      of<<temp<<endl;
    }
  }
}
void makeP2n(TString name,TString infile,TString outfile, Int_t isUpdate){
  ifstream f(infile);
  ofstream of;
  if (isUpdate==1) of.open(outfile,ofstream::out | ofstream::app);
  else of.open(outfile);
  Int_t num;
  f>>num;
  TString temp;
  for (Int_t i=0;i<num;i++){
    of<<"single_rate\n"<<name<<"\n1\n";
    f>>temp;
    of<<temp;
    of<<"\n5\n";
    for (Int_t j=0;j<6;j++){
      f>>temp;
      of<<temp<<endl;
    }
  }
}
void makeP3n(TString name,TString infile,TString outfile, Int_t isUpdate){
  ifstream f(infile);
  ofstream of;
  if (isUpdate==1) of.open(outfile,ofstream::out | ofstream::app);
  else of.open(outfile);
  Int_t num;
  f>>num;
  TString temp;
  for (Int_t i=0;i<num;i++){
    of<<"single_rate\n"<<name<<"\n1\n";
    f>>temp;
    of<<temp;
    of<<"\n6\n";
    for (Int_t j=0;j<7;j++){
      f>>temp;
      of<<temp<<endl;
    }
  }
}
