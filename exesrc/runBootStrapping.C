#include "BootStrapping.h"
#include "ArgParser.hh"
#include "ReadConfig.hh"

#include <iostream>


int main(){
  using namespace std;
  /*ArgParser a(argc,argv);
  a.addArgument("WorkspaceFile",ArgParser::required,"path to the workspace file");

  string ret;
  if(a.process(ret) !=0){
    cout << "Invalid Options:  " << ret <<endl;
    a.printOptions(argv[0]);
    return 0;
  }

  TString List = a.getArgument("WorkspaceFile");
  std::ifstream *st = new ifstream(List.Data(),std::ifstream::in);
  char line[800];
  string wsFile;
  while(st->good()){
    st->getline(line,800);
    wsFile = line;
  }
  delete st;
  */  
  string wsFile = "/mnt/hadoop/store/user/amott/Hgg2013/Hgg/workspaces/hgg_22Jan2013_R9_CIC.root";
  BootStrapping bs(wsFile);
  bs.run();
  return 0;
}
