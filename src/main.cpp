#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include "utilities.h"
#include "term.h"
#include "orbital.h"
#include "finput.h"
#include "globals.h"
#include "work.h"
#include "tensor.h"

using namespace std;

int main(int argc, char **argv) {
  // handle input and output
  string inputfile, outputfile,
    exePath = exepath();
  int iarg=1;
  bool algo = false;
  while ( iarg<argc && argv[iarg][0]=='-') {// handle options
    if (strcmp(argv[iarg],"-h")==0 || strcmp(argv[iarg],"--help")==0) {
      cout << "quantwo <input-file> [<output-file>]" << endl;
      // print README file if exists
      ifstream readme;
      readme.open((exePath+"README").c_str());
      if (readme.is_open()) {
        string line;
        while (readme.good()) {
          getline (readme,line);
          cout << line << endl;
        }
      }
      return 0;
    } else if (strcmp(argv[iarg],"-v")==0 || strcmp(argv[iarg],"--verbose")==0) {
      if ( iarg == argc-1 || !str2num<int>(Input::verbose,argv[iarg+1],std::dec)){
        Input::verbose = 1;
      } else {
        ++iarg;
      }
    } else if (strcmp(argv[iarg],"-a")==0 || strcmp(argv[iarg],"--algo")==0) {
      // the input file is an algofile
      algo = true;
    }
    ++iarg;
  }
  if (iarg >= argc) error("Please provide an input file!");
  inputfile=argv[iarg];
  if (argc>iarg+1)
    outputfile=argv[iarg+1];
  else 
    outputfile = FileName(inputfile,true)+".tex";
  // read input
  Finput finput(exePath);
  ifstream fin;
  fin.open(inputfile.c_str());
  // save input file
  std::vector< std::string > inp;
  if (fin.is_open())
  {
    string line;
    while (fin.good()) {
      getline (fin,line);
      _xout1(line << endl);
      inp.push_back(line);
    }
  }
  else
    error("Bad input file!");
  fin.close();
  ofstream fout;
  fout.open(outputfile.c_str());
  Output myfout(fout);
  // set current ouput to fout
  MyOut::pcurout = &myfout;
  MyOut::pcurout->nlines += Input::iPars["output"]["emptylines"];
  if ( inp.size() == 0 ){
    say("Empty input file!");
    return 1;
  }
  for ( lui il = 0; il < inp.size(); ++il ){
    if ( algo ) {
    } else if ( finput.addline(inp[il]) ){
      if ( finput.sumterms().size() == 0 ){
        say("Empty equation!");
        continue;
      }
//       Input::verbose = 2;
      TermSum sum_finp(finput.sumterms());
      TermSum sum_NO;
      if ( Input::iPars["prog"]["wick"] == 0 )
        sum_NO = Q2::normalOrderPH(sum_finp);
      else 
        sum_NO = Q2::wick(sum_finp);
      TermSum sum_final1(Q2::reduceSum(sum_NO)),
        sum_final(Q2::postaction(sum_final1));
      _xout1(finput << endl);
      _xout2(" = " << sum_finp << endl);
      _xout2(" = " << sum_NO << endl);
      _xout1(" = " << sum_final << endl);
      
      // input 
      const std::vector<std::string> & finlines = finput.inlines();
      for ( unsigned int i = 0; i < finlines.size(); ++i ){
        MyOut::pcurout->buf << finlines[i] << endl;
        MyOut::pcurout->flushbuf();
      }
      // write to a file
      MyOut::pcurout->beq();
      // input-equation
      const std::vector<std::string> & fineq = finput.ineq();
      for ( unsigned int i = 0; i < fineq.size(); ++i ){
        MyOut::pcurout->buf << fineq[i] << endl;
        if ( i + 1 == fineq.size() ) {
          MyOut::pcurout->buf << "=";
          MyOut::pcurout->flushbuf();
          MyOut::pcurout->newlineeqn();
        }
      }
      MyOut::pcurout->buf <<sum_final << endl;
      MyOut::pcurout->eeq();
      if ( Input::iPars["prog"]["diagrams"] > 0 )
        Q2::printdiags(MyOut::pcurout ,sum_final);
      if ( Input::iPars["prog"]["algo"] > 0 ){
        ofstream falgout;
        std::string algofile;
        if ( outputfile.substr(outputfile.size()-4).compare(".tex") == 0 )
          algofile = outputfile.substr(0,outputfile.size()-4);
        else
          algofile = outputfile;
        algofile += ".alg";
        falgout.open(algofile.c_str());
        Q2::printalgo(falgout,sum_final);
      }
      finput.clear();
    }
  }
  // set current ouput back to default
  MyOut::pcurout = &MyOut::defout;
  fout.close();
  
//   // test nextwordpos
//   std::string
//     test("\\+\\cmd_ab_{cd}ef");
// 
//   lui ipos = 0;
//   for (uint ii = 0; ii < 10 && ipos < test.size(); ++ii ){
//     lui ipos1 = IL::nextwordpos(test,ipos,true,true);
//     xout << test.substr(ipos,ipos1-ipos) << "  ";
//     ipos = ipos1;
//   }
//   xout << std::endl;
//   
//   ipos = 0;
//   for (uint ii = 0; ii < 10 && ipos < test.size(); ++ii ){
//     lui ipos1 = IL::nextwordpos(test,ipos,false,true);
//     xout << test.substr(ipos,ipos1-ipos) << "  ";
//     ipos = ipos1;
//   }
//   xout << std::endl;
//   
//   ipos = 0;
//   for (uint ii = 0; ii < 10 && ipos < test.size(); ++ii ){
//     lui ipos1 = IL::nextwordpos(test,ipos,true,false);
//     xout << test.substr(ipos,ipos1-ipos) << "  ";
//     ipos = ipos1;
//   }
//   xout << std::endl;
//   
//   ipos = 0;
//   for (uint ii = 0; ii < 10 && ipos < test.size(); ++ii ){
//     lui ipos1 = IL::nextwordpos(test,ipos,false,false);
//     xout << test.substr(ipos,ipos1-ipos) << "  ";
//     ipos = ipos1;
//   }
//   xout << std::endl;
  
//  // test permutation multiplication
//  Orbital i("i"),j("j"),k("k");
//  Permut p1,p2;
//  p1 += Permut(j,k);
//  p1 += Permut(k,j);
//  p2 += Permut(i,j);
//  p2 += Permut(j,i);
//  xout << "P1, P2" << std::endl;
//  xout << p1 << p2 << std::endl;
//  p1 *= p2;
//  xout << p1 << std::endl;
  return 0;
}

