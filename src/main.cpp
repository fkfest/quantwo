#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include "utilities.h"
#include "term.h"
#include "orbital.h"
#include "finput.h"
#include "globals.h"

using namespace std;

int main(int argc, char **argv) {
  // version of the program
  string version("1.0");
  // handle input and output
  string inputfile, outputfile,
    exePath = exepath(), outPath;
  if (argc>1) 
  {
    unsigned int i=1;
    while (argv[i][0]=='-')
    {// handle options
      if (strcmp(argv[i],"-h")==0 || strcmp(argv[i],"--help")==0)
      {
        cout << "quantwo <input-file> [<output-file>]" << endl;
        // print README file if exists
        ifstream readme;
        readme.open((exePath+"README").c_str());
        if (readme.is_open())
        {
          string line;
          while (readme.good())
          {
            getline (readme,line);
            cout << line << endl;
          }
        }
        return 0;
      }
      else if (strcmp(argv[i],"-v")==0 || strcmp(argv[i],"--version")==0)
      {
        cout << "Quantwo, Version " << version << endl;
        return 0;
      }
      ++i;
    }
    if (i<argc) {
      inputfile=argv[i];
      if (argc>i+1)
        outputfile=argv[i+1];
      else
        outputfile=inputfile+".tex";
    } else {
      error("Please provide an input file!");
    }
    outPath = DirName(outputfile);
  } else
    error("Please provide an input file!");
  // read input
  Finput finput;
  ifstream fin;
  fin.open(inputfile.c_str());
  if (fin.is_open())
  {
    string line;
    while (fin.good())
    {
      getline (fin,line);
      cout << line << endl;
      finput+=line;
    }
  }
  else
    error("Bad input file!");
  fin.close();
  
  if (finput.sumterms().size()==0)
    say("No equation in input file!");
  else
  {
    Sum<Term,double> sum_finp(finput.sumterms());
//    Sum<Term,double> sum_NO(Q2::normalOrderPH(sum_finp));
    Sum<Term,double> sum_NO(Q2::wick(sum_finp));
    Sum<Term,double> sum_final(Q2::reduceSum(sum_NO));
    cout << finput << endl;
//    cout << " = " << sum_finp << endl;
//    cout << " = " << sum_NO << endl;
    cout << " = " << sum_final << endl;
    // write to a file
    ofstream fout;
    fout.open(outputfile.c_str());
    Output myfout(fout);
    // set current ouput to fout
    MyOut::pcurout = &myfout;
    // "save" 15 lines for titel
    MyOut::pcurout->nlines=15;
    MyOut::pcurout->beq();
    MyOut::pcurout->buf <<sum_final << endl;
    MyOut::pcurout->eeq();
    // set current ouput back to default
    MyOut::pcurout = &MyOut::defout;
    fout.close();
  }
  return 0;
}

