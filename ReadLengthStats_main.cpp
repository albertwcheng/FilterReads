/***************************************************************************
 Copyright 2011 Wu Albert Cheng <albertwcheng@gmail.com>
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 *******************************************************************************/ 

//Filter paired reads, and or read quality stuff etc
#include <string>
#include <iostream>
#include <fstream>
#include <set>	
#include <StringUtil.h>
#include "AdvGetOptCpp/AdvGetOpt.h"
using namespace std;

class ReadLengthStats_opts
{
  public:
    vector<string> infiles;
	string ofile;

};


void printUsage(string programName){
	cerr<<programName<<" [--out ofile ] file1 file2 ... fileN"<<endl;
	
}

void addStat(vector<int>& readLengthStats,int readLength){
	while(readLengthStats.size()<readLength){
		readLengthStats.push_back(0);	
	}
	
	readLengthStats[readLength-1]++;
}

bool statFastqReadLengths(string filename,vector<int>& readLengthStats,int& totalNumReads){
	
	ifstream fil(filename.c_str());
	unsigned int lino=0; //store the number of line fastq content
	unsigned int rlino=0; //store the real number of line

	
	int length;
	
	while(fil.good()){
	
		
		string buffer;
		
		
		getline(fil,buffer);
		
		rlino++;
		
		if(buffer.length()<1)
			continue;
		
		lino++;
		switch(lino%4){
			case 1:
				if(buffer[0]!='@'){
					cerr<<"fastq formating error! @line "<<rlino<<" : "<<buffer<<endl;
					cerr<<"line 1 of each fastq record should start with @ sign"<<endl;
					fil.close();
					return false;
				}
				
				break;
			case 2:

				length=buffer.length();
				
				totalNumReads++;
				
				addStat(readLengthStats,length);
				
				break;
			case 3:
				if(buffer[0]!='+'){
					cerr<<"fastq formatting error! @line "<<rlino<<" : "<<buffer<<endl;
					cerr<<"line 3 of each fastq record should start with + sign"<<endl;
					fil.close();
					return false;
				}

				break;
			case 0:

				break;
		}
		
	}
	
	fil.close();
	return true;

}



int runReadLengthStats( ReadLengthStats_opts &opts){
	vector<int> readLengthStats;
	int totalNumReads=0;
	
	for(vector<string>::iterator i=opts.infiles.begin();i!=opts.infiles.end();i++)
		statFastqReadLengths(*i,readLengthStats,totalNumReads);
	
	ofstream ofil(opts.ofile.c_str());
	
	ofil<<"X\t#Len=X\t\%Len=X\t#Len<=X\t\%Len<=X\t#Len>=X\t\%Len>=X"<<endl;
	
	int cumNumX=0;
	
	for(int i=0;i<readLengthStats.size();i++){
		int curLength=i+1;
		int numX=readLengthStats[i];
		
		cumNumX+=numX;
		int invCumNumX=totalNumReads-cumNumX+numX; 
		double perNumX=double(numX)/totalNumReads*100;
		double perCumX=double(cumNumX)/totalNumReads*100;
		double perInvCumNumX=double(invCumNumX)/totalNumReads*100;
		
		ofil<<curLength<<"\t";
		ofil<<numX<<"\t";
		ofil<<perNumX<<"\t";
		ofil<<cumNumX<<"\t";
		ofil<<perCumX<<"\t";
		ofil<<invCumNumX<<"\t";
		ofil<<perInvCumNumX<<endl;
	}
	
	ofil.close();
	
	return 0;
}


int main(int argc,char*argv[])
{
	
	vector<string> long_options;

	long_options.push_back("ofile=");
	
	
	ReadLengthStats_opts opts;
	
	map<string,string> optmap;
	
	EasyAdvGetOptOut argsFinal=easyAdvGetOpt(argc,argv,"",&long_options);
	if(argsFinal.success){
		argsFinal.print(cerr);
		parseOptsIntoMap(argsFinal.opts,optmap);
		
	}
	else{
		printUsage(argsFinal.programName);
		return 1;
	}
	
	if(argsFinal.args.size()<1){
		printUsage(argsFinal.programName);
		return 1;
	
	}
	
	
	opts.ofile=getOptValue(optmap,"--ofile","/dev/stdout");
	
	for(vector<string>::iterator i=argsFinal.args.begin();i!=argsFinal.args.end();i++){
		opts.infiles.push_back(*i);	
	}
	
	
	return runReadLengthStats(opts);
	
	return 0;

}