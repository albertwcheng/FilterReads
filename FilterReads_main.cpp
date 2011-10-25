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
#include <StringUtil.h>
#include "AdvGetOptCpp/AdvGetOpt.h"
using namespace std;

class filterReads_opts
{
  public:
    bool pairedReads;
    string file1;
    string file2;
    string outfile1;
    string outfile2;
  
  filterReads_opts():pairedReads(false)
  {
  
  }

};


void printUsage(string programName){
	cerr<<programName<<" [options] file1 file2 outfile1 outfile2"<<endl;
	
	cerr<<"Options:"<<endl;
	cerr<<"--paired. retain only paired reads"<<endl;
	
	
}

class FastqRecord{
	public:
	string line1;
	string line2;
	string line3;
	string line4;
};

bool readFastqIntoNameKeyedRecord(string filename,map<string, FastqRecord>& records){
	
	ifstream fil(filename.c_str());
	unsigned int lino=0;
	
	FastqRecord fqr;
	string name;
	vector<string> splits;
	
	while(fil.good()){
	
		string buffer;
		getline(fil,buffer);
		
		if(buffer.length()<1)
			continue;
		
		lino++;
		switch(lino%4){
			case 1:
				fqr.line1=buffer;
				
				StringUtil::split(buffer,"/",splits);
				name=splits[0];
				//cerr<<"name="<<name<<endl;
				
				break;
			case 2:
				fqr.line2=buffer;
				break;
			case 3:
				fqr.line3=buffer;
				break;
			case 0:
				fqr.line4=buffer;
				records.insert(map<string, FastqRecord>::value_type(name,fqr));
				break;
		}
		
	}
	
	fil.close();
	return true;

}

void outputFastqRecord(ostream &os, FastqRecord &fqr){
	os<<fqr.line1<<endl;
	os<<fqr.line2<<endl;
	os<<fqr.line3<<endl;
	os<<fqr.line4<<endl;
}

int runFilterReads( filterReads_opts &opts){

	map<string, FastqRecord> file1r;
	map<string, FastqRecord> file2r;
	
	unsigned int readsOut1=0;
	unsigned int readsOut2=0;
	
	
	cerr<<"Reading file1 "<<opts.file1<<endl;
	readFastqIntoNameKeyedRecord(opts.file1,file1r);
	cerr<<file1r.size()<<" reads read"<<endl;
	cerr<<"Reading file2 "<<opts.file2<<endl;
	readFastqIntoNameKeyedRecord(opts.file2,file2r);
	cerr<<file2r.size()<<" reads read"<<endl;
	cerr<<"Finished reading files. Start filtering"<<endl;
	
	ofstream of1(opts.outfile1.c_str());
	ofstream of2(opts.outfile2.c_str());
	map<string, FastqRecord>::iterator i;
	
	map<string, FastqRecord>::iterator j;
	
	if(opts.pairedReads){
		cerr<<"filter paired"<<endl;
		for(i=file1r.begin();i!=file1r.end();i++){
			if((j=file2r.find(i->first))!=file2r.end()){
				outputFastqRecord(of1,i->second);
				outputFastqRecord(of2,j->second);
				readsOut1++;
				readsOut2++;
			}
		}
	}
	
	cerr<<readsOut1<<" reads output to "<<opts.outfile1<<endl;
	
	cerr<<readsOut2<<" reads output to "<<opts.outfile2<<endl;
	
	of1.close();
	of2.close();
	
	return 0;
}


int main(int argc,char*argv[])
{
	
	vector<string> long_options;
	//required:
	long_options.push_back("paired");

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
	
	filterReads_opts opts;
	
	if(argsFinal.args.size()<4)
	{
		printUsage(argsFinal.programName);
		return 1;
	}
	
	opts.file1=argsFinal.args[0];
	opts.file2=argsFinal.args[1];
	opts.outfile1=argsFinal.args[2];
	opts.outfile2=argsFinal.args[3];
	
	
	//clustergramOpts.bamfile=getOptValue(optmap,"--target-read-file");
	
	opts.pairedReads=hasOpt(optmap,"--paired");
	
	return runFilterReads(opts);
	
	return 0;

}