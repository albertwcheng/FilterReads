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
	unsigned int line1;
	unsigned int line2;
	unsigned int line3;
	unsigned int line4;
	inline bool operator < (const FastqRecord& fqr) const{
		return line1<fqr.line1;
	}
};

bool readFastqIntoNameKeyedRecord(string filename,map<string, FastqRecord>& records){
	
	ifstream fil(filename.c_str());
	unsigned int lino=0; //store the number of line fastq content
	unsigned int rlino=0; //store the real number of line
	
	
	string name;
	//unsigned int firstlineseek;
	vector<string> splits;
	FastqRecord fqr;
	
	while(fil.good()){
	
		
		string buffer;
		unsigned int cseek=fil.tellg();
		
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
				
				
				StringUtil::split(buffer,"/",splits);
				name=splits[0];
				fqr.line1=cseek;
				//cerr<<"name="<<name<<endl;
				
				break;
			case 2:
				fqr.line2=cseek;
				break;
			case 3:
				if(buffer[0]!='+'){
					cerr<<"fastq formatting error! @line "<<rlino<<" : "<<buffer<<endl;
					cerr<<"line 3 of each fastq record should start with + sign"<<endl;
					fil.close();
					return false;
				}
				fqr.line3=cseek;
				break;
			case 0:
				fqr.line4=cseek;
				records.insert(map<string, FastqRecord>::value_type(name,fqr));
				break;
		}
		
	}
	
	fil.close();
	return true;

}

void copyLines(string fileSrc,string fileDst,set<FastqRecord>&linos){
	
	
	ifstream fil(fileSrc.c_str());
	ofstream fout(fileDst.c_str());
	
	for(set<FastqRecord>::iterator linoI=linos.begin();linoI!=linos.end();linoI++){
		
		
		
		string buffer;
		
		fil.seekg(linoI->line1);
		getline(fil,buffer);
		fout<<buffer<<endl;
		
		fil.seekg(linoI->line2);
		getline(fil,buffer);
		fout<<buffer<<endl;	
		
		fil.seekg(linoI->line3);
		getline(fil,buffer);
		fout<<buffer<<endl;
		
		fil.seekg(linoI->line4);
		getline(fil,buffer);
		fout<<buffer<<endl;
			
	}
	
	
	
	fil.close();
	fout.close();

	
}

int runFilterReads( filterReads_opts &opts){

	map<string, FastqRecord> file1r; //name -> seek
	map<string, FastqRecord> file2r; //name -> seek
	
	unsigned int readsOut1=0;
	unsigned int readsOut2=0;
	
	cerr<<"First Pass"<<endl;
	cerr<<"Reading file1 "<<opts.file1<<endl;
	if(!readFastqIntoNameKeyedRecord(opts.file1,file1r)){
		cerr<<"Failed. Abort"<<endl;
		return 1;
	}
	cerr<<file1r.size()<<" reads read"<<endl;
	cerr<<"Reading file2 "<<opts.file2<<endl;
	if(!readFastqIntoNameKeyedRecord(opts.file2,file2r)){
		cerr<<"Failed. Abort"<<endl;
		return 1;
	}
	cerr<<file2r.size()<<" reads read"<<endl;
	cerr<<"Finished first-passing files. Start filtering"<<endl;
	
	set<FastqRecord> f1linesToOutput;
	set<FastqRecord> f2linesToOutput;
	
	map<string, FastqRecord>::iterator i;
	map<string, FastqRecord>::iterator j;
	
	
	if(opts.pairedReads){
		cerr<<"filter paired"<<endl;
		for(i=file1r.begin();i!=file1r.end();i++){
			if((j=file2r.find(i->first))!=file2r.end()){
				f1linesToOutput.insert(i->second);
				f2linesToOutput.insert(j->second);
				readsOut1++;
				readsOut2++;
			}
		}
	}
	
	cerr<<"Second Pass"<<endl;
	
	cerr<<"Copying filtered lines from file1 to outfile1"<<endl;
	copyLines(opts.file1,opts.outfile1,f1linesToOutput);
	cerr<<readsOut1<<" reads output to "<<opts.outfile1<<endl;
	
	cerr<<"Copying filtered lines from file2 to outfile2"<<endl;
	copyLines(opts.file2,opts.outfile2,f2linesToOutput);
	cerr<<readsOut2<<" reads output to "<<opts.outfile2<<endl;
	
	cerr<<"<Done>"<<endl;
	

	
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
	
	if(!opts.pairedReads){
		cerr<<"currently we only do filtering for paired reads. please use --paired option"<<endl;
		printUsage(argsFinal.programName);
		return 1;
	}
	
	return runFilterReads(opts);
	
	return 0;

}