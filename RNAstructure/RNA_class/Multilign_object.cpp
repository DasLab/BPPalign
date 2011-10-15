#include "Multilign_object.hpp"
#include <sstream>
#include "../src/random.h" //random number generator
#include <time.h>
#include "Dynalign_object.h"

using std::string;

//constructor
Multilign_object::Multilign_object(std::vector <string> seq_list, std::vector <string> ct_list, const bool random, const bool dsv, const bool ali):\
  Seq_List(seq_list), Ct_List(ct_list), Dsv(dsv), Ali(ali){
  if( ErrorCode = generate_seq_ct() )
    return;
  size_t Seq_Num = Seq_List.size();
  if( random ){
    randomnumber rand;//using DHM's random number generator
    rand.seed( time( NULL ) );
    for (size_t i = 1; i < Seq_Num; ++i){
      size_t r = i + size_t ( rand.roll() * Seq_Num ) % ( Seq_Num - i );
      pstr tmp = *(seq_ct + i); *(seq_ct + i) = *(seq_ct + r); *(seq_ct + r) = tmp; // randomize the order except the first pair of seq_ct, i.e. the main pair 
    }
  }

  if( ErrorCode = generate_seq_pair() )
    return;
  ErrorCode=0;
}


//Destructor
Multilign_object::~Multilign_object(){
  delete [] seq_pair;
  delete [] seq_ct;
}



// initialize data members Seq_List and Ct_List
int Multilign_object::generate_seq_ct(){
  if(Seq_List.size()!=Ct_List.size()){
    return 5001; //if not equal, then the program should be exited
  }
  
  seq_ct = new pstr [Seq_List.size()];
  
  for(vstr_index i = 0; i != Seq_List.size(); ++i){
    *(seq_ct + i) = std::make_pair(Seq_List[i], Ct_List[i]); // add the pair into seq_ct array
    //std::cout << (seq_ct+j)->second << " " << j <<" 44l" << std::endl;
  }
  return 0;
}



// initialize data member seq_pair
int Multilign_object::generate_seq_pair(){
  Seq_Pair_Num = Seq_List.size() - 1;
  if (Seq_Pair_Num < 1 ){
    return 5002; //if not enough pairs of calculation, exit
  }
  seq_pair = new pstr [Seq_Pair_Num];
  for(vstr_index i = 0; i != Seq_Pair_Num; ++i){
    *(seq_pair + i) = std::make_pair(seq_ct->first, (seq_ct+i+1)->first);
    //std::cout << (seq_pair+i)->first << " " << (seq_pair+i)->second << std::endl;
  }
  return 0;
}



//return error messages based on the code from GetErrorCode and other error codes
string Multilign_object::GetErrorMessage(int error){
  if (error==0) return "No Error.\n";
  else if (error==5001) return "Different number of seq filenames in ct filenames in the conf file.\n";
  else if (error==5002) return "Less than 2 seq filenames in conf file.\n";
  //else if (error<5000) return Dynalign_object::GetErrorMessage(error);
  else if (error < 5000) return "Error occurs in Dynalign_object class.\n";
  else return "Unknown Error\n";
}



// the core function doing dynalign calculation and templating
int Multilign_object::calculate_dynalign(const short int maxtrace, 
					 const short int bpwin, const short int awin, 
					 const short int percent, 
					 const short int imaxseparation, 
					 const float gap, 
					 const bool singleinsert, 
					 const bool optimalonly, 
					 const short int singlefold_subopt_percent, 
					 const bool local, 
					 const short int numProcessors,
					 int maxpairs,
					 const float maxdsvchange,
					 const short int cycles){
  maxpairs = get_maxpairs(maxpairs);
  std::stringstream ss;
  std::string **dsv_file = NULL, **ali_file = NULL;
  // if Dsv is true, .dsv files will be stored as j.i_<seq1_name>_<seq2_name>.aout, where j is the cycle number and i is the order of calculation in the cycle
  if (Dsv){
    dsv_file= new std::string * [cycles];
    for ( int j =0; j< cycles; ++j){
      dsv_file[j] = new std::string [Seq_Pair_Num];
      for ( size_t i = 0; i < Seq_Pair_Num; ++i){
        std::string one, two;
	//std::cout << seq_pair[i].first << " " << seq_pair[i].second << std::endl;
        size_t index = std::string::npos;	
        one = seq_pair[i].first;
        if ( ( index = one.rfind(".seq") ) == (one.length() - 4) )
          one=one.substr(0, index);
	
        two = seq_pair[i].second;
        if ( ( index = two.rfind(".seq") ) == (two.length() -4) )
          two = two.substr(0, index);
	
        index = one.rfind('/');
        while(index!=std::string::npos && one.substr(index-1, 2)=="\\/")
          index=one.rfind('/', index);
        if(index!=std::string::npos)
          dsv_file[j][i] = one.substr(index+1);
	else dsv_file[j][i] = one;
	
        index = two.rfind('/');
        while(index!= std::string::npos && two.substr(index-1, 2)=="\\/")
          index = two.rfind('/', index);
        if(index!=std::string::npos)
          dsv_file[j][i] = dsv_file[j][i] + "_" + two.substr(index+1);	
	else dsv_file[j][i] = dsv_file[j][i] + "_" + two;
	
        ss.str(""); ss << i+1; 
        dsv_file[j][i] = ss.str() + "_" + dsv_file[j][i] ;
        ss.str(""); ss << j+1;
        dsv_file[j][i] = ss.str() + "." + dsv_file[j][i]+".dsv";
	//cout << dsv_file[j][i] << std::endl;
      }
    }
  }
  // if Ali is true, .aout files will be stored as j.i_<seq1_name>_<seq2_name>.aout, where j is the cycle number and i is the order of calculation in the cycle
  if (Ali){
    ali_file= new std::string * [cycles];
    for ( int j =0; j< cycles; ++j){
      ali_file[j] = new std::string [Seq_Pair_Num];
      for ( size_t i = 0; i < Seq_Pair_Num; ++i){
        std::string one, two;
        size_t index = std::string::npos;
        one = seq_pair[i].first;
        if ( ( index = one.rfind(".seq") ) == (one.length() - 4) )
          one=one.substr(0, index);
	
        two = seq_pair[i].second;
        if ( ( index = two.rfind(".seq") ) == (two.length() -4) )
          two = two.substr(0, index);

        index = one.rfind('/');
        while(index!=std::string::npos && one.substr(index-1, 2)=="\\/")
          index=one.rfind('/', index);
        if(index!=std::string::npos)
          ali_file[j][i] = one.substr(index+1);
	else ali_file[j][i] = one;
	
        index = two.rfind('/');
        while(index!= std::string::npos && two.substr(index-1, 2)=="\\/")
          index = two.rfind('/', index);
        if(index!=std::string::npos)
          ali_file[j][i] = ali_file[j][i] + "_" + two.substr(index+1);
	else ali_file[j][i] = ali_file[j][i] + "_" + two;
	
        ss.str(""); ss << i+1;
        ali_file[j][i] = ss.str() + "_" + ali_file[j][i] ;
        ss.str(""); ss << j+1;
        ali_file[j][i] = ss.str() + "." + ali_file[j][i]+".aout";
	//cout << ali_file[j][i] << std::endl;
      }
    }
  } 
  
  int bp_num = 0;

  for ( int j =0 ; j < cycles; ++j){
    for ( size_t i = 0; i < Seq_Pair_Num; ++i){
      Dynalign_object *dynalign_instance = new Dynalign_object((seq_pair+i)->first.c_str(), 2, (seq_pair+i)->second.c_str(), 2);
      //cout << (seq_pair+i)->first << " " << (seq_pair+i)->second << " " << i << " " << j << " " << cycles <<endl;
      if (i==0 && j==0) {
	// doing the core dynalign calculation
	if (ErrorCode = dynalign_instance->Dynalign( maxtrace, bpwin, awin, 
						     percent, imaxseparation, 
						     gap, singleinsert, 
						     dsv_file[j][i].c_str(), 
						     optimalonly, 
						     singlefold_subopt_percent,
						     local, numProcessors,
						     maxpairs, cycles ) ){
	  //std::cout << "ERROR(cycle " << ++j << ' ' << (seq_pair+i)->first << ' ' << (seq_pair+i)->second << "): " << dynalign_instance->GetErrorMessage(ErrorCode);
	  return ErrorCode;
	}
	else{
	  // templating the dsv file
	  Dynalign_object *dynalign_tmp = new Dynalign_object(dsv_file[j][i].c_str());
	  // set a cutoff to count the number of base pairs with lowest free energy below this cutoff
	  int cutoff = int(dynalign_tmp->GetLowestEnergy()*0.8);
	  //std::cout << cutoff << " 198l" << endl;
	  // lenght is the length of main sequences
	  int length = RNA( (seq_pair+i)->first.c_str(), 2 ).GetSequenceLength();
	  //cout << length << " 200l" << endl;
	  // count the number of base pairs with lowest free energy below this cutoff
	  for (int i=1; i<=length; ++i){
	    for (int j=i; j<=length; ++j){
	      if (dynalign_tmp->GetBestPairEnergy(1, i, j) < cutoff )
		++bp_num;
	    }
	  }
	  //cout << bp_num;
	  bp_num -= length;
	  if (Seq_Pair_Num > 1)
	    // now the bp_num is the number gradually disallowed for each dynalign calculation in 1st cycle
	    bp_num /= (Seq_Pair_Num-1);
	  // cout << bp_num << " zane" << endl;
	  delete dynalign_tmp;
	}
      }
  
      else{
        if (i==0 && j!=0) {
	  if (ErrorCode = dynalign_instance->Templatefromdsv(dsv_file[j-1][Seq_Pair_Num-1].c_str(), maxdsvchange)){
	    //std::cout << "ERROR(cycle " << ++j << ' '<< (seq_pair+i)->first << ' ' << (seq_pair+i)->second << "): "  <<  dynalign_instance->GetErrorMessage(ErrorCode);
	    return ErrorCode;
	  }
	}
        else {
	  if (ErrorCode = dynalign_instance->Templatefromdsv(dsv_file[j][i-1].c_str(), maxdsvchange)){
	    //std::cout <<"ERROR(cycle " << ++j << ' ' << (seq_pair+i)->first << ' ' << (seq_pair+i)->second << "): "  <<  dynalign_instance->GetErrorMessage(ErrorCode);
	    return ErrorCode;
	  }
	}
	if (j==0){
	  if (ErrorCode = dynalign_instance->Dynalign( maxtrace, bpwin, awin,
						       percent, imaxseparation,
						       gap, singleinsert,
						       dsv_file[j][i].c_str(),
						       optimalonly,
						       singlefold_subopt_percent,
						       local, numProcessors,
						       maxpairs+bp_num*(Seq_Pair_Num-i), cycles ) ){
	    //std::cout << "ERROR(cycle " << ++j << ' ' << (seq_pair+i)->first << ' ' << (seq_pair+i)->second << "): "  << dynalign_instance->GetErrorMessage(ErrorCode);
	    return ErrorCode;
	  }
	}
	else{
	  if (ErrorCode = dynalign_instance->Dynalign( maxtrace, bpwin, awin,
						       percent, imaxseparation,
						       gap, singleinsert,
						       dsv_file[j][i].c_str(),
						       optimalonly,
						       singlefold_subopt_percent,
						       local, numProcessors,
						       maxpairs, cycles)){
	    //std::cout << "ERROR(cycle " << ++j << ' ' << (seq_pair+i)->first << ' ' << (seq_pair+i)->second << "): "  << dynalign_instance->GetErrorMessage(ErrorCode);
	    return ErrorCode;
	  }
	}
      }
      dynalign_instance->WriteAlignment( ali_file[j][i].c_str() );
      if (j==cycles-1){ // this is the last cycle of Dynalign calculation, output ct files
	//std::cout << (seq_ct+i+1)->second.c_str() << std::endl;
        if (i==Seq_Pair_Num-1){
          if (ErrorCode = dynalign_instance->GetRNA1()->WriteCt( seq_ct->second.c_str() )) return ErrorCode;
	  if (ErrorCode = dynalign_instance->GetRNA2()->WriteCt( (seq_ct+i+1)->second.c_str() )) return ErrorCode;
	}
	else 
	  if (ErrorCode = dynalign_instance->GetRNA2()->WriteCt( (seq_ct+i+1)->second.c_str() )) return ErrorCode;
      }
      delete dynalign_instance;
    }
  }
  

  for (int i = 0; i < cycles; ++i){
    if (Dsv)
      delete  [] dsv_file[i];
    if (Ali)
      delete [] ali_file[i];
  }
  if (Dsv) delete [] dsv_file;
  if (Ali) delete [] ali_file;
}


int Multilign_object::get_maxpairs(int maxpairs) {
  if(maxpairs >= 0 ) return maxpairs;
  else {
    int nt=0;
    for ( size_t i = 0; i < Seq_List.size(); ++i){
      RNA *seq = new RNA(Seq_List[i].c_str(), 2);
      nt = nt + seq->GetSequenceLength();
    }
    return nt/Seq_List.size(); //return the average length of all the input sequences
  }
}
