#if !defined MULTILIGN_H
#define MULTILIGN_H

#include <utility> //needed for library type pair
#include <vector>
#include <string>

using std::string;

class Multilign_object {
public:
  //!Default constructor:
  Multilign_object();
  //!Constructor:
  //!\param seq_list is a 2-D c string array storing the name of the input .seq files
  //!\param ct_list is a 2-D c string array storing the name of the output .ct files
  //!\param main_seq is a size_t number indicating which one in the seq_list(the 1st one by default) is the main sequence
  //!\param random is a boolean value indicating whether to randomize the sequence order or not (no randomization by default)
  //!\param dsv is a boolean value indicating whether to generate .dsv files or not (not by default)
  //!\param ali is a boolean value indicating whether to generate .ali files or not (not by default)
  Multilign_object(std::vector <string> seq_list, std::vector <string> ct_list, const bool random=false, const bool dsv=true, const bool ali=false);
  //Destructor
  ~Multilign_object();
  //functions
  //! Generate pairs of sequences for dynalign calculation
  //! \return the value of ErrorCode
  int generate_seq_pair();
  //! Generate pairs of the sequences and their corresponding cts
  //! \return the value of ErrorCode
  int generate_seq_ct();
  //! The core function doing dynalign calculation and templating
  //! In case of error, the function returns a non-zero that can be parsed by GetErrorMessage() or GetErrorMessageString().
  //! \param maxtrace is the maximum number of common structures to be determined.
  //! \param bpwin the the base pair window parameter, where 0 allows the structures to have similar pairs and larger windows make the structures more diverse. 
  //! \param awin is the alignment window parameter, where 0 allows the alignments to be similar and larger values make the alignments more diverse.
  //! \param percent is the maximum percent difference in total folding free energy change above the lowest for suboptimal common structures.
  //! \param imaxseparation is the maximum separation between aligned nucleotides.  Values >= 0 are the traditional parameter, those below zero trigger the HMM alignment method, which is now prefered.
  //! \param gap is the cost of adding gap nucleotides in the alignment in kcal/mol.
  //! \param singleinsert is whether single basepair inserts are allowed in one sequence vs the other.
  //! \param savefile is c-string with the name of a dynalign savefile (*.dsv) to be created.
  //! \param optimalonly can be used to turn on a calculation of only the energy (when true) and not the structures.
  //! \param singlefold_subopt_percent is the maximum % difference of folding energy above the lowest free energy structure for pairs in single sequence folding that will be allowed in the dynalign calculation.
  //! \param local is whether Dynalign is being run in local (true) or global mode (false).
  //! \param numProcessors is the number of processors to use for the calculation.  This requires a compilation for SMP.
  //! \param maxpairs is under development for multiple sequence folding.  Use -1 (default) for now.
  //! \Param maxdsvchange is the percentage above the lowest free energy for the base pair templating. Use 1 by default
  //! \param cycles is under development for multiple sequence folding.  Use 2 (default) for now.
  //! \return an int that indicates an error code (0 = no error, non-zero = error occurred).

  int calculate_dynalign(const short int maxtrace, 
			 const short int bpwin, const short int awin, 
			 const short int percent, const short int imaxseparation=-99, 
			 const float gap=0.4, 
			 const bool singleinsert=true, 
			 //const char savefile[]=NULL, 
			 const bool optimalonly=false, 
			 const short int singlefold_subopt_percent=30, 
			 const bool local=false, 
			 const short int numProcessors=1,
			 int maxpairs=-1,
			 const float maxdsvchange=1,
			 const short int cycles=2);


  //******************************************************
  //Functions to return error information
  //******************************************************
  
  //!Return an error code, where a return of zero is no error.
  
  //!  This function returns and error flag that is generated during construction by RNA(const char &filename, const int type, const bool IsRNA=true) or from CalculateFreeEnergy().
  //!    An error of zero is always no error.  Other codes are errors and a c-string can be fetched for the error with GetErrorMessage().
  //!\return An integer that provides the error code.
  int GetErrorCode(){return ErrorCode;};
  //!  Return error messages based on code from GetErrorCode and other error codes.    
  
  //!    0 = no error
  //!    1000 = Error associated with sequence 1 or with a procedure, function will get message from sequence 1 (the inherited RNA class).
  //!    2000 = Error associated with sequence 2, function will get message from sequence 2 (the RNA2 class).
  //!    3000 = Errors with each sequence, function will get messages from each.
  //!\param error is the integer error code provided by GetErrorCode().
  //!\return A pointer to a c string that provides an error message or from other functions that return integer error codes.
  string GetErrorMessage(int error);

  //This function is use to get the default maxpairs for users if maxpairs is set to -1
  //! \Param maxpairs is the parameters read from the configuration file
  //! \return the maxpairs: if maxpairs is set to negative, return the average length of all the input sequences; otherwise return the value in configuration file set by users
  int get_maxpairs(int maxpairs);
protected:
  int ErrorCode;
private:
  typedef std::pair<string, string>  pstr;
  typedef std::vector<string>::size_type vstr_index;
  std::vector <string> Seq_List; // the list of .seq filenames
  std::vector <string> Ct_List;
  pstr *seq_pair; // a list of pairs of .seq filename for pairwise dynalign calculation
  pstr *seq_ct; // a list of pairs of .seq and its corresponding .ct
  size_t Seq_Pair_Num; // the size of seq_pair array
  bool Dsv;// generate .dsv files or not (yes by default)
  bool Ali;// generate .ali files or not (not by default)
};

#endif
