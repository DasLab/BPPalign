

//Use the precompiler to make sure the class definition is not included more than once.
#if !defined(TWORNA_H)
#define TWORNA_H

#include "../RNA_class/RNA.h"
#include "thermodynamics.h"

//Use a precompiler definition for the maximum length of an error message (used in GetErrorMessage, below).
#define COMPOUNDMESSAGELENGTH 200

//! TwoRNA Class.
/*!
	The TwoRNA class provides an entry point for all the two sequence prediction routines of RNAstructure.  
	This contains two instances of the RNA class to provide the functionality of RNA.
*/

//Note the stylized comments provide facility for automatic documentation via doxygen.

class TwoRNA  {

	public:
		//******************************************************
		//Constructors
		//******************************************************

		//! Constructor - user provides a sequences as c strings.

		//!	Input sequences should contain A,C,G,T,U,a,c,g,t,u,x,X.
		//!	Capitalization makes no difference.
		//!	T=t=u=U.  If IsRNA is true, the backbone is RNA, so U is assumed.  If IsRNA is false, the backbone is DNA, so T is assumed.
		//!	x=X= nucleotide that neither stacks nor pairs.
		//!	For now, any unknown nuc is considered 'X'.
		//! Both sequences are passed to underlying RNA classes for each sequence.
		//!	\param sequence1 is a NULL terminated c string for sequence 1.
		//!	\param sequence2 is a NULL terminated c string for sequence 2.
		//!	\param IsRNA is a bool that indicates whether these sequences are RNA or DNA.  true=RNA.  false=DNA.  Default is true.  Both sequences must have the same backbone.
		TwoRNA(const char sequence1[], const char sequence2[], bool IsRNA=true); 


		//! Constructor - user provides a filenames for existing files as a c string.

		//!	The existing files, specified by filenames, can either be a ct file, a sequence, or an RNAstructure save file. 
		//!	Therefore, the user provides a flag for the file type: 
		//!		type = 1 => .ct file, type = 2 => .seq file, type = 3 => partition function save (.pfs) file, type = 4 => folding save file (.sav).
		//! The file opening is performed by the constructors for the RNA classes that underlie each sequence.
		//!	This constructor generates internal error codes that can be accessed by GetErrorCode() after the constructor is called.  0 = no error.
		//! The errorcode can be resolved to a c string using GetErrorMessage.		
		//!	Note that the contructor needs to be explicitly told, via IsRNA, what the backbone is because files do not store this information.
		//! Note also that save files explicitly store the thermodynamic parameters, therefore changing the backbone type as compaared to the original calculation will not change structure predictions.
		//! \param filename1 is a null terminated c string and refers to sequence 1.
		//! \param filename2 is a null terminated c string and refers to sequence 2.
		//! \param type1 is an integer that indicates the file type for sequence 1.
		//! \param type2 is an integer that indicates the file type for sequence 2.
		//!	\param IsRNA is a bool that indicates whether these sequences are RNA or DNA.  true=RNA.  false=DNA.  Default is true.  Only one backbone is allowed for both sequences.
		TwoRNA(const char filename1[], const int type1, const char filename2[], const int type2, bool IsRNA=true);

		//! Constructor
		//! Default constructor that requires no parameters.
		TwoRNA();

		

		//******************************************************
		//Functions to return error information
		//******************************************************

		//! Return an error code, where a return of zero is no error.

		//!	This function returns and error flag that is generated during construction by RNA(const char &filename, const int type, const bool IsRNA=true) or from CalculateFreeEnergy().
		//!		An error of zero is always no error.  Other codes are errors and a c-string can be fetched for the error with GetErrorMessage().
		//! \return An integer that provides the error code.
		int GetErrorCode();

		//!	Return error messages based on code from GetErrorCode and other error codes.		

		//!		0 = no error
		//!		1000 = Error associated with sequence 1 or with a procedure, function will get message from sequence 1 (the inherited RNA class).
		//!		2000 = Error associated with sequence 2, function will get message from sequence 2 (the RNA2 class).
		//!		3000 = Errors with each sequence, function will get messages from each.
		//! \param error is the integer error code provided by GetErrorCode().
		//! \return A pointer to a c string that provides an error message or from other functions that return integer error codes.
		char* GetErrorMessage(const int error);

		//!	Return error messages based on code from GetErrorCode and other error codes.		

		//!	Although RNA generally uses c strings, this member function returns a string that is suitable for interfacing with JAVA, etc.
		//!		See the error list in the GetErrorMessage() entry.
		//!\param error is the integer error code provided by GetErrorCode() or from other functions that return integer error codes.
		//!\return A string that provides an error message.
		std::string GetErrorMessageString(const int error);

		//******************************************************************
		//Functions that access underlying classes:
		//******************************************************************
		
		//!Access the underlying RNA class.
		//!This is provided for use with two sequence methods.
		//!Generally, there is no need for end users to use this function.

		//!\return A pointer to the underlying RNA class for sequence 1.
		RNA *GetRNA1();

		//!Access the underlying RNA class.
		//!This is provided for use with two sequence methods.
		//!Generally, there is no need for end users to use this function.

		//!\return A pointer to the underlying RNA class for sequence 2.
		RNA *GetRNA2();


		//******************************************************
		//Destructor
		//******************************************************

		//!Destructor

		//!The destructor cleans up all memory allocation.
		~TwoRNA();

		//! Provide facility for storing and accessing error messages.
		char compoundmessage[COMPOUNDMESSAGELENGTH];

	protected:
		//Integer to keep track of error codes.
		//These errors result on file i/o problems during construction.
		//The errors can be accessed using GetErrorCode().
		int ErrorCodeTwo;

	private:
		
		//Use two RNA classes because this encapsulates two sequences.
		RNA *rna1,*rna2;

		


};

#endif //TWORNA_H defined