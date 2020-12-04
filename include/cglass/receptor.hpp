#ifndef _CGLASS_RECEPTOR_H_
#define _CGLASS_RECEPTOR_H_

#include "site.hpp"
#include "cortex.hpp"

/* receptor.hpp declares the Receptor class. The Receptor class represents
 * biological receptors that can attach to crosslinking and motor proteins. 
 * In this code, Receptor's are modelled as sites on a mesh that can attach 
 * to the Anchor of a Crosslinker.
*/

class Receptor: public Site {
private:
  receptor_parameters *sparams_; // Input parameters unique to receptors
public:
  Receptor(unsigned long seed); // Constructor with RNG seed input

  // Save yaml input file parameters locally
  void Init(receptor_parameters *sparams);

  // Initialize Receptor members with yaml input file parameters
  void SetParameters();

  // Read/write binaries
  void WriteSpec(std::fstream &ospec);
  void ReadSpec(std::fstream &ispec);

  // Convert binary data to text. Static to avoid needing to istantiate
  // species members.
  static void ConvertSpec(std::fstream &ispec, std::fstream &otext);
  static void WriteSpecTextHeader(std::fstream &otext);
};

#endif
