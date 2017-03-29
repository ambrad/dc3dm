/* dc3dm: Software to form and apply a 3D DDM operator for a nonuniformly
 *        discretized rectangular fault
 *   Version 0.3
 *   Andrew M. Bradley
 *   ambrad@cs.stanford.edu
 *   CDFM Group, Geophysics, Stanford
 * https://pangea.stanford.edu/research/CDFM/software
 * dc3dm is licensed as follows:
 *   Open Source Initiative OSI - Eclipse Public License 1.0
 *   http://www.opensource.org/licenses/eclipse-1.0
 */
 
// Key-value-file interface to dc3dm's functionality.

#ifndef INCLUDE_DC3DM_DC3DM
#define INCLUDE_DC3DM_DC3DM

#include <string>
#include "util/include/KeyValueFile.hpp"
#include "dc3dmHelp.hpp"

namespace dc3dm {
struct Inputs {
  bool GetInputs (const std::string& kvf_fn);
  bool GetInputs (const util::KeyValueFile* kvf);
  virtual bool ProcessKvf(const util::KeyValueFile* kvf, std::string& err) = 0;
};

namespace mesh {
void PrintHelp();
int Run(const std::string& kvf_fn);
int Run(const util::KeyValueFile* kvf);
}
namespace build {
void PrintHelp();
int Run(const std::string& kvf_fn);
int Run(const util::KeyValueFile* kvf);
}
namespace compress {
void PrintHelp();
int Run(const std::string& kvf_fn);
int Run(const util::KeyValueFile* kvf);
}
namespace compresskxk {
void PrintHelp();
int Run(const std::string& kvf_fn);
int Run(const util::KeyValueFile* kvf);
}
}

#endif
