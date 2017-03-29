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
 
#include <stdio.h>
#include "dc3dm.hpp"
using namespace std;
using namespace util;

struct CommandInput : public dc3dm::Inputs {
  string command;

  virtual bool ProcessKvf (const KeyValueFile* kvf, string& err) {
    const string* s;
    if (!kvf->GetString("command", s)) {
      err = "Missing <command>.";
      return false;
    } else {
      command = *s;
      return true;
    }
  }
};

int main (int argc, char** argv) {
  if (argc == 1 || (argc == 2 && string(argv[1]) == "help")) {
    printf(_dc3dm_all_PrintHelp_text_);
    printf(_dc3dm_Header_text_);
    return 0;
  } else if (argc == 3 && string(argv[1]) == "help") {
    string cmd = argv[2];
    if (cmd == "mesh") {
      dc3dm::mesh::PrintHelp();
      return 0;
    } else if (cmd == "build") {
      dc3dm::build::PrintHelp();
      return 0;
    } else if (cmd == "compress") {
      dc3dm::compress::PrintHelp();
      return 0;
    } else if (cmd == "compresskxk") {
      dc3dm::compresskxk::PrintHelp();
      return 0;
    } else {
      printf("%s is not a valid command.\n\n", cmd.c_str());
      printf(_dc3dm_all_PrintHelp_text_);
      printf(_dc3dm_Header_text_);
      return -1;
    }
  } else if (argc == 2) {
    CommandInput in;
    if (!in.GetInputs(argv[1])) return -1;
    if (in.command == "mesh")
      return dc3dm::mesh::Run(argv[1]);
    else if (in.command == "build")
      return dc3dm::build::Run(argv[1]);
    else if (in.command == "compress")
      return dc3dm::compress::Run(argv[1]);
    else if (in.command == "compresskxk")
      return dc3dm::compresskxk::Run(argv[1]);
    else {
      printf("%s is not a valid command.\n\n", in.command.c_str());
      printf(_dc3dm_all_PrintHelp_text_);
      printf(_dc3dm_Header_text_);
      return -1;
    }
  }
}
