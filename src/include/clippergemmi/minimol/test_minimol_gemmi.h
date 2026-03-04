#include "../core/clipper_test.h"
#include "minimol_io_gemmi.h"
#include <gemmi/symmetry.hpp>

namespace clipper {
  class Test_minimol_gemmi : public Test_base {
    public:
      bool run(const String &filename);
    protected:
      bool test(const String &id, const Spacegroup &sg1, const gemmi::SpaceGroup &sg2);
      bool test(const String &id, const MiniMol &mmol, const GemmiModel &gmol);
  };
}