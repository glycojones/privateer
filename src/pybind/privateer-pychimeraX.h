#include <privateer-lib.h>
#include <privateer-json.h>

struct TorsionEntry {
    std::string sugar_1;
    std::string sugar_2;
    std::string atom_number_1;
    std::string atom_number_2;
    float phi;
    float psi;
};

struct TableEntry
{
    std::string svg;
    std::string wurcs;
    std::string chain;
    std::string glyconnect_id = "NotFound";
    std::string glytoucan_id = "NotFound";
    std::string id;
    int torsion_err = 0;
    int conformation_err = 0;
    int anomer_err = 0;
    int puckering_err = 0;
    int chirality_err = 0;

    std::vector<TorsionEntry> torsions;
};

std::vector<TableEntry> validate(const std::string &file, const std::string &name);