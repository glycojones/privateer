/*! \file clipper-glyco_data.cpp
Implementation file for the sugar database */

// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York


#include "clipper-glyco_data.h"

namespace clipper
{

namespace data
{

const clipper::String conformational_landscape[] =
{
    "Pln", "4c1", "1c4",
    "3Ob", "b25", "14b", "b3O", "25b", "b14",
    "Oev", "ev5", "4ev", "ev3", "2ev", "ev1", "3ev", "ev2", "1ev", "evO", "5ev", "ev4",
    "Oh5", "4h5", "4h3", "2h3", "2h1", "Oh1", "3h2", "1h2", "1hO", "5hO", "5h4", "3h4",
    "Os2", "1s5", "1s3", "2sO", "5s1", "3s1",
    "3t2", "3ev", "3t4", "ev4", "Ot4", "Oev", "Ot1", "ev1", "2t1", "2ev", "2t3", "ev3", "4t3", "4ev", "4tO", "evO", "1tO", "1ev", "1t2", "ev2"
};

const clipper::String iupac_conformational_landscape[] =
{
    "Planar ring", "<sup>4</sup>C<sub>1</sub>", "<sup>1</sup>C<sub>4",
    "<sup>3,O</sup>B", "B<sub>2,5</sub>", "<sup>1,4</sup>B", "B<sub>3,O</sub>", "<sup>2,5</sup>B", "B<sub>1,4</sub>",
    "<sup>O</sup>E", "E<sub>5</sub>", "<sup>4</sup>E", "E<sub>3</sub>", "<sup>2</sup>E", "E<sub>1</sub>", "<sup>3</sup>E",
    "E<sub>2</sub>", "<sup>1</sup>E", "E<sub>O</sub>", "<sup>5</sup>E", "E<sub>4</sub>",
    "<sup>O</sup>H<sub>5</sub>", "<sup>4</sup>H<sub>5</sub>", "<sup>4</sup>H<sub>3</sub>", "<sup>2</sup>H<sub>3</sub>",
    "<sup>2</sup>H<sub>1</sub>",
    "<sup>O</sup>H<sub>1</sub>", "<sup>3</sup>H<sub>2</sub>", "<sup>1</sup>H<sub>2</sub>", "<sup>1</sup>H<sub>O</sub>",
    "<sup>5</sup>H<sub>O</sub>",
    "<sup>5</sup>H<sub>4</sub>", "<sup>3</sup>H<sub>4</sub>",
    "<sup>O</sup>S<sub>2</sub>", "<sup>1</sup>S<sub>5</sub>", "<sup>1</sup>S<sub>3</sub>", "<sup>2</sup>S<sub>O</sub>",
    "<sup>5</sup>S<sub>1</sub>", "<sup>3</sup>S<sub>1</sub>",
    "<sup>3</sup>T<sub>2</sub>", "<sup>3</sup>E","<sup>3</sup>T<sub>4</sub>", "E<sub>4</sub>","<sup>O</sup>T<sub>4</sub>",
    "<sup>O</sup>E","<sup>O</sup>T<sub>1</sub>", "E<sub>1</sub>","<sup>2</sup>T<sub>1</sub>", "<sup>2</sup>E",
    "<sup>2</sup>T<sub>3</sub>", "E<sub>3</sub>","<sup>4</sup>T<sub>3</sub>", "<sup>4</sup>E","<sup>4</sup>T<sub>O</sub>",
    "E<sub>O</sub>","<sup>1</sup>T<sub>O</sub>", "<sup>1</sup>E","<sup>1</sup>T<sub>2</sub>", "E<sub>2</sub>"
};

const disaccharide disaccharide_database[] =
{
    { "SUC" ,  "Sucrose",
        { "FRU" ,	 "B", 	 "D", 	 "FRUCTOSE" ,                                     "O2' C2' C3' C4' C5'", 0.377, "ev4", 0.016, 5.815 },
        { "GLC" ,	 "A", 	 "D", 	 "ALPHA-D-GLUCOSE" ,                                "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.358 }
    },

    { "LAT" ,  "Beta-lactose",
        { "GAL" ,        "B",    "D",   "Beta-D-galactopyranose" ,                    "O5' C1' C2' C3' C4' C5'", 0.621, "4c1", 0.003, 2.036 },
        { "BGC" ,        "B",    "D",   "Beta-D-glucopyranose" ,                            "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.036 }
    },

    { "LBT" ,  "Alpha-lactose",
        { "GAL" ,        "B",    "D",   "Beta-D-galactopyranose" ,                          "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.036 },
        { "GLC" ,        "A",    "D",   "Alpha-D-glucopyranose" ,                     "O5' C1' C2' C3' C4' C5'", 0.621, "4c1", 0.003, 2.036 }
    },

    { "CBI" ,  "Cellobiose",
        { "BGC" ,        "B",    "D",   "Beta-D-glucopyranose" ,                     "O5' C1' C2' C3' C4' C5'", 0.621, "4c1", 0.003, 2.036 },
        { "BGC" ,        "B",    "D",   "Beta-D-glucopyranose" ,                           "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.036 }
    },

    { "MAL" ,  "Maltose",
        { "GLC" ,        "A",    "D",   "Alpha-D-glucopyranose" ,                     "O5' C1' C2' C3' C4' C5'", 0.621, "4c1", 0.003, 2.036 },
        { "GLC" ,        "A",    "D",   "Alpha-D-glucopyranose" ,                           "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.036 }
    },

    { "TRE" ,  "Trehalose",
        { "GLC" ,        "A",    "D",   "Alpha-D-glucopyranose" ,                     "O5P C1P C2P C3P C4P C5P", 0.621, "4c1", 0.003, 2.036 },
        { "GLC" ,        "A",    "D",   "Alpha-D-glucopyranose" ,                           "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.036 }
    },

    { "BXP" ,  "Xylobiose",
        { "XYP" ,        "B",    "D",   "Beta-D-xylopyranose" ,                       "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.036 },
        { "XYP" ,        "B",    "D",   "Beta-D-xylopyranose" ,                       "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.036 }
    }
};

const int disaccharide_database_size = sizeof( disaccharide_database ) / sizeof( disaccharide_database[0] );

const sugar_database_entry sugar_database[] =
{
    // TO DO: ADD BM7, currently does not exist in CCP4 monomer library
    { "13A" ,    "B",    "L",    "7-(3,4-DIHYDROXY-5R-HYDROXYMETHYLTETRAHYDROFU" ,  "O1 C2 C3 C4 C5",    0.380, "4ev", 0.017, 5.882 },
    { "145" ,	 "B", 	 "D", 	 "1-O-[O-NITROPHENYL]-BETA-D-GALACTOPYRANOSE" ,     "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.051 },
    { "147" ,	 "B", 	 "D", 	 "1-O-[P-NITROPHENYL]-BETA-D-GALACTOSE" ,           "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.051 },
    { "149" ,	 "A", 	 "D", 	 "D-GALACTONOLACTONE" ,				    "O5 C1 C2 C3 C4 C5", 0.553, "4h3", 0.038, 7.168 },
    { "15L" ,	 "A", 	 "D", 	 "GALACTARO-1,5-LACTONE" , 			    "O5 C1 C2 C3 C4 C5", 0.552, "4c1", 0.018, 4.147 },
    { "16G" ,	 "A", 	 "D", 	 "N-ACETYL-D-GLUCOSAMINE-6-PHOSPHATE" , 	    "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.405 },
    { "18D" ,	 "A", 	 "L", 	 "3,5-DIDEOXY-5-(PROPANOYLAMINO)-D-GLYCERO-ALPH" ,  "O6 C2 C3 C4 C5 C6", 0.592, "1c4", 0.001, 1.404 },
    { "1GL" ,	 "A", 	 "D", 	 "4-O-METHYL-2,6-DIDEOXY-ALPHA-D-GALACTO-HEXOPY" ,  "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.037 },
    { "1GN" ,	 "B", 	 "D", 	 "2-DEOXY-2-AMINOGALACTOSE" , 			    "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.004, 2.022 },
    { "289" ,	 "B", 	 "D", 	 "D-GLYCERO-ALPHA-D-MANNO-HEPTOPYRANOSE" , 	    "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.387 },
    { "291" ,	 "A", 	 "D", 	 "PROP-2-EN-1-YL7-O-CARBAMOYL-L-GLYCERO-ALPHA-" ,   "O5 C1 C2 C3 C4 C5", 0.592, "4c1", 0.001, 1.363 },
    { "293" ,	 "A", 	 "D", 	 "(2S,4R,5S,6R)-6-((S)-1,2-DIHYDROXYETHYL)TETRAH" , "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.386 },
    { "2DG" ,	 "A", 	 "D", 	 "2-DEOXY-BETA-D-GALACTOSE" , 	 		    "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.036 },
    { "2DR" ,	 "B", 	 "D", 	 "2-DEOXY-BETA-D-ERYTHRO-PENTOFURANOSE" , 	    "O4 C1 C2 C3 C4",    0.379, "ev1", 0.014, 5.461 },
    { "2FG" ,	 "B", 	 "D", 	 "2-FLUORO-2-DEOXY-BETA-D-GALACTOPYRANOSE" , 	    "O5 C1 C2 C3 C4 C5", 0.594, "4c1", 0.002, 1.388 },
    { "2GL" ,	 "A", 	 "D", 	 "4-O-ACETYL-2,6-DIDEOXY-ALPHA-D-GALACTO-HEXOPY" ,  "O5 C1 C2 C3 C4 C5", 0.533, "1c4", 0.004, 2.549 },
    { "2GS" ,	 "A", 	 "D", 	 "2-O-METHYL-ALPHA-D-GALACTOPYRANOSE" , 	    "O5 C1 C2 C3 C4 C5", 0.532, "1c4", 0.004, 2.485 },
    { "2M5" ,	 "A", 	 "D", 	 "METHYL7-DEOXY-L-GLYCERO-ALPHA-D-MANNO-HEPTOP" ,   "O5 C1 C2 C3 C4 C5", 0.592, "4c1", 0.001, 1.373 },
    { "3FM" ,	 "A", 	 "D", 	 "3-O-FORMAMIDO-ALPHA-D-MANNOPYRANOSIDE" ,          "O5 C1 C2 C3 C4 C5", 0.619, "4c1", 0.003, 1.995 },
    { "3HD" ,	 "N", 	 "D", 	 "3-O-METHYL-O-ALPHA-D-MANNOPYRANOSYL" ,            "O5 C1 C2 C3 C4 C5", 0.619, "4c1", 0.003, 1.982 },
    { "3MG" ,	 "B", 	 "D", 	 "3-O-METHYL-BETA-D-GLUCOPYRANOSE" ,                "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.374 },
    { "42D" ,	 "A", 	 "L", 	 "3,5-DIDEOXY-5-[(METHOXYCARBONYL)AMINO]-D-GLYC" ,  "O6 C2 C3 C4 C5 C6", 0.593, "1c4", 0.002, 1.374 },
    { "445" ,	 "B", 	 "D", 	 "N-[OXO(PHENYLAMINO)ACETYL]-BETA-D-GLUCOPYRANO" ,  "O5 C1 C2 C3 C4 C5", 0.594, "4c1", 0.001, 1.352 },
    { "46M" ,	 "N", 	 "D", 	 "(4AR,6R,7S,8R,8AS)-HEXAHYDRO-6,7,8-TRIHYDROXY" ,  "O5 C1 C2 C3 C4 C5", 0.626, "4c1", 0.004, 2.322 },
    { "475" ,	 "B", 	 "D", 	 "N-[OXO(PYRIDIN-2-YLAMINO)ACETYL]-BETA-D-GLUCO" ,  "O5 C1 C2 C3 C4 C5", 0.594, "4c1", 0.001, 1.380 },
    { "49A" ,	 "N", 	 "L", 	 "4,9-AMINO-2,4-DEOXY-2,3-DEHYDRO-N-ACETYL-NEUR" ,  "O6 C2 C3 C4 C5 C6", 0.340, "5h4", 0.087, 9.082 },
    { "4AM" ,	 "N", 	 "L", 	 "N4-AMINO-2-DEOXY-2,3-DEHYDRO-N-NEURAMINICACID" ,  "O6 C2 C3 C4 C5 C6", 0.340, "5h4", 0.087, 9.071 },
    { "4GP" ,	 "B", 	 "D", 	 "N-(BETA-D-GLUCOPYRANOSYL)OXAMICACID" ,            "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.075 },
    { "4GL" ,	 "A", 	 "D", 	 "alpha-D-gulopyranose" ,                           "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.382 },
    { "4NN" ,	 "N", 	 "D", 	 "N-[(5S,6R)-5-hydroxy-6-(hydroxymethyl)-2-oxo-5,6-dihydro-2H-pyran-3-yl]acetamide" ,    "O5 C1 C2 C3 C4 C5", 0.621, "5ev", 0.003, 2.016 },
    { "5N6" ,	 "A", 	 "D", 	 "9-O-acetyl-5-acetamido-3,5-dideoxy-D-glycero-alpha-D-galacto-non-2-ulopyranosonic acid" , "O6 C2 C3 C4 C5 C6", 0.593, "1c4", 0.001, 1.379 }, // new entry 01/12/2021
    { "6GP" ,	 "B", 	 "D", 	 "METHYL-N-(BETA-D-GLUCOPYRANOSYL)OXAMATE" ,        "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.001 },
    { "6LW" ,	 "N", 	 "D", 	 "(Z)-L-Arabinonhydroximo-1,4-lactone" , 			"O4 C1 C2 C3 C4", 0.552, "ev3", 0.018, 4.147 },
    { "6MN" ,	 "A", 	 "D", 	 "2-AMINO-2-DEOXY-6-O-PHOSPHONO-ALPHA-D-MANNOPY" ,  "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.369 },
    { "6PZ" ,	 "B", 	 "D", 	 "5,7-bisacetamido-3,5,7,9-tetradeoxy-L-glycero-alpha-L-manno-non-2-ulopyranosonic acid" , "O6 C2 C3 C4 C5 C6", 0.593, "1c4", 0.001, 1.379 }, // new entry 01/12/2021
    { "7CV" ,    "A",    "L",    "6-deoxy-2,3-di-O-methyl-alpha-L-mannopyranose",   "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.369 },
    { "7D1" ,    "N",    "D",    "1,2-Dideoxy-D-mannose",                           "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.369 },
    { "7LQ" ,    "N",    "D",    "GLC-CHEX",                                        "C7 C1 C2 C3 C4 C5", 0.593, "4h5", 0.003, 1.304 },
    { "7GP" ,	 "B", 	 "D", 	 "ETHYL-N-(BETA-D-GLUCOPYRANOSYL)OXAMATE" ,         "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.011 },
    { "7JZ" ,	 "B", 	 "D", 	 "2-DEOXY-2,2-DIFLUORO-BETA-D-LYXO-HEXOPYRANOSE" ,  "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.342 },
    { "83Y" ,	 "A", 	 "L", 	 "3-O-sulfo-alpha-L-rhamnopyranose" ,               "O5 C1 C2 C3 C4 C5", 0.620, "1c4", 0.003, 2.041 },
    { "8EX" ,	 "B", 	 "D", 	 "2-acetamido-2-deoxy-4,6-di-O-sulfo-beta-D-galactopyranose", "O5 C1 C2 C3 C4 C5", 0.563, "4c1", 0.002, 1.715 },
    { "8GP" ,	 "B", 	 "D", 	 "N-(BETA-D-GLUCOPYRANOSYL)-N'-CYCLOPROPYLOXAL" ,   "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.006 },
    { "9AM" ,	 "A", 	 "D", 	 "9-AMINO-2-DEOXY-2,3-DEHYDRO-N-ACETYL-NEURAMIN" ,  "O6 C2 C3 C4 C5 C6", 0.541, "5h4", 0.087, 9.084 },
    { "9RN" ,	 "A", 	 "D", 	 "3,6-anhydro-D-galactose" ,                        "O5 C1 C2 C3 C4 C5", 0.695, "4c1", 0.016, 4.879 },
    { "9WJ" ,	 "B", 	 "D", 	 "4-acetamido-4,6-dideoxy-beta-D-mannopyranose" ,   "O5 C1 C2 C3 C4 C5", 0.548, "4c1", 0.002, 1.803 },
    { "A2G" ,    "A",    "D",    "N-ACETYL-2-DEOXY-2-AMINO-GALACTOSE",              "O5 C1 C2 C3 C4 C5" , 0.593, "4c1", 0.002, 1.350 },
    { "A6P" ,	 "A", 	 "D", 	 "6-O-PHOSPHONO-ALPHA-D-ALLOPYRANOSE" ,             "O5 C1 C2 C3 C4 C5", 0.592, "4c1", 0.002, 1.380 },
    { "AAL" ,	 "A", 	 "L", 	 "3,6-ANHYDRO-L-GALACTOSE" ,                        "O5 C1 C2 C3 C4 C5", 0.695, "4c1", 0.016, 4.879 },
    { "ABE" ,	 "A", 	 "D", 	 "ABEQUOSE" ,                                       "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.023 },
    { "AC1" ,	 "A", 	 "D", 	 "4,6-dideoxy-4-{[(1S,4R,5S,6S)-4,5,6-trihydroxy-3-(hydroxymethyl)cyclohex-2-en-1-yl]amino}-alpha-D-glucopyranose", "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.023 },
    { "ADA" ,	 "A", 	 "D", 	 "ALPHA-D-GALACTOPYRANURONICACID" ,                 "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 1.994 },
    { "AFP" ,	 "A", 	 "D", 	 "ALPHAFRUCTOSE1,6-DIPHOSPHATE" ,                   "O5 C2 C3 C4 C5"   , 0.380, "Oev", 0.015, 5.686 },
    { "AGL" ,	 "A", 	 "D", 	 "4,6-DIDEOXY-4-AMINO-ALPHA-D-GLUCOSE" ,            "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 2.027 },
    { "AH8" ,	 "A", 	 "L", 	 "(2R,3R,4R,5S)-2-AZIDO-5-(HYDROXYMETHYL)TETRAH" ,  "O4 C1 C2 C3 C4"   , 0.382, "evO", 0.014, 5.711 },
    { "AHG" ,	 "B", 	 "D", 	 "2,5-ANHYDROGLUCITOL-1,6-BIPHOSPHATE" ,            "O5 C2 C3 C4 C5"   , 0.379, "ev1", 0.013, 5.447 },
    { "AHR" ,	 "A", 	 "L", 	 "ALPHA-L-ARABINOFURANOSE" ,                        "O4 C1 C2 C3 C4", 0.371, "ev3", 0.013, 5.447 },
    { "AIG" ,	 "B", 	 "D", 	 "hexyl 3-amino-3-deoxy-beta-D-galactopyranoside" , "O5 C1 C2 C3 C4 C5", 0.623, "4c1", 0.003, 2.051 },
    { "AOG" ,	 "B", 	 "D", 	 "octyl 3-amino-3-deoxy-beta-D-galactopyranoside" , "O5 C1 C2 C3 C4 C5", 0.623, "4c1", 0.003, 2.051 },
    { "AMG" ,	 "A", 	 "D", 	 "ALPHA-METHYL-D-GALACTOSIDE" ,                     "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 2.053 },
    { "AMN" ,	 "A", 	 "L", 	 "9-DEOXY-9-AMINO-2-O-METHYL-5-N-ACETYL-ALPHA-D" ,  "O6 C2 C3 C4 C5 C6", 0.620, "1c4", 0.003, 2.024 },
    { "AMP" ,    "B",    "D",    "ADENOSINE MONOPHOSPHATE",                       "O4' C1' C2' C3' C4'", 0.450, "3t2", 0.031, 4.000 },
    { "AMU" ,	 "B", 	 "N", 	 "BETA-N-ACETYLMURAMICACID" ,                       "O5 C1 C2 C3 C4 C5", 0.619, "4c1", 0.003, 2.023 },
    { "AMV" ,	 "B", 	 "D", 	 "METHYL2-(ACETYLAMINO)-3-O-[(1R)-1-CARBOXYETH" ,   "O5 C1 C2 C3 C4 C5", 0.558, "4c1", 0.002, 1.909 },
    { "ANA" ,	 "A", 	 "L", 	 "4-O-ACETYL-ALPHA-2-OMETHYL-5-N-ACETYL-ALPHA-D" ,  "O6 C2 C3 C4 C5 C6", 0.622, "1c4", 0.003, 2.051 },
    { "AQA" ,	 "B", 	 "L", 	 "4-deoxy-beta-L-threo-hex-4-enopyranuronic acid" , "O5 C1 C2 C3 C4 C5", 0.512, "2h1", 0.095, 7.411 },
    { "ARA" ,	 "A", 	 "L", 	 "ALPHA-L-ARABINOSE" ,                              "O5 C1 C2 C3 C4 C5", 0.623, "4c1", 0.003, 2.036 },
    { "ARB" ,	 "B", 	 "L", 	 "BETA-L-ARABINOSE" ,                               "O5 C1 C2 C3 C4 C5", 0.624, "4c1", 0.003, 2.029 },
    { "ARI" ,	 "B", 	 "D", 	 "(2R,3R,6R)-6-hydroxy-2-methyltetrahydro-2H-pyran-3-yl acetate",  "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.037 }, // new entry 01/12/2021
    { "ARW" ,	 "B", 	 "D", 	 "METHYLBETA-D-ARABINOPYRANOSIDE" ,                 "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.031 },
    { "ASG" ,	 "B", 	 "D", 	 "2-DEOXY-2-ACETAMIDO-BETA-D-GALACTOSE-4-SULFATE" , "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.031 },
    { "ASO" ,	 "N", 	 "D", 	 "1,5-ANHYDROSORBITOL" ,                            "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.032 },
    { "AXP" ,	 "A", 	 "L", 	 "4-ACETAMIDO-2,4-DIDEXOY-D-GLYCERO-BETA-D-GALA" ,  "O6 C2 C3 C4 C5 C6", 0.621, "4c1", 0.002, 2.023 },
    { "AXR" ,	 "A", 	 "D", 	 "METHYLALPHA-D-ARABINOFURANOSIDE" ,                "O4 C1 C2 C3 C4"   , 0.379, "1ev", 0.014, 5.457 },
    { "B16" ,	 "B", 	 "D", 	 "1,6-DI-O-PHOSPHONO-BETA-D-GLUCOPYRANOSE" ,        "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.377 },
    { "B6D" ,	 "B", 	 "D", 	 "2,4-bisacetamido-2,4,6-trideoxy-beta-D-glucopyranose", "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.022 },
    { "B7G" ,    "B",    "D",    "HEPTYL-BETA-D-GLUCOPYRANOSIDE",                   "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.033 },
    { "BBV" ,	 "A", 	 "D", 	 "benzyl 2-acetamido-2-deoxy-alpha-D-glucopyranoside","O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.377 },
    { "BDF" ,	 "B", 	 "D", 	 "BETA-D-FRUCTOPYRANOSE" ,                          "O6 C2 C3 C4 C5 C6", 0.621, "4c1", 0.003, 2.036 },
    { "BDG" ,	 "A", 	 "D", 	 "O-2,6-DIAMINO-2,6-DIDEOXY-ALPHA-D-GLUCOPYRANOSE" ,"O5 C1 C2 C3 C4 C5", 0.625, "4c1", 0.003, 2.036 },
    { "BDP" ,	 "B", 	 "D", 	 "BETA-D-GLUCOPYRANURONICACID" ,                    "O5 C1 C2 C3 C4 C5", 0.592, "4c1", 0.001, 1.376 },
    { "BDR" ,	 "B", 	 "D", 	 "beta-D-ribofuranose" ,                            "O4 C1 C2 C3 C4"   , 0.372, "3ev", 0.015, 5.717 },
    { "BEM" ,	 "B", 	 "D", 	 "BETA-D-MANNURONICACID" ,                          "O5 C1 C2 C3 C4 C5", 0.594, "4c1", 0.002, 1.366 },
    { "BG6" ,	 "B", 	 "D", 	 "BETA-D-GLUCOSE-6-PHOSPHATE" ,                     "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.019 },
    { "BGC" ,	 "B", 	 "D", 	 "BETA-D-GLUCOSE" ,                                 "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.036 },
    { "BGN" ,	 "B", 	 "D", 	 "N-BUTANOYL-2-AMINO-2-DEOXY-GLUCOPYRANOSIDE" ,     "O5 C1 C2 C3 C4 C5", 0.619, "4c1", 0.003, 2.026 },
    { "BGP" ,	 "B", 	 "D", 	 "BETA-GALACTOSE-6-PHOSPHATE" ,                     "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.003 },
    { "BGL" ,    "B",    "D",    "B-2-OCTYLGLUCOSIDE",                              "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.000 },
    { "BGS" ,	 "B", 	 "D", 	 "BETA-D-GLUCOPYRANOSYLSULFONYLETHANE" ,            "O5 C1 C2 C3 C4 C5", 0.623, "4c1", 0.003, 2.051 },
    { "BHG" ,	 "B", 	 "D", 	 "2-HEXYLOXY-6-HYDROXYMETHYL-TETRAHYDRO-PYRA" ,     "O5 C1 C2 C3 C4 C5", 0.623, "4c1", 0.003, 2.051 },
    { "BM3" ,	 "A", 	 "D", 	 "2-(ACETYLAMINO)-2-DEOXY-ALPHA-D-MANNOPYRANOSE" ,  "O5 C1 C2 C3 C4 C5", 0.548, "4c1", 0.002, 1.803 },
    { "BMA" ,	 "B", 	 "D", 	 "BETA-D-MANNOSE" ,                                 "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.029 },
    { "BMX" ,	 "A", 	 "D", 	 "2-(ACETYLAMINO)-2-DEOXY-6-O-PHOSPHONO-ALPHA-D" ,  "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.382 },
    { "BNG" ,    "B",    "D",    "B-NONYLGLUCOSIDE",                                "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 2.003 },
    { "BOG" ,    "B",    "D",    "B-OCTYLGLUCOSIDE",                                "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 1.803 },
    { "BXX" ,	 "B", 	 "D", 	 "BETA-D-ARABINOFURANOSE" ,                         "O4 C1 C2 C3 C4"   , 0.382, "Oev", 0.014, 5.717 },
    { "BXY" ,	 "A", 	 "D", 	 "ALPHA-D-ARABINOFURANOSE" ,                        "O4 C1 C2 C3 C4"   , 0.380, "1ev", 0.014, 5.461 },
    { "C3X" ,	 "B", 	 "D", 	 "2,3-EPOXYPROPYL-BETA-D-XYLOSIDE" ,                "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.028 },
    { "C4X" ,	 "B", 	 "D", 	 "3,4-EPOXYBUTYL-BETA-D-XYLOSIDE" ,                 "O5 C1 C2 C3 C4 C5", 0.619, "4c1", 0.003, 2.031 },
    { "C5X" ,	 "B", 	 "D", 	 "4,5-EPOXYPENTYL-BETA-D-XYLOSIDE" ,                "O5 C1 C2 C3 C4 C5", 0.619, "4c1", 0.003, 2.031 },
    { "CBF" ,	 "A", 	 "D", 	 "C-(1-HYDROGYL-BETA-D-GLUCOPYRANOSYL)FORMAMID" ,   "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.400 },
    { "CDG" ,	 "B", 	 "D", 	 "METHYL4,6-O-[(1R)-1-CARBOXYETHYLIDENE]-BETA-" ,   "O5 C1 C2 C3 C4 C5", 0.612, "4c1", 0.006, 1.980 },
    { "CDR" ,	 "B", 	 "D", 	 "(2R,5S,6R)-6-methyltetrahydro-2H-pyran-2,5-diol",  "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.037 }, // new entry 01/12/2021
    { "CEG" ,	 "B", 	 "D", 	 "4,6-O-(1-CARBOXYETHYLIDENE)-BETA-D-GLUCOSE" ,     "O5 C1 C2 C3 C4 C5", 0.626, "4c1", 0.004, 2.299 },
    { "CGF" ,	 "B", 	 "D", 	 "C-(1-AZIDO-ALPHA-D-GLUCOPYRANOSYL)FORMAMIDE" ,    "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.016 },
    { "CNP" ,	 "N", 	 "N", 	 "2-PROPENYL-N-ACETYL-NEURAMICACID" ,               "O6 C2 C3 C4 C5 C6", 0.620, "1c4", 0.003, 2.003 },
    { "CR1" ,	 "B", 	 "D", 	 "1-DEOXY-1-METHOXYCARBAMIDO-BETA-D-GLUCOPYRANO" ,  "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.010 },
    { "CR6" ,	 "B", 	 "D", 	 "1-DEOXY-1-ACETYLAMINO-BETA-D-GLUCO-2-HEPTULOP" ,  "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.009 },
    { "CRA" ,	 "B", 	 "D", 	 "1-DEOXY-1-METHOXYCARBAMIDO-BETA-D-GLUCO-2-HEP" ,  "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.035 },
    { "CYH" ,	 "N", 	 "N", 	 "CYCLOHEXANONE"                                 ,  "C6 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.035 },
    { "D6G" ,	 "A", 	 "D", 	 "2-DEOXY-6-O-PHOSPHONO-ALPHA-D-ARABINO-HEXOPYR" ,  "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.061 },
    { "DAG" ,	 "B", 	 "D", 	 "4,6-DIDEOXY-4-AMINO-BETA-D-GLUCOPYRANOSIDE" ,     "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.043 },
    { "DAN" ,	 "N", 	 "L", 	 "2-DEOXY-2,3-DEHYDRO-N-ACETYL-NEURAMINICACID" ,    "O6 C2 C3 C4 C5 C6", 0.540, "5h4", 0.087, 9.083 },
    { "DDA" ,	 "B", 	 "D", 	 "2,6-DIDEOXY-BETA-D-GLUCOSE" ,                     "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.052 },
    { "DDL" ,	 "B", 	 "D", 	 "2,6-DIDEOXY-BETA-D-GALACTOSE" ,                   "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 2.036 },
    { "DEG" ,	 "A", 	 "D", 	 "BUTYLALPHA-D-MANNOPYRANOSIDE" ,                   "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.001, 2.060 },
    { "DFX" ,	 "N", 	 "D", 	 "1,5-anhydro-2-deoxy-2-fluoro-D-xylitol" ,         "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.032 },
    { "DGO" ,	 "N", 	 "D", 	 "D-glucal" ,                                       "O5 C1 C2 C3 C4 C5", 0.620, "4ev", 0.003, 2.032 },
    { "DGS" ,	 "A", 	 "D", 	 "3,6-ANHYDRO-D-GALACTOSE-2-SULFATE" ,              "O5 C1 C2 C3 C4 C5", 0.696, "1c4", 0.017, 4.890 },
    { "DK4" ,	 "B", 	 "D", 	 "1-(3-DEOXY-3-FLUORO-BETA-D-GLUCOPYRANOSYL)-5-" ,  "O5 C1 C2 C3 C4 C5", 0.592, "4c1", 0.001, 1.372 },
    { "DK5" ,	 "B", 	 "D", 	 "1-(2,3-DIDEOXY-3-FLUORO-BETA-D-ARABINO-HEXOPY" ,  "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.394 },
    { "DKX" ,	 "B", 	 "D", 	 "1-(3-DEOXY-3-FLUORO-BETA-D-GLUCOPYRANOSYL)PYR" ,  "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.353 },
    { "DKY" ,	 "B", 	 "D", 	 "1-(3-DEOXY-3-FLUORO-BETA-D-GLUCOPYRANOSYL)-4-" ,  "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.385 },
    { "DKZ" ,	 "B", 	 "D", 	 "4-AMINO-1-(3-DEOXY-3-FLUORO-BETA-D-GLUCOPYRAN" ,  "O5 C1 C2 C3 C4 C5", 0.594, "4c1", 0.002, 1.390 },
    { "DL6" ,	 "B", 	 "D", 	 "2-AZIDO-N-((2R,3R,4S,5S,6R)-3,4,5-TRIHYDROXY-" ,  "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.002, 2.072 },
    { "DLF" ,	 "A", 	 "L", 	 "2-DEOXY-ALPHA-L-FUCOPYRANOSIDE" ,                 "O5 C1 C2 C3 C4 C5", 0.621, "1c4", 0.004, 2.048 },
    { "DLG" ,	 "B", 	 "D", 	 "HEXYL3-DEOXY-BETA-D-GALACTOPYRANOSE" ,            "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.012 },
    { "DO8" ,	 "A", 	 "D", 	 "3-DEOXY-D-MANNO-2-OCTULOSONATE-8-PHOSPHATE" ,     "O6 C2 C3 C4 C5 C6", 0.622, "4c1", 0.003, 2.046 },
    { "DRI" ,	 "B", 	 "D", 	 "4-O-METHYL-2,6-DIDEOXY-BETA-D-GLUCOSE" ,          "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.002, 2.045 },
    { "DSR" ,	 "B", 	 "D", 	 "2,6-DIDEOXY-4-THIO-BETA-D-ALLOSEPYRANOSIDE" ,     "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.015 },
    { "DT6" ,	 "B", 	 "D", 	 "2,4-bisacetamido-2,4-dideoxy-beta-D-glucopyranose","O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.022 },
    { "DVC" ,	 "N", 	 "N", 	 "(2R,4S,6S)-4-AZANYL-4,6-DIMETHYL-OXANE-2,5,5-" ,  "O5 C1 C2 C3 C4 C5", 0.592, "1c4", 0.001, 1.368 },
    { "E5G" ,	 "A", 	 "D", 	 "5-HYDROXYPENTYLALPHA-D-GLUCOPYRANOSIDE" ,         "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.377 },
    { "EAG" ,	 "B", 	 "D", 	 "2-AMINOETHYL2-(ACETYLAMINO)-2-DEOXY-BETA-D-G" ,   "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.355 },
    { "EAT" ,    "A",    "L",    "EAT",                                           "NAH CAE CAG CAA CAF", 0.300, "ev3", 0.030, 7.980 },
    { "EBG" ,	 "A", 	 "D", 	 "2-HYDROXYMETHYL-6-(2-OXIRANYL-ETHOXY)-TETRAHY" ,  "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.054 },
    { "EBQ" ,	 "A", 	 "D", 	 "2-[(2R)-OXIRAN-2-YL]ETHYLALPHA-D-GLUCOPYRANO" ,   "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.377 },
    { "EGA" ,	 "B", 	 "D", 	 "ethyl beta-D-galactopyranoside" ,                 "O5 C1 C2 C3 C4 C5", 0.623, "4c1", 0.003, 2.051 },
    { "EMP" ,	 "A", 	 "L", 	 "2,4-DIDEOXY-4-(ETHYLAMINO)-3-O-METHYLALPHA-L" ,   "O5 C1 C2 C3 C4 C5", 0.549, "4c1", 0.003, 1.746 },
    { "EPG" ,	 "A", 	 "D", 	 "2-HYDROXYMETHYL-6-OXIRANYLMETHOXY-TETRAHYDRO-" ,  "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.068 },
    { "EQP" ,	 "A", 	 "L", 	 "(4-ACETAMIDO-2,4-DIDEOXY-D-GLYCERO-ALPHA-D-GA" ,  "O6 C2 C3 C4 C5 C6", 0.620, "1c4", 0.003, 1.970 },
    { "ERE" ,	 "A", 	 "L", 	 "(1R,3S,4R,5S)-3-AMINO-2,3,6-TRIDEOXY-3-METHYL" ,  "O5 C1 C2 C3 C4 C5", 0.593, "1c4", 0.001, 1.397 },
    { "ERI" ,	 "A", 	 "L", 	 "3-C-methyl-4-O-acetyl-alpha-L-Olivopyranose",     "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.037 }, // new entry 01/12/2021
    { "F1P" ,	 "B", 	 "D", 	 "1-O-PHOSPHONO-BETA-D-FRUCTOPYRANOSE" ,            "O6 C2 C3 C4 C5 C6", 0.621, "1c4", 0.003, 2.045 },
    { "F1X" ,	 "B", 	 "D", 	 "1-O-PHOSPHONO-BETA-D-FRUCTOFURANOSE" ,            "O5 C2 C3 C4 C5"   , 0.372, "3ev", 0.015, 5.727 },
    { "F6P" ,	 "B", 	 "D", 	 "FRUCTOSE-6-PHOSPHATE" , 	 		                "O5 C2 C3 C4 C5"   , 0.372, "3ev", 0.015, 5.714 },
    { "FBP" ,	 "B", 	 "D", 	 "BETA-FRUCTOSE-1,6-DIPHOSPHATE" ,                  "O5 C2 C3 C4 C5"   , 0.372, "3ev", 0.015, 5.718 },
    { "FCA" ,	 "A", 	 "D", 	 "ALPHA-D-FUCOSE" ,                                 "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.034 },
    { "FCB" ,	 "B", 	 "D", 	 "BETA-D-FUCOSE" ,                                  "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.021 },
    { "FDP" ,	 "B", 	 "D", 	 "FRUCTOSE-2,6-DIPHOSPHATE" ,                       "O5 C2 C3 C4 C5"   , 0.372, "3ev", 0.015, 5.714 },
    { "FRU" ,	 "B", 	 "D", 	 "FRUCTOSE" ,                                       "O5 C2 C3 C4 C5"   , 0.377, "ev4", 0.016, 5.815 },
    { "FUB" ,	 "B", 	 "L", 	 "beta-L-arabinofuranose" ,                         "O4 C1 C2 C3 C4", 0.371, "ev3", 0.013, 5.447 },
    { "FUC" ,	 "A", 	 "L", 	 "ALPHA-L-FUCOSE" ,                                 "O5 C1 C2 C3 C4 C5", 0.621, "1c4", 0.003, 2.034 },
    { "FUL" ,	 "B", 	 "L", 	 "BETA-L-FUCOSE" ,                                  "O5 C1 C2 C3 C4 C5", 0.621, "1c4", 0.003, 2.044 },
    { "FUY" ,	 "B", 	 "L", 	 "BETA-L-FUCOSYL-AZIDE" ,                           "O5 C1 C2 C3 C4 C5", 0.592, "1c4", 0.001, 1.376 },
    { "G0S" ,	 "B", 	 "D", 	 "3-(BETA-D-GALACTOPYRANOSYLTHIO)PROPANOICACID" ,   "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.374 },
    { "G16" ,	 "A", 	 "D", 	 "ALPHA-D-GLUCOSE1,6-BISPHOSPHATE" ,                "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.051 },
    { "G1P" ,	 "A", 	 "D", 	 "ALPHA-D-GLUCOSE-1-PHOSPHATE" ,                    "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.045 },
    { "G2F" ,	 "A", 	 "D", 	 "2-DEOXY-2-FLUORO-ALPHA-D-GLUCOPYRANOSE" ,         "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.045 },
    { "G3F" ,	 "B", 	 "D", 	 "3-DEOXY-3-FLUORO-BETA-D-GLUCOPYRANOSE" ,          "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.036 },
    { "G4D" ,	 "A", 	 "D", 	 "4-DEOXY-ALPHA-D-GLUCOSE" ,                        "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.047 },
    { "G4S" ,	 "B", 	 "D", 	 "4-O-SULFO-BETA-D-GALACTOPYRANOSE" ,               "O5 C1 C2 C3 C4 C5", 0.563, "4c1", 0.002, 1.715 },
    { "G6D" ,	 "A", 	 "D", 	 "alpha-D-quinovopyranose" ,                        "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.022 },
    { "G6P" ,	 "A", 	 "D", 	 "ALPHA-D-GLUCOSE-6-PHOSPHATE" ,                    "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.022 },
    { "G6S" ,	 "B", 	 "D", 	 "D-GALACTOSE-6-SULFATE" ,                          "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.363 },
    { "G7P" ,	 "B", 	 "D", 	 "6,7-DIDEOXY-7-PHOSPHONO-BETA-D-GLUCO-HEPTOPYR" ,  "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.377 },
    { "GAA" ,	 "A", 	 "D", 	 "METANITROPHENYL-ALPHA-D-GALACTOSIDE" ,            "O5 C1 C2 C3 C4 C5", 0.623, "4c1", 0.003, 2.051 },
    { "GAD" ,	 "N", 	 "L", 	 "2,6-anhydro-3-deoxy-L-threo-hex-2-enonic acid" ,  "O5 C1 C2 C3 C4 C5", 0.511, "2h1", 0.095, 7.421 }, // new entry 01/12/2021
    { "GAF" ,	 "A", 	 "D", 	 "2-DEOXY-2-FLUORO-ALPHA-D-GALACTOPYRANOSE" ,       "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.038 },
    { "GAL" ,	 "B", 	 "D", 	 "BETA-D-GALACTOSE" ,                               "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.015 },
    { "GAT" ,	 "A", 	 "D", 	 "4'-AMINOPHENYL-ALPHA-D-GALACTOPYRANOSIDE" ,       "O5 C1 C2 C3 C4 C5", 0.623, "4c1", 0.003, 2.051 },
    { "GLA" ,	 "A", 	 "D", 	 "ALPHAD-GALACTOSE" ,                               "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.023 },
    { "GC1" ,	 "A", 	 "L", 	 "2,6-anhydro-L-gulonic acid" ,                     "O5 C1 C2 C3 C4 C5", 0.340, "4c1", 0.087, 9.071 },
    { "GC4" ,	 "B", 	 "D", 	 "4-DEOXY-D-GLUCURONICACID" ,                       "O5 C1 C2 C3 C4 C5", 0.619, "4c1", 0.003, 2.028 },
    { "GCD" ,	 "A", 	 "L", 	 "4,5-DEHYDRO-D-GLUCURONICACID" ,                   "O5 C1 C2 C3 C4 C5", 0.511, "2h1", 0.095, 7.421 },
    { "GCN" ,	 "A", 	 "D", 	 "3-DEOXY-D-GLUCOSAMINE" ,                          "O5 C1 C2 C3 C4 C5", 0.619, "4c1", 0.003, 1.983 },
    { "GCS" ,	 "B", 	 "D", 	 "D-GLUCOSAMINE" ,                                  "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.034 },
    { "GCU" ,	 "A", 	 "D", 	 "D-GLUCURONICACID" ,                               "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.007 },
    { "GCV" ,	 "A", 	 "D", 	 "4-O-METHYL-ALPHA-D-GLUCURONICACID" , 	            "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 1.996 },
    { "GCW" ,	 "B", 	 "D", 	 "4-O-METHYL-BETA-D-GLUCURONICACID" ,               "O5 C1 C2 C3 C4 C5", 0.619, "4c1", 0.003, 2.014 },
    { "GDA" ,	 "B", 	 "D", 	 "4-DEOXY-4-AMINO-BETA-D-GLUCOSE" ,                 "O5 C1 C2 C3 C4 C5", 0.555, "4c1", 0.006, 1.810 },
    { "GDL" ,	 "A", 	 "D", 	 "2-ACETAMIDO-2-DEOXY-D-GLUCONO-1,5-LACTONE" ,      "O5 C1 C2 C3 C4 C5", 0.538, "4c1", 0.044, 5.629 },
    { "GFP" ,	 "A", 	 "D", 	 "2-DEOXY-2-FLUORO-ALPHA-D-GLUCOSE-1-PHOSPHATE" ,   "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.002, 2.044 },
    { "GL0" ,	 "B", 	 "D", 	 "BETA-D-GULOPYRANOSE" ,                            "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.382 },
    { "GL1" ,	 "A", 	 "D", 	 "1-O-PHOSPHONO-ALPHA-D-GALACTOPYRANOSE" ,          "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.031 },
    { "GL2" ,	 "N", 	 "D", 	 "3-AMINO-8,9,10-TRIHYDROXY-7-HYDROXYMETHYL-6-O" ,  "O5 C1 C2 C3 C4 C5", 0.612, "4c1", 0.003, 2.113 },
    { "GLA" ,	 "A", 	 "D", 	 "ALPHAD-GALACTOSE" ,                               "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.023 },
    { "GLC" ,	 "A", 	 "D", 	 "ALPHA-D-GLUCOSE" ,                                "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.358 },
    { "GLD" ,	 "A", 	 "D", 	 "4,6-DIDEOXYGLUCOSE" ,                             "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.044 },
    { "GLF" ,	 "A", 	 "D", 	 "1-FLUORO-ALPHA-1-DEOXY-D-GLUCOSE" ,               "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.044 },
    { "GLG" ,	 "A", 	 "D", 	 "ALPHA-D-GLUCOPYRANOSYL-2-CARBOXYLICACIDAMID" ,    "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.054 },
    { "GLP" ,	 "A", 	 "D", 	 "GLUCOSAMINE6-PHOSPHATE" ,                         "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.028 },
    { "GLS" ,	 "B", 	 "D", 	 "BETA-D-GLUCOPYRANOSESPIROHYDANTOIN" ,             "O5 C1 C2 C3 C4 C5", 0.613, "4c1", 0.003, 2.134 },
    { "GM0" ,	 "A", 	 "D", 	 "4-O-phosphono-L-glycero-alpha-D-manno-heptopyranose","O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.016 },
    { "GMB" ,	 "B", 	 "D", 	 "1,7-DI-O-PHOSPHONO-L-GLYCERO-BETA-D-MANNO-HEP" ,  "O5 C1 C2 C3 C4 C5", 0.594, "4c1", 0.002, 1.365 },
    { "GMH" ,	 "A", 	 "D", 	 "L-GLYCERO-D-MANNO-HEPTOPYRANOSE" ,                "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.016 },
    { "GNA" ,	 "B", 	 "L", 	 "2,4-DEOXY-4-GUANIDINO-5-N-ACETYL-NEURAMINICA" ,   "O6 C2 C3 C4 C5 C6", 0.620, "1c4", 0.003, 2.009 },
    { "GNS" ,	 "A", 	 "D", 	 "2-deoxy-2-(sulfoamino)-alpha-D-glucopyranose" ,   "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.045 },
    { "GNX" ,	 "A", 	 "D", 	 "2-deoxy-3-O-sulfo-2-(sulfoamino)-alpha-D-glucopyranose" , "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.045 },
    { "GP1" ,	 "A", 	 "D", 	 "2-amino-2-deoxy-1-O-phosphono-alpha-D-glucopyranose","O5 C1 C2 C3 C4 C5", 0.625, "4c1", 0.003, 2.036 },
    { "GP4" ,	 "A", 	 "D", 	 "GLUCOSAMINE4-PHOSPHATE" ,                   "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.028 },
    { "GPM" ,	 "N", 	 "D", 	 "(1S)-1,5-anhydro-1-(phosphonomethyl)-D-glucitol" ,"O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.377 },
    { "GS1" ,	 "B", 	 "D", 	 "1-THIO-BETA-D-GLUCOPYRANOSE" ,                    "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.029 },
    { "GTM" ,	 "B", 	 "D", 	 "O1-METHYL-4-DEOXY-4-THIO-BETA-D-GLUCOSE" ,        "O5 C1 C2 C3 C4 C5", 0.619, "4c1", 0.003, 1.992 },
    { "GTR" ,	 "B", 	 "D", 	 "GALACTURONIC ACID",                               "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.026 },
    { "GU0" ,	 "B", 	 "D", 	 "2,3,6-TRI-O-SULFONATO-ALPHA-L-GALACTOPYRANOSE" ,  "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.335 },
    { "GU1" ,	 "B", 	 "D", 	 "2,3-di-O-methyl-beta-D-glucopyranuronic acid" ,   "O5 C1 C2 C3 C4 C5", 0.557, "1c4", 0.002, 1.853 },
    { "GU3" ,	 "A", 	 "D", 	 "METHYL3-O-METHYL-2,6-DI-O-SULFO-ALPHA-D-GLUC" ,   "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.387 },
    { "GU4" ,	 "A", 	 "D", 	 "2,3,4,6-TETRA-O-SULFONATO-ALPHA-D-GLUCOPYRANO" ,  "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 2.039 },
    { "GU5" ,	 "A", 	 "D", 	 "2,3-DI-O-METHYL-6-O-SULFONATO-ALPHA-D-GLUCOPY" ,  "O5 C1 C2 C3 C4 C5", 0.518, "1c4", 0.004, 2.739 },
    { "GU6" ,	 "B", 	 "L", 	 "2,3,6-TRI-O-SULFONATO-ALPHA-D-GLUCOPYRANOSE" ,    "O5 C1 C2 C3 C4 C5", 0.557, "1c4", 0.002, 1.853 },
    { "GU8" ,	 "B", 	 "D", 	 "2,3,6-TRI-O-METHYL-BETA-D-GLUCOPYRANOSE" ,        "O5 C1 C2 C3 C4 C5", 0.557, "4c1", 0.001, 1.830 },
    { "GU9" ,	 "A", 	 "D", 	 "2,3,6-TRI-O-METHYL-ALPHA-D-GLUCOPYRANOSE" ,       "O5 C1 C2 C3 C4 C5", 0.518, "1c4", 0.004, 2.739 },
    { "GUP" ,	 "A", 	 "L", 	 "ALPHA-L-GULOPYRANOSIDE" ,                         "O5 C1 C2 C3 C4 C5", 0.621, "1c4", 0.003, 2.028 },
    { "GXL" ,	 "A", 	 "L", 	 "ALPHA-L-GALACTOPYRANOSE" ,                        "O5 C1 C2 C3 C4 C5", 0.621, "1c4", 0.003, 2.033 },
    { "GYP" ,	 "A", 	 "D", 	 "METHYL-ALPHA-D-GLUCOPYRANOSIDE" ,                 "O5 C1 C2 C3 C4 C5", 0.621, "1c4", 0.003, 2.031 },
    { "GZL" ,	 "B", 	 "D", 	 "beta-D-galactofuranose" ,                         "O4 C1 C2 C3 C4", 0.372, "ev3", 0.002, 2.051 }, // new entry 01/12/2021
    { "H1M" ,	 "A", 	 "D", 	 "METHYL2-DEOXY-2-(2-HYDROXYETHYL)-ALPHA-D-MAN" ,   "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.364 },
    { "H2P" ,	 "A", 	 "D", 	 "HEPTULOSE-2-PHOSPHATE" ,                          "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.072 },
    { "HSQ" ,	 "A", 	 "L", 	 "2-acetylamino-2-deoxy-alpha-L-idopyranose" ,      "O5 C1 C2 C3 C4 C5", 0.621, "1c4", 0.003, 2.037 },
    { "IDG" ,	 "B", 	 "L", 	 "O-2,6-DIAMINO-2,6-DIDEOXY-BETA-L-IDOPYRANOSE" ,   "O5 C1 C2 C3 C4 C5", 0.621, "1c4", 0.003, 2.037 },
    { "IDR" ,	 "A", 	 "L", 	 "L-IDURONICACID" ,                                 "O5 C1 C2 C3 C4 C5", 0.621, "1c4", 0.003, 2.021 },
    { "IDS" ,	 "A", 	 "L", 	 "2-O-SULFO-ALPHA-L-IDOPYRANURONICACID" ,           "O5 C1 C2 C3 C4 C5", 0.620, "1c4", 0.003, 2.016 },
    { "IDY" ,	 "A", 	 "L", 	 "1-O-methyl-2-O-sulfo-alpha-L-idopyranuronic acid","O5 C1 C2 C3 C4 C5", 0.620, "1c4", 0.003, 2.016 },
    { "IMK" ,	 "B", 	 "D", 	 "2-(BETA-D-GLUCOPYRANOSYL)-5-METHYL-1-BENZIMID" ,  "O5 C1 C2 C3 C4 C5", 0.555, "4c1", 0.002, 1.855 },
    { "IPT" ,	 "B", 	 "D", 	 "ISOPROPYL-1-BETA-D-THIOGALACTOSIDE" ,             "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.042 },
    { "IXD" ,	 "N", 	 "N", 	 "4-DEOXY-2-O-SULFO-BETA-D-ERYTHRO-HEX-4-ENOPYR" ,  "O5 C1 C2 C3 C4 C5", 0.511, "2h1", 0.095, 7.412 },
    { "JHM" ,	 "A", 	 "D", 	 "2-DEOXY-6-O-SULFO-ALPHA-D-ARABINO-HEXOPYRANOS" ,  "O5 C1 C2 C3 C4 C5", 0.592, "4c1", 0.001, 1.370 },
    { "JZR" ,	 "B", 	 "D", 	 "HEXYLBETA-D-GLUCOPYRANOSIDE" ,                    "O5 C1 C2 C3 C4 C5", 0.594, "4c1", 0.001, 1.393 },
    { "K5B" ,	 "N", 	 "D", 	 "4,7-anhydro-3-deoxy-D-gluco-oct-2-ulosonic acid" ,"O7 C4 C5 C6 C7", 0.621, "2ev", 0.003, 2.048 },
    { "KDA" ,	 "A", 	 "D", 	 "(3-DEOXY-D-MANNO-OCT-2-ULOSONICACID)-2-O-ALL" ,   "O6 C2 C3 C4 C5 C6", 0.622, "4c1", 0.003, 2.027 },
    { "KD5" ,	 "N", 	 "D", 	 "4,7-anhydro-3-deoxy-D-manno-oct-2-ulosonic acid", "O7 C4 C5 C6 C7", 0.500, "Oh5", 0.095, 8.063 },
    { "KDB" ,	 "B", 	 "D", 	 "3,4,5-TRIDEOXY-ALPHA-D-ERYTHRO-OCT-3-EN-2-ULO" ,  "O6 C2 C3 C4 C5 C6", 0.500, "Oh5", 0.095, 8.063 },
    { "KDE" ,	 "B", 	 "L", 	 "PROP-2-EN-1-YL3-DEOXY-BETA-L-GULO-OCT-2-ULOP" ,   "O6 C2 C3 C4 C5 C6", 0.592, "4c1", 0.001, 1.366 },
    { "KDO" ,	 "B", 	 "D", 	 "3-DEOXY-D-MANNO-OCT-2-ULOSONICACID" ,             "O6 C2 C3 C4 C5 C6", 0.621, "4c1", 0.003, 2.048 },
    { "KDR" ,	 "B", 	 "D", 	 "PROP-2-EN-1-YL3-DEOXY-ALPHA-D-MANNO-OCTOS-2-" ,   "O6 C2 C3 C4 C5 C6", 0.592, "4c1", 0.001, 1.348 },
    { "KHP" ,	 "A", 	 "L", 	 "2-HYDROXYMETHYL-5-(4-NITRO-PHENOXY)-TETRAH" ,   "O4' C1B C2B C3B C4B", 0.379, "ev3", 0.013, 5.447 },
    { "KME" ,	 "A", 	 "D", 	 "(1E)-PROP-1-EN-1-YL3-DEOXY-7-O-METHYL-ALPHA-" ,   "O6 C2 C3 C4 C5 C6", 0.593, "4c1", 0.001, 1.348 },
    { "KO1" ,	 "A", 	 "D", 	 "D-GLYCERO-ALPHA-D-TALO-OCT-2-ULOPYRANOSONICA" ,   "O6 C2 C3 C4 C5 C6", 0.593, "4c1", 0.001, 1.379 },
    { "KO2" ,	 "A", 	 "D", 	 "PROP-2-EN-1-YLD-GLYCERO-ALPHA-D-TALO-OCT-2-U" ,   "O6 C2 C3 C4 C5 C6", 0.593, "4c1", 0.002, 1.365 },
    { "KOT" ,	 "B", 	 "D", 	 "1-BETA-D-GLUCOPYRANOSYL-4-PHENYL-1H-1,2,3-TRI" ,  "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.349 },
    { "L6S" ,	 "A", 	 "L", 	 "6-O-SULFO-ALPHA-L-GALACTOSE" ,                    "O5 C1 C2 C3 C4 C5", 0.593, "1c4", 0.002, 1.407 },
    { "LDY" ,	 "A", 	 "D", 	 "ALPHA-D-LYXOPYRANOSE" ,                           "O5 C1 C2 C3 C4 C5", 0.594, "1c4", 0.001, 1.370 },
    { "LFR" ,	 "B", 	 "L", 	 "BETA-L-FRUCTOFURANOSE" ,                          "O5 C2 C3 C4 C5",    0.372, "ev3", 0.016, 5.700 },
    { "LGC" ,	 "A", 	 "D", 	 "(3S,4R,5R,6S)-3,4,5-TRIHYDROXY-6-(HYDROXYMETH" ,  "O5 C1 C2 C3 C4 C5", 0.537, "4c1", 0.044, 5.611 },
    { "LGU" ,	 "A", 	 "L", 	 "ALPHA-L-GULURONATE" ,                             "O5 C1 C2 C3 C4 C5", 0.593, "1c4", 0.001, 1.375 },
    { "LOG" ,	 "A", 	 "D", 	 "LOGNAC" ,                                         "O5 C1 C2 C3 C4 C5", 0.550, "1c4", 0.001, 1.375 },
    { "LRH" ,	 "B", 	 "L", 	 "6-DEOXY-BETA-L-FRUCTOFURANOSE" ,                  "O5 C2 C3 C4 C5",    0.313, "1tO", 0.010, 4.410 },
    { "LXB" ,	 "B", 	 "D", 	 "2-(ACETYLAMINO)-2-DEOXY-BETA-D-GULOPYRANOSE" ,    "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.415 },
    { "LXC" ,	 "B", 	 "L", 	 "BETA-L-XYLOPYRANOSE" ,                            "O5 C1 C2 C3 C4 C5", 0.619, "1c4", 0.003, 2.003 },
    { "LXZ" ,	 "A", 	 "D", 	 "2-(ACETYLAMINO)-2-DEOXY-ALPHA-D-IDOPYRANOSE" ,    "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.384 },
    { "LZ0" ,	 "A", 	 "L", 	 "[1-(2-OXOETHYL)-1H-1,2,3-TRIAZOL-5-YL]METHYL" ,   "O5 C1 C2 C3 C4 C5", 0.537, "1c4", 0.011, 2.592 },
    { "M07" ,	 "N", 	 "D", 	 "(5R,7R,8S,9S,10R)-7-(HYDROXYMETHYL)-3-(4-METH" ,  "O5 C1 C2 C3 C4 C5", 0.583, "4c1", 0.003, 1.678 },
    { "M08" ,	 "N", 	 "D", 	 "(5R,7R,8S,9S,10R)-7-(HYDROXYMETHYL)-3-PHENYL-" ,  "O5 C1 C2 C3 C4 C5", 0.583, "4c1", 0.004, 1.695 },
    { "M09" ,	 "N", 	 "D", 	 "(3S,5R,7R,8S,9S,10R)-7-(HYDROXYMETHYL)-3-(4-N" ,  "O5 C1 C2 C3 C4 C5", 0.586, "4c1", 0.002, 1.584 },
    { "M1P" ,	 "A", 	 "D", 	 "ALPHA-D-MANNOSE1-PHOSPHATE" ,                     "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.002, 2.060 },
    { "M6D" ,	 "B", 	 "D", 	 "6-O-PHOSPHONO-BETA-D-MANNOPYRANOSE" ,             "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.390 },
    { "M6P" ,	 "A", 	 "D", 	 "ALPHA-D-MANNOSE-6-PHOSPHATE" ,                    "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.027 },
    { "M7P" ,	 "B", 	 "D", 	 "D-GLYCERO-D-MANNOPYRANOSE-7-PHOSPHATE" ,          "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.028 },
    { "M8C" ,	 "A", 	 "D", 	 "METHYLALPHA-D-GALACTOPYRANURONATE" ,              "O5 C1 C2 C3 C4 C5", 0.564, "4c1", 0.003, 1.223 },
    { "MA1" ,	 "A", 	 "D", 	 "1,4-DITHIO-ALPHA-D-GLUCOPYRANOSE" ,               "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.049 },
    { "MA2" ,	 "A", 	 "D", 	 "4-S-METHYL-4-THIO-ALPHA-D-GLUCOPYRANOSE" ,        "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.035 },
    { "MA3" ,	 "A", 	 "D", 	 "O1-METHYL-4-DEOXY-4-THIO-ALPHA-D-GLUCOSE" ,       "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.024 },
    { "MAG" ,	 "B", 	 "D", 	 "BETA-METHYL-N-ACETYL-D-GLUCOSAMINE" ,             "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.020 },
    { "MAN" ,	 "A", 	 "D", 	 "ALPHA-D-MANNOSE" ,                                "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.039 },
    { "MAT" ,	 "A", 	 "L", 	 "2,4-DIDEOXY-4-[2-(PROPYL)AMINO]-3-O-METHYLAL" ,   "O5 C1 C2 C3 C4 C5", 0.592, "1c4", 0.002, 1.377 },
    { "MAV" ,	 "A", 	 "D", 	 "ALPHA-D-MANNOPYRANURONICACID" ,                   "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.393 },
    { "MAW" ,	 "A", 	 "L", 	 "4-DEOXY-ALPHA-L-ERYTHRO-HEX-4-ENOPYRANURONIC" ,   "O5 C1 C2 C3 C4 C5", 0.512, "2h1", 0.095, 7.411 },
    { "MBF" ,	 "B", 	 "D", 	 "2-DEOXY-2-FLUORO-BETA-D-MANNOSE" ,                "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.004, 2.024 },
    { "MBG" ,	 "B", 	 "D", 	 "METHYL-BETA-GALACTOSE" ,                          "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.026 },
    { "MDA" ,	 "B", 	 "D", 	 "2,6-DIDEOXY-3C-METHYL-D-RIBOPYRANOSIDE" ,         "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.033 },
    { "MDP" ,	 "B", 	 "D", 	 "N-CARBOXY-N-METHYL-MURAMICACID" ,                 "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.380 },
    { "MFB" ,	 "B", 	 "L", 	 "BETA-L-METHYL-FUCOSE" ,                           "O5 C1 C2 C3 C4 C5", 0.564, "1c4", 0.002, 1.726 },
    { "MFU" ,	 "A", 	 "L", 	 "ALPHA-L-METHYL-FUCOSE" ,                          "O5 C1 C2 C3 C4 C5", 0.620, "1c4", 0.003, 2.024 },
    { "MGC" ,	 "A", 	 "D", 	 "ALPHA-METHYL-N-ACETYL-D-GALACTOSAMINE" ,          "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.031 },
    { "MGL" ,	 "B", 	 "D", 	 "O1-METHYL-GLUCOSE" ,                              "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.016 },
    { "MGS" ,	 "A", 	 "D", 	 "1,2-O-DIMETHYL-4-[2,4-DIHYDROXY-BUTYRAMIDO]-4" ,  "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.038 },
    { "MMA" ,	 "A", 	 "D", 	 "O1-METHYL-MANNOSE" ,                              "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.009 },
    { "MN0" ,	 "A", 	 "L", 	 "METHYL3,5-DIDEOXY-5-[(HYDROXYACETYL)AMINO]-D" ,   "O6 C2 C3 C4 C5 C6", 0.593, "1c4", 0.001, 1.376 },
    { "MNA" ,	 "A", 	 "L", 	 "2-O-METHYL-5-N-ACETYL-ALPHA-D-NEURAMINICACI" ,    "O6 C2 C3 C4 C5 C6", 0.620, "1c4", 0.003, 2.025 },
    { "MQT" ,	 "B", 	 "D", 	 "METHYL2-O-ACETYL-3-O-(4-METHYLBENZOYL)-BETA-" ,   "O5 C1 C2 C3 C4 C5", 0.594, "4c1", 0.001, 1.352 },
    { "MRP" ,	 "A", 	 "L", 	 "3-O-METHYL-ALPHA-L-RHAMNOPYRANOSIDE" ,            "O5 C1 C2 C3 C4 C5", 0.619, "1c4", 0.003, 2.019 },
    { "MUB" ,	 "A", 	 "D", 	 "N-ACETYLMURAMICACID" ,                            "O5 C1 C2 C3 C4 C5", 0.518, "1c4", 0.005, 2.794 },
    { "MUR" ,	 "B", 	 "D", 	 "MURAMICACID" ,                                    "O5 C1 C2 C3 C4 C5", 0.619, "4c1", 0.003, 1.981 },
    { "MXY" ,	 "B", 	 "L", 	 "2-O-methyl-beta-L-fucopyranose" ,                 "O5 C1 C2 C3 C4 C5", 0.621, "1c4", 0.003, 2.042 },
    { "MXZ" ,	 "A", 	 "L", 	 "6-DEOXY-2-O-METHYL-ALPHA-L-GALACTOPYRANOSE" ,     "O5 C1 C2 C3 C4 C5", 0.621, "1c4", 0.003, 2.042 },
    { "N1L" ,	 "B", 	 "D", 	 "2-AMINO-2-DEOXY-BETA-D-GLUCOPYRANURONICACID" ,    "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.389 },
    { "NA1" ,	 "A", 	 "D", 	 "METHYLN-ACETYLALLOSAMINE" ,                       "O5 C1 C2 C3 C4 C5", 0.537, "1c4", 0.004, 2.320 },
    { "NAA" ,	 "B", 	 "D", 	 "N-ACETYL-D-ALLOSAMINE" ,                          "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.014 },
    { "NAG" ,	 "B", 	 "D", 	 "N-ACETYL-D-GLUCOSAMINE" ,                         "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.391 },
    { "NBG" ,	 "B", 	 "D", 	 "1-N-ACETYL-BETA-D-GLUCOSAMINE" ,                  "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.004, 1.995 },
    { "NDG" ,	 "A", 	 "D", 	 "2-(ACETYLAMINO)-2-DEOXY-A-D-GLUCOPYRANOSE" ,      "O5 C1 C2 C3 C4 C5" , 0.593, "4c1", 0.001, 1.378 },
    { "NFG" ,	 "B", 	 "D", 	 "2,4-dinitrophenyl 2-deoxy-2-fluoro-beta-D-glucopyranoside", "O5 C1 C2 C3 C4 C5" , 0.593, "4c1", 0.001, 1.378 }, // new entry 01/12/2021
    { "NGS" ,	 "B", 	 "D", 	 "2-acetamido-2-deoxy-6-O-sulfo-beta-D-glucopyranose", "O5 C1 C2 C3 C4 C5", 0.540, "4c1", 0.010, 1.500 }, // new entry 01/12/2021
    { "NG1" ,	 "A", 	 "D", 	 "N-ACETYL-ALPHA-D-GALACTOSAMINE1-PHOSPHATE" ,      "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.002, 2.010 },
    { "NG6" ,	 "B", 	 "D", 	 "N-ACETYL-D-GALACTOSAMINE 6-SULFATE" ,             "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.002, 2.010 },
    { "NGA" ,	 "B", 	 "D", 	 "N-ACETYL-D-GALACTOSAMINE" , 	 		            "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.347 },
    { "NGC" ,	 "A", 	 "D", 	 "N-glycolyl-alpha-neuraminic acid" , 	 		    "O6 C2 C3 C4 C5 C6", 0.593, "4c1", 0.002, 1.347 },
    { "NGK" ,	 "A", 	 "D", 	 "2-(ACETYLAMINO)-2-DEOXY-4-O-SULFO-ALPHA-D-GALAC" ,"O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.347 },
    { "NGZ" ,	 "A", 	 "L", 	 "2-(ACETYLAMINO)-2-DEOXY-ALPHA-L-GLUCOPYRANOSE" ,  "O5 C1 C2 C3 C4 C5", 0.620, "1c4", 0.003, 2.015 },
    { "NNG" ,	 "A", 	 "D", 	 "2-DEOXY-2-{[(S)-HYDROXY(METHYL)PHOSPHORYL]AMI" ,  "O5 C1 C2 C3 C4 C5", 0.557, "4c1", 0.002, 1.599 },
    { "NOK" ,	 "N", 	 "D", 	 "2-ACETAMIDO-1,2-DIDEOXYNOJIRMYCIN" ,              "N5 C1 C2 C3 C4 C5", 0.557, "4c1", 0.010, 2.599 },
    { "NTF" ,	 "B", 	 "D", 	 "N-TRIFLURO-ACETYL-BETA-D-GLUCOPYRANOSYLAMINE" ,   "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.002, 2.084 },
    { "NXD" ,	 "B", 	 "L", 	 "METHYL5-(ACETYLAMINO)-9-{[AMINO(OXO)ACETYL]A" ,   "O6 C2 C3 C4 C5 C6", 0.621, "1c4", 0.002, 2.064 },
    { "OAK" ,	 "B", 	 "D", 	 "N-(PHENYLCARBONYL)-BETA-D-GLUCOPYRANOSYLAMINE" ,  "O5 C1 C2 C3 C4 C5", 0.592, "4c1", 0.001, 1.368 },
    { "OI7" ,	 "A", 	 "D", 	 "1,7-DI-O-PHOSPHONO-BETA-D-ALTRO-HEPT-2-ULOFUR" ,  "O5 C2 C3 C4 C5"   , 0.379, "ev4", 0.016, 5.864 },
    { "OPM" ,	 "A", 	 "D", 	 "O1-PENTYL-MANNOSE" ,                              "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.023 },
    { "OTG" ,	 "A", 	 "D", 	 "ORTHO-TOLUOYLGLUCOSAMINE" ,                       "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.032 },
    { "OX2" ,	 "B", 	 "D", 	 "2-(BETA-D-GLUCOPYRANOSYL)-5-METHYL-1,3,4-OXAD" ,  "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.014 },
    // { "OXZ" ,	 "N", 	 "N", 	 "TETRAHYDROOXAZINE" ,                              "O5 N1 C2 C3 C4 C5", 0.622, "4c1", 0.001, 2.060 }, // is it worth validating ring monomers like this
    { "P6P" ,	 "A", 	 "D", 	 "6-O-PHOSPHONO-ALPHA-D-FRUCTOFURANOSE" ,           "O5 C2 C3 C4 C5"   , 0.381, "Oev", 0.015, 5.689 },
    { "PA1" ,	 "A", 	 "D", 	 "2-AMINO-2-DEOXY-ALPHA-D-GLUCOPYRANOSE" ,          "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.023 },
    { "PDX" ,	 "A", 	 "D", 	 "2,3-DI-O-SULFO-ALPHA-D-GLUCOPYRANOSE" ,           "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 2.037 },
    { "PH5" ,	 "A", 	 "L", 	 "2-PHENYL-PROP5AC" ,                               "O6 C2 C3 C4 C5 C6", 0.621, "1c4", 0.003, 2.005 },
    { "PKM" ,	 "A", 	 "D", 	 "4-O-acetyl-5-acetamido-3,5-dideoxy-D-glycero-alpha-D-galacto-non-2-ulopyranosonic acid" , "O6 C2 C3 C4 C5 C6", 0.593, "1c4", 0.001, 1.379 },
    { "PNA" ,	 "A", 	 "D", 	 "4'-NITROPHENYL-ALPHA-D-MANNOPYRANOSIDE" ,         "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.027 },
    { "PNG" ,	 "A", 	 "D", 	 "4'-NITROPHENYL-ALPHA-D-GLUCOPYRANOSIDE" ,         "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.002, 2.063 },
    { "PNJ" ,	 "B", 	 "D", 	 "PNP-BETA-D-GLUCOSAMINE" ,                         "O5 C1 C2 C3 C4 C5", 0.618, "4c1", 0.018, 4.026 },
    { "PNW" ,	 "B", 	 "D", 	 "4-NITROPHENYLBETA-D-GLUCOPYRANOSIDE" ,            "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.339 },
    { "PPC" ,	 "A", 	 "D", 	 "5-PHOSPHORIBOSYL-1-(BETA-METHYLENE)PYROPHOSP" ,   "O4 C1 C2 C3 C4"   , 0.381, "Oev", 0.015, 5.688 },
    { "PRP" ,	 "A", 	 "D", 	 "ALPHA-PHOSPHORIBOSYLPYROPHOSPHORICACID" ,         "O4 C1 C2 C3 C4"   , 0.381, "Oev", 0.015, 5.696 },
    { "PSG" ,	 "B", 	 "D", 	 "PARA-NITROPHENYL1-THIO-BETA-D-GLUCOPYRANOSID" ,   "O5 C1 C2 C3 C4 C5", 0.558, "4c1", 0.002, 1.830 },
    { "PSV" ,	 "A", 	 "D", 	 "ALPHA-D-PSICOFURANOSE" ,                          "O5 C2 C3 C4 C5"   , 0.378, "ev4", 0.017, 5.851 },
    { "RAM" ,	 "A", 	 "L", 	 "ALPHA-L-RHAMNOSE" ,                               "O5 C1 C2 C3 C4 C5", 0.620, "1c4", 0.003, 2.041 },
    { "RAO" ,	 "A", 	 "L", 	 "1-O-METHYL-ALPHA-RHAMNOSE" ,                      "O5 C1 C2 C3 C4 C5", 0.621, "1c4", 0.002, 2.046 },
    { "RER" ,	 "A", 	 "L", 	 "(1R,3S,4S,5S)-3-AMINO-2,3,6-TRIDEOXY-3-METHYL" ,  "O5 C1 C2 C3 C4 C5", 0.592, "1c4", 0.002, 1.366 },
    { "RGG" ,	 "B", 	 "D", 	 "(2R)-2,3-DIHYDROXYPROPYLBETA-D-GALACTOPYRANO" ,   "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.378 },
    { "RHC" ,	 "B", 	 "D", 	 "5-(3-AMINO-4,4-DIHYROXY-BUTYLSULFANYLMETHYL)-" ,  "O4 C1 C2 C3 C4"   , 0.380, "ev1", 0.014, 5.459 },
    { "RIB" ,	 "A", 	 "D", 	 "alpha-D-ribofuranose" ,                           "O4 C1 C2 C3 C4"   , 0.372, "3ev", 0.015, 5.717 },
    { "RI2" ,	 "A", 	 "D", 	 "1,5-DI-O-PHOSPHONO-ALPHA-D-RIBOFURANOSE" ,        "O4 C1 C2 C3 C4"   , 0.372, "3ev", 0.015, 5.717 },
    { "RIP" ,	 "B", 	 "D", 	 "RIBOSE(PYRANOSEFORM)" ,                           "O5 C1 C2 C3 C4 C5", 0.619, "4c1", 0.004, 1.990 },
    { "RR7" ,	 "B", 	 "D", 	 "2-deoxy-beta-D-arabino-hexopyranose" ,            "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 2.051 }, // new entry 01/12/2021
    { "RY7" ,	 "A", 	 "D", 	 "4,6-dideoxy-4-{[(1S,2S,3S,4R,5R)-2,3,4-trihydroxy-5-(hydroxymethyl)cyclohexyl]amino}-alpha-D-glucopyranose", "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.023 },
    { "RM4" ,	 "B", 	 "L", 	 "6-DEOXY-BETA-L-MANNOPYRANOSE" ,                   "O5 C1 C2 C3 C4 C5", 0.593, "1c4", 0.001, 1.370 },
    { "RST" ,	 "A", 	 "L", 	 "RISTOSAMINE" ,                                    "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.388 },
    { "RTV" ,	 "N", 	 "D", 	 "2-(acetylamino)-1,5-anhydro-2-deoxy-D-mannitol" , "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.388 },
    { "RUG" ,	 "B", 	 "D", 	 "1-BETA-D-GLUCOPYRANOSYL-4-(HYDROXYMETHYL)-1H-" ,  "O5 C1 C2 C3 C4 C5", 0.594, "4c1", 0.002, 1.365 },
    { "S06" ,	 "N", 	 "D", 	 "(3S,5R,7R,8S,9S,10R)-7-(HYDROXYMETHYL)-3-(2-N" ,  "O5 C1 C2 C3 C4 C5", 0.585, "4c1", 0.003, 1.603 },
    { "S13" ,	 "N", 	 "D", 	 "(3S,5R,7R,8S,9S,10R)-7-(HYDROXYMETHYL)-3-(4-M" ,  "O5 C1 C2 C3 C4 C5", 0.585, "4c1", 0.002, 1.588 },
    { "SFU" ,	 "A", 	 "L", 	 "METHYL6-DEOXY-1-SELENO-ALPHA-L-GALACTOPYRANO" ,   "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.365 },
    { "SG4" ,	 "A", 	 "D", 	 "3,4-DI-O-ACETYL-6-O-SULFAMOYL-ALPHA-D-GLUCOPY" ,  "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.397 },
    { "SGA" ,	 "B", 	 "D", 	 "O3-SULFONYLGALACTOSE" ,                           "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.039 },
    { "SGC" ,	 "B", 	 "D", 	 "4-DEOXY-4-THIO-BETA-D-GLUCOPYRANOSE" ,            "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.022 },
    { "SGN" ,	 "A", 	 "D", 	 "N,O6-DISULFO-GLUCOSAMINE" ,                       "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 1.997 },
    { "SHB" ,	 "B", 	 "D", 	 "methyl beta-D-galactopyranuronate" ,              "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.370 },
    { "SHD" ,	 "A", 	 "D", 	 "alpha-D-altropyranose" ,                          "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.051 },
    { "SHG" ,	 "B", 	 "D", 	 "2-DEOXY-2-FLUORO-BETA-D-GLUCOPYRANOSE" ,          "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.370 },
    { "SIA" ,	 "A", 	 "L", 	 "O-SIALIC ACID" ,                                  "O6 C2 C3 C4 C5 C6", 0.621, "1c4", 0.003, 1.996 },
    { "SID" ,	 "A", 	 "L", 	 "METHYL9-S-ACETYL-5-(ACETYLAMINO)-3,5-DIDEOXY" ,   "O6 C2 C3 C4 C5 C6", 0.593, "1c4", 0.001, 1.347 },
    { "SLB" ,	 "B", 	 "L", 	 "5-N-ACETYL-BETA-D-NEURAMINICACID" ,               "O6 C2 C3 C4 C5 C6", 0.620, "1c4", 0.003, 2.001 },
    { "SLM" ,	 "B", 	 "L", 	 "SIALYLAMIDE" ,                                    "O6 C2 C3 C4 C5 C6", 0.592, "1c4", 0.001, 1.375 },
    { "SN5" ,	 "B", 	 "D", 	 "2-DEOXY-2-(ETHANETHIOYLAMINO)-BETA-D-GLUCOPYR" ,  "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.372 },
    { "SNG" ,	 "B", 	 "D", 	 "METHYL2-ACETAMIDO-1,2-DIDEOXY-1-SELENO-BETA-" ,   "O5 C1 C2 C3 C4 C5", 0.606, "4c1", 0.009, 2.290 },
    { "SOG" ,    "B",    "D",    "2-HYDROXYMETHYL-6-OCTYLSULFANYL-TETRAHYDRO-PYRAN","O5 C1 C2 C3 C4 C5", 0.606, "4c1", 0.003, 2.090 },
    { "SOE" ,	 "A", 	 "L", 	 "ALPHA-L-SORBOPYRANOSE" ,                          "O6 C2 C3 C4 C5 C6", 0.593, "1c4", 0.001, 1.366 },
    { "SSG" ,	 "B", 	 "D", 	 "1,4-DEOXY-1,4-DITHIO-BETA-D-GLUCOPYRANOSE" ,      "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.049 },
    { "STZ" ,	 "B", 	 "D", 	 "STREPTOZOTOCIN" ,                                 "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.369 },
    { "SUS" ,	 "A", 	 "D", 	 "2-DEOXY-3,6-DI-O-SULFO-2-(SULFOAMINO)-ALPHA-D" ,  "O5 C1 C2 C3 C4 C5", 0.584, "4c1", 0.004, 2.350 },
    { "V3P" ,	 "B", 	 "D", 	 "4-iodophenyl 1,4-dithio-beta-D-glucopyranoside" , "O5 C1 C2 C3 C4 C5", 0.594, "4c1", 0.001, 1.393 },
    { "TA6" ,	 "B", 	 "D", 	 "6-O-PHOSPHONO-BETA-D-TAGATOFURANOSE" ,            "O5 C2 C3 C4 C5"   , 0.372, "3ev", 0.016, 5.716 },
    { "TGA" ,	 "B", 	 "D", 	 "METHANETHIOSULFONYL-GALACTOSIDE" ,                "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.368 },
    { "TM9" ,    "B",    "D",    "WRONG INTERPRETATION OF GCS (GLUCOSAMINE)",       "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.368 },
    { "TMR" ,	 "B", 	 "D", 	 "2,6-DIDEOXY-4-THIOMETHYL-BETA-D-RIBOHEXOPYRAN" ,  "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.020 },
    { "TMX" ,	 "B", 	 "D", 	 "2-DEOXY-2-(TRIMETHYLAMMONIO)-BETA-D-GLUCOPYRA" ,  "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.370 },
    { "TNR" ,	 "A", 	 "D", 	 "O-(2-ACETAMIDO-2-DEOXY-ALPHA-D-GALACTOPYRANOSE" , "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.001 },
    { "TRV" ,	 "B", 	 "D", 	 "6-O-octanoyl-beta-D-fructofuranose" ,             "O5 C2 C3 C4 C5"   , 0.372, "3ev", 0.015, 5.727 },
    { "TOA" ,	 "A", 	 "D", 	 "3-DEOXY-3-AMINOGLUCOSE" ,                         "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.002, 2.021 },
    { "TOC" ,	 "A", 	 "D", 	 "2,3,6-TRIDEOXY-2,6-DIAMINOGLUCOSE" , 	            "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 2.047 },
    { "TQY" ,	 "A", 	 "D", 	 "6-O-octanoyl-alpha-D-glucopyranose" , 	        "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 2.047 },
    { "TT7" ,	 "B", 	 "D", 	 "4-O-phosphono-beta-D-fructofuranose" ,            "O5 C2 C3 C4 C5"   , 0.372, "3ev", 0.016, 5.716 }, // new entry 01/12/2021
    { "TUP" ,	 "A", 	 "D", 	 "3-deoxy-3-fluoro-alpha-D-glucopyranose" ,         "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.047 },
    { "TVG" ,	 "B", 	 "D", 	 "propyl beta-D-galactofuranoside" ,                "O5 C1 C2 C3 C4 C5", 0.540, "4c1", 0.010, 1.500 }, // new entry 01/12/2021
    { "TVS" ,	 "B", 	 "D", 	 "prop-2-en-1-yl 2-(acetylamino)-2-deoxy-beta-D-glucopyranoside" , "O5 C1 C2 C3 C4 C5", 0.540, "4c1", 0.010, 1.500 }, // new entry 01/12/2021
    { "TVV" ,	 "B", 	 "D", 	 "	3-O-prop-2-yn-1-yl-beta-D-galactopyranose" ,    "O5 C1 C2 C3 C4 C5", 0.540, "4c1", 0.010, 1.500 }, // new entry 01/12/2021
    { "TWA" ,	 "B", 	 "D", 	 "2,3,4-tri-O-sulfo-beta-D-altropyranose" ,         "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.051 },
    { "TWD" ,	 "B", 	 "D", 	 "2,3-di-O-sulfo-alpha-L-glucopyranose" ,           "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.051 },
    { "TXB" ,	 "A", 	 "D", 	 "4-deoxy-4-thio-alpha-D-xylopyranose" ,            "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.002, 2.039 },
    { "TYV" ,	 "A", 	 "D", 	 "TYVELOSE" ,                                       "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.047 },
    { "U1Y" ,	 "B", 	 "D", 	 "methyl 6-thio-beta-D-glucopyranoside",            "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 2.051 }, // new entry 01/12/2021
    { "U2A" ,	 "B", 	 "D", 	 "methyl 2-thio-beta-D-glucopyranoside",            "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 2.051 }, // new entry 01/12/2021
    { "U2D" ,	 "A", 	 "D", 	 "6-O-decanoyl-alpha-D-glucopyranose",              "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 2.051 }, // new entry 01/12/2021
    { "UAP" ,	 "A", 	 "D", 	 "1,4-DIDEOXY-5-DEHYDRO-O2-SULFO-GLUCURONICACI" ,   "O5 C1 C2 C3 C4 C5", 0.442, "1h2", 0.100, 8.090 },
    { "VG1" ,	 "A", 	 "D", 	 "ALPHA-D-GLUCOSE-1-PHOSPHATE-6-VANADATE" ,         "O5 C1 C2 C3 C4 C5", 0.619, "4c1", 0.003, 2.033 },
    { "VJ1" ,	 "A", 	 "D", 	 "2-deoxy-3-O-[(1R,3R)-1,3-dihydroxytetradecyl]-2-{[(3R)-3-hydroxytetradecanoyl]amino}-1-O-phosphono-alpha-D-glucopyranose" , "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.045 },
    { "VJ4" ,	 "A", 	 "D", 	 "2-deoxy-2-{[(1S,3R)-1-hydroxy-3-(pentanoyloxy)undecyl]amino}-4-O-phosphono-alpha-D-glucopyranose" , "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.045 },
    { "WIA" ,	 "B", 	 "D", 	 "methyl 6-thio-beta-D-galactopyranoside" ,         "O5 C1 C2 C3 C4 C5", 0.623, "4c1", 0.003, 2.051 },
    { "X1P" ,	 "A", 	 "D", 	 "1-O-PHOSPHONO-ALPHA-D-XYLOPYRANOSE" ,             "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.002, 2.039 },
    { "X2F" ,	 "A", 	 "D", 	 "2-DEOXY-2-FLUOROXYLOPYRANOSE" ,                   "O5 C1 C2 C3 C4 C5", 0.620, "4c1", 0.003, 2.070 },
    { "YIO" ,	 "B", 	 "D", 	 "1-thio-beta-D-galactopyranose" ,                  "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.003, 2.051 },
    { "XYF" ,	 "B", 	 "D", 	 "5(R)-5-FLUORO-BETA-D-XYLOPYRANOSYL-ENZYMEINT" ,   "O5 C1 C2 C3 C4 C5", 0.557, "4c1", 0.012, 1.634 },
    { "XYP" ,	 "B",	 "D",	 "BETA-D-XYLOPYRANOSE" ,                            "O5 C1 C2 C3 C4 C5", 0.619, "4c1", 0.003, 2.003 },
    { "YYQ" ,	 "A", 	 "L", 	 "2-acetamido-2-deoxy-alpha-L-galactopyranose" , 	"O5 C1 C2 C3 C4 C5", 0.594, "4c1", 0.002, 1.388 },
    { "XYS" ,	 "A", 	 "D", 	 "XYLOPYRANOSE" ,                                   "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 2.051 },
    { "XXR" ,	 "A", 	 "D", 	 "alpha-D-rhamnopyranose" ,                         "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 2.051 },
    { "YYJ" ,	 "B", 	 "D", 	 "1,3,4,6-tetra-O-sulfo-beta-D-fructofuranose" ,    "O5 C2 C3 C4 C5"   , 0.372, "3ev", 0.016, 5.716 }, // new entry 01/12/2021
    { "YX0" ,	 "A", 	 "L", 	 "[(3E)-3-(1-HYDROXYETHYLIDENE)-2,3-DIHYDROISOX" ,  "O5 C1 C2 C3 C4 C5", 0.594, "1c4", 0.002, 1.379 },
    { "YX1" ,	 "B", 	 "D", 	 "2-DEOXY-2-{[(2-HYDROXY-1-METHYLHYDRAZINO)CARB" ,  "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.001, 1.385 },
    { "YZ0" ,	 "B", 	 "D", 	 "methyl beta-D-mannopyranoside" ,                  "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.001, 2.060 },
    { "Z3Q" ,	 "B", 	 "D", 	 "2-azidoethyl 2-acetamido-2-deoxy-beta-D-glucopyranoside", "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 2.051 }, // new entry 01/12/2021
    { "Z4R" ,	 "A", 	 "D", 	 "methyl 3-thio-alpha-D-mannopyranoside" ,          "O5 C1 C2 C3 C4 C5", 0.622, "4c1", 0.001, 2.060 },
    { "Z4Y" ,	 "A", 	 "D", 	 "6-thio-alpha-D-mannopyranose" ,                   "O5 C1 C2 C3 C4 C5", 0.548, "4c1", 0.002, 1.803 },
    { "Z5L" ,	 "A", 	 "D", 	 "methyl 2-thio-alpha-D-mannopyranoside" ,          "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.003, 2.027 },
    { "Z61" ,	 "A", 	 "D", 	 "2-deoxy-alpha-D-arabino-hexopyranose",            "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 2.051 }, // new entry 01/12/2021
    { "Z6W" ,	 "B", 	 "D", 	 "(5R)-5-[(2R)-2-hydroxynonyl]-beta-D-xylulofuranose", "O5 C2 C3 C4 C5", 0.372, "ev2", 0.002, 2.051 }, // new entry 01/12/2021
    { "Z9D" ,	 "B", 	 "D", 	 "4-deoxy-beta-D-xylo-hexopyranose" ,               "O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 2.051 }, // new entry 01/12/2021
    { "Z9H" ,	 "A", 	 "D", 	 "3,4-di-O-methyl-2,6-di-O-sulfo-alpha-D-glucopyranose","O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 2.051 }, // new entry 01/12/2021
    { "Z9K" ,	 "A", 	 "L", 	 "3-O-methyl-2-O-sulfo-alpha-L-idopyranuronic acid","O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 2.051 }, // new entry 01/12/2021
    { "Z9L" ,	 "A", 	 "D", 	 "methyl 2,3,6-tri-O-sulfo-alpha-D-glucopyranoside","O5 C1 C2 C3 C4 C5", 0.621, "4c1", 0.002, 2.051 }, // new entry 01/12/2021
    { "Z9M" ,	 "B", 	 "D", 	 "2-amino-2-deoxy-4-O-phosphono-beta-D-glucopyranose" ,"O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.374 },
    { "Z9N" ,	 "A", 	 "D", 	 "alpha-D-fructofuranose" ,                         "O5 C2 C3 C4 C5", 0.372, "3ev", 0.002, 2.051 }, // new entry 01/12/2021
    { "ZD0" ,	 "A", 	 "D", 	 "methyl 4,6-dideoxy-4-{[(2R)-2,4-dihydroxybutanoyl]amino}-alpha-D-mannopyranoside", "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.377 },
    { "ZDO" ,	 "A", 	 "D", 	 "methyl 2-deoxy-6-O-sulfo-2-(sulfoamino)-alpha-D-glucopyranoside", "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.377 },
    { "ZCZ" ,	 "A", 	 "D", 	 "4,6-dideoxy-4-{[(2R)-2,4-dihydroxybutanoyl]amino}-2-O-methyl-alpha-D-mannopyranose", "O5 C1 C2 C3 C4 C5", 0.593, "4c1", 0.002, 1.377 },
    { "  A" ,	 "B",	 "D", 	 "ADENOSINE-5'-MONOPHOSPHATE" ,                   "O4' C1' C2' C3' C4'", 0.236, "Oev", 0.016, 3.834 },   // RNA/DNA
    { "  U" ,    "B",    "D",    "URIDINE-5'-MONOPHOSPHATE" ,                     "O4' C1' C2' C3' C4'", 0.291, "Oev", 0.019, 4.394 },   //
    { "  G" ,    "B",    "D",    "GUANOSINE-5'-MONOPHOSPHATE" ,                   "O4' C1' C2' C3' C4'", 0.241, "Oev", 0.016, 3.844 },   //
    { "  C" ,    "B",    "D",    "CYTIDINE-5'-MONOPHOSPHATE" ,                    "O4' C1' C2' C3' C4'", 0.272, "Oev", 0.017, 4.170 },   //
    { " DA" ,	 "B",	 "D", 	 "2'-DEOXYADENOSINE-5'-MONOPHOSPHATE" ,           "O4' C1' C2' C3' C4'", 0.193, "ev4", 0.015, 3.993 },   //
    { " DG" ,    "B",    "D",    "2'-DEOXYGUANOSINE-5'-MONOPHOSPHATE" ,           "O4' C1' C2' C3' C4'", 0.181, "ev4", 0.019, 4.432 },   //
    { " DC" ,    "B",    "D",    "2'-DEOXYCYTIDINE-5'-MONOPHOSPHATE" ,            "O4' C1' C2' C3' C4'", 0.189, "Ot4", 0.018, 3.623 },   //
    { " DT" ,    "B",    "D",    "2'-DEOXYTHYMIDINE-5'-MONOPHOSPHATE" ,           "O4' C1' C2' C3' C4'", 0.245, "ev3", 0.017, 4.176 } } ;//

    const int sugar_database_size = sizeof( sugar_database ) / sizeof( sugar_database[0] );

    bool found_in_database ( clipper::String name )
    {
        for (int i = 0; i < sugar_database_size ; i++)
            if (name.trim() == sugar_database[i].name_short.trim())
                return true;
        return false;
    } //!< returns true if found

    bool found_in_database ( std::string name )
    {
        clipper::String name_clipper = name.c_str();
        for (int i = 0; i < sugar_database_size ; i++)
            if ( name_clipper.trim() == sugar_database[i].name_short.trim() )
                return true;
        return false;
    } //!< returns true if found

    bool is_amino_acid ( std::string name )
    {
        std::vector<std::string> amino_acids = {
                                                "ALA",
                                                "ARG",
                                                "ASN",
                                                "ASP",
                                                "CYS",
                                                "GLN",
                                                "GLU",
                                                "GLY",
                                                "HIS",
                                                "ILE",
                                                "LEU",
                                                "LYS",
                                                "MET",
                                                "PHE",
                                                "PRO",
                                                "SER",
                                                "THR",
                                                "TRP",
                                                "TYR",
                                                "VAL",
                                                "HYP", // hydroxyproline
                                                "LYZ", // hydroxylysine
                                                "SEP" // phosphoserine
        };

        if ( std::find(amino_acids.begin(), amino_acids.end(), name) != amino_acids.end() )
            return true;
        else
            return false;
    }

    std::string carbname_of( std::string name )
    {
        clipper::String new_name;

        // codes for hexoses

        if      ( name == "GLC" ) new_name = "Glc"   ; // alpha
        else if ( name == "BGC" ) new_name = "Glc"   ; // beta
        else if ( name == "MAN" ) new_name = "Man"   ; // alpha
        else if ( name == "BMA" ) new_name = "Man"   ; // beta
        else if ( name == "GLA" ) new_name = "Gal"   ; // alpha
        else if ( name == "GAL" ) new_name = "Gal"   ; // beta
        else if ( name == "FUC" ) new_name = "Fuc"   ; // alpha - l - fucose
        else if ( name == "FCB" ) new_name = "Fuc"   ; // beta - d - fucose
        else if ( name == "FUL" ) new_name = "Fuc"   ; // beta - l - fucose
        else if ( name == "XYS" ) new_name = "Xyl"   ; // alpha
        else if ( name == "XYP" ) new_name = "Xyl"   ; // beta

        // codes for hexosamines
        // couldn't find codes for: ManN (either), GalN (either)

        else if ( name == "GCS" ) new_name = "GlcN"  ; // beta
        else if ( name == "PA1" ) new_name = "GlcN"  ; // alpha

        // codes for N-acetyl hexosamines
        // couldn't find codes for: ManNAc (beta)

        else if ( name == "NAG" ) new_name = "GlcNAc"; // beta
        else if ( name == "NDG" ) new_name = "GlcNAc"; // alpha
        else if ( name == "NGA" ) new_name = "GalNAc"; // beta
        else if ( name == "A2G" ) new_name = "GalNAc"; // alpha
        else if ( name == "BM3" ) new_name = "ManNAc"; // alpha
        else if ( name == "BM7" ) new_name = "ManNAc"; // beta

        // codes for acidic sugars
        // couldn't find codes for: Neu5Gc (either)
        // Neu5Gc is NGC(alpha) and NGE(beta)

        else if ( name == "SIA" ) new_name = "Neu5Ac" ; // alpha
        else if ( name == "SLB" ) new_name = "Neu5Ac" ; // beta
        else if ( name == "IDR" ) new_name = "IdoA"   ; // alpha
        else if ( name == "KDM" ) new_name = "KDN"    ; // alpha
        else if ( name == "KDN" ) new_name = "KDN"    ; // beta
        else if ( name == "BDP" ) new_name = "GlcA"   ; // beta
        else if ( name == "GCU" ) new_name = "GlcA"   ; // alpha
        else if ( name == "MAV" ) new_name = "ManA"   ; // alpha
        else if ( name == "BEM" ) new_name = "ManA"   ; // beta
        else if ( name == "GTR" ) new_name = "GalA"   ; // beta
        else if ( name == "ADA" ) new_name = "GalA"   ; // alpha

        else if ( name == "DAN" ) new_name = "NeuAc" ; // Undetermined.
        else new_name = name;

        return new_name;
    }


    bool residue_has_alternate_monomer( std::string name, bool glucose_only )
    {
        std::unordered_set<std::string> residues_with_alternate_monomer;
        if (glucose_only) 
        {
            residues_with_alternate_monomer = 
            {
            "GLC", "BGC", // Glc
            };
        }
        else
        {
            residues_with_alternate_monomer = 
            {
            "GLC", "BGC", // Glc
            "MAN", "BMA", // Man
            "GLA", "GAL", // Gal
            "NDG", "NAG", // GlcNAc
            "BM3", "BM7", // ManNAc
            "A2G", "NGA", // GalNAc
            "GCU", "BDP", // GlcA
            "MAV", "BEM", // ManA
            "ADA", "GTR" // GalA
            };
        }
        
        if (residues_with_alternate_monomer.find(name) != residues_with_alternate_monomer.end())
		    return true;
        else
            return false; 
    }

    bool residue_has_alternate_anomer( std::string name )
    {
        std::unordered_set<std::string> residues_with_alternate_anomer = {
            "GLC", "BGC",
            "MAN", "BMA",
            "GLA", "GAL",
            "FUC", "FUL",
            "FCB", "FCA",
            "XYS", "XYP",
            "GCS", "PA1",
            "NDG", "NAG",
            "A2G", "NGA",
            "BM3", "BM7",
            "SIA", "SLB",
            "KDM", "KDN",
            "GCU", "BDP",
            "MAV", "BEM",
            "ADA", "GTR"
        };

        if (residues_with_alternate_anomer.find(name) != residues_with_alternate_anomer.end())
		    return true;
        else
            return false; 
    }


    std::string alternative_anomer( std::string name )
    {
        clipper::String new_name;

        if      ( name == "GLC" ) new_name = "BGC"   ; 
        else if ( name == "BGC" ) new_name = "GLC"   ; 
        else if ( name == "MAN" ) new_name = "BMA"   ; 
        else if ( name == "BMA" ) new_name = "MAN"   ; 
        else if ( name == "GLA" ) new_name = "GAL"   ; 
        else if ( name == "GAL" ) new_name = "GLA"   ; 
        else if ( name == "FUC" ) new_name = "FUL"   ;    
        else if ( name == "FUL" ) new_name = "FUC"   ;
        else if ( name == "FCB" ) new_name = "FCA"   ;
        else if ( name == "FCA" ) new_name = "FCB"   ;  
        else if ( name == "XYS" ) new_name = "XYP"   ; 
        else if ( name == "XYP" ) new_name = "XYS"   ; 

        // codes for hexosamines
        // couldn't find codes for: ManN (either), GalN (either)

        else if ( name == "GCS" ) new_name = "PA1"  ; 
        else if ( name == "PA1" ) new_name = "GCS"  ; 

        // codes for N-acetyl hexosamines
        // couldn't find codes for: ManNAc (beta)

        else if ( name == "NAG" ) new_name = "NDG"; 
        else if ( name == "NDG" ) new_name = "NAG"; 
        else if ( name == "NGA" ) new_name = "A2G"; 
        else if ( name == "A2G" ) new_name = "NGA"; 
        else if ( name == "BM3" ) new_name = "BM7";
        else if ( name == "BM7" ) new_name = "BM3"; 

        // codes for acidic sugars
        // couldn't find codes for: Neu5Gc (either)
        // Neu5Gc is NGC(alpha) and NGE(beta)

        else if ( name == "SIA" ) new_name = "SLB" ; 
        else if ( name == "SLB" ) new_name = "SIA" ; 
        else if ( name == "KDM" ) new_name = "KDN"    ; 
        else if ( name == "KDN" ) new_name = "KDM"    ; 
        else if ( name == "BDP" ) new_name = "GCU"   ; 
        else if ( name == "GCU" ) new_name = "BDP"   ; 
        else if ( name == "MAV" ) new_name = "BEM"   ; 
        else if ( name == "BEM" ) new_name = "MVA"   ; 
        else if ( name == "GTR" ) new_name = "ADA"   ; 
        else if ( name == "ADA" ) new_name = "GTR"   ;

        else new_name = "Error";

        return new_name;
    }

    
    std::vector<std::string> alternative_monomer( std::string name )
    {
        std::set<std::string> monomer_bin;

        std::vector<std::set<std::string>> bins_alternate_monomers = {
            { std::set<std::string>( { "GLC", "MAN", "GLA" } ) }, // alpha ucoses.
            { std::set<std::string>( { "BGC", "BMA", "GAL" } ) }, // beta ucoses. 
            { std::set<std::string>( { "NDG", "BM3", "A2G" } ) }, // alpha NAcs
            { std::set<std::string>( { "NAG", "BM7", "NGA" } ) }, // beta NAcs
            { std::set<std::string>( { "GCU", "MAV", "ADA" } ) }, // alpha As
            { std::set<std::string>( { "BDP", "BEM", "GTR" } ) }, // beta As
        };

        for(int i = 0; i < bins_alternate_monomers.size(); i++)
        {  
            if (bins_alternate_monomers[i].find(name) != bins_alternate_monomers[i].end())
                monomer_bin = bins_alternate_monomers[i];
        }

        std::vector<std::string> alternative_monomers;
        for (auto element : monomer_bin)
        {
            if(element != name) alternative_monomers.push_back(element);
        }

        return alternative_monomers;
    }



    // Long term task: expand this list after running on entire pdb database. 
    std::string convert_to_wurcs_residue_code( std::string name )
    {
        clipper::String wurcs_residue_code;

        // codes for hexoses

        if      ( name == "GLC" ) wurcs_residue_code = "a2122h-1a_1-5"   ; // alpha
        else if ( name == "BGC" ) wurcs_residue_code = "a2122h-1b_1-5"   ; // beta
        else if ( name == "MAN" ) wurcs_residue_code = "a1122h-1a_1-5"   ; // alpha
        else if ( name == "BMA" ) wurcs_residue_code = "a1122h-1b_1-5"   ; // beta
        else if ( name == "GLA" ) wurcs_residue_code = "a2112h-1a_1-5"   ; // alpha
        else if ( name == "GAL" ) wurcs_residue_code = "a2112h-1b_1-5"   ; // beta
        else if ( name == "FCA" ) wurcs_residue_code = "a2112m-1a_1-5"   ; // alpha - d - fucose
        else if ( name == "FCB" ) wurcs_residue_code = "a2112m-1b_1-5"   ; // beta - d - fucose
        else if ( name == "FUC" ) wurcs_residue_code = "a1221m-1a_1-5"   ; // alpha - l - fucose
        else if ( name == "FUL" ) wurcs_residue_code = "a1221m-1b_1-5"   ; // beta - l - fucose
        else if ( name == "XYS" ) wurcs_residue_code = "a212h-1a_1-5"   ; // alpha
        else if ( name == "XYP" ) wurcs_residue_code = "a212h-1b_1-5"   ; // beta
        else if ( name == "Z9D" ) wurcs_residue_code = "a21d2h-1b_1-5"   ; // beta

        // codes for hexosamines
        // couldn't find codes for: ManN (either), GalN (either)

        else if ( name == "GCS" ) wurcs_residue_code = "a2122h-1b_1-5_2*N"  ; // beta
        else if ( name == "PA1" ) wurcs_residue_code = "a2122h-1a_1-5_2*N"  ; // alpha

        // codes for N-acetyl hexosamines
        // couldn't find codes for: ManNAc (beta)

        else if ( name == "NAG" ) wurcs_residue_code = "a2122h-1b_1-5_2*NCC/3=O"; // beta
        else if ( name == "NDG" ) wurcs_residue_code = "a2122h-1a_1-5_2*NCC/3=O"; // alpha
        else if ( name == "NGA" ) wurcs_residue_code = "a2112h-1b_1-5_2*NCC/3=O"; // beta
        else if ( name == "A2G" ) wurcs_residue_code = "a2112h-1a_1-5_2*NCC/3=O"; // alpha
        else if ( name == "BM3" ) wurcs_residue_code = "a1122h-1a_1-5_2*NCC/3=O"; // alpha
        else if ( name == "BM7" ) wurcs_residue_code = "a1122h-1b_1-5_2*NCC/3=O"; // beta

        // codes for acidic sugars
        // couldn't find codes for: Neu5Gc (either)
        // Neu5Gc is NGC(alpha) and NGE(beta)

        else if ( name == "SIA" ) wurcs_residue_code = "Aad21122h-2a_2-6_5*NCC/3=O" ; // alpha
        else if ( name == "SLB" ) wurcs_residue_code = "Aad21122h-2b_2-6_5*NCC/3=O" ; // beta
        else if ( name == "IDR" ) wurcs_residue_code = "a2121A-1a_1-5"   ; // alpha
        else if ( name == "KDM" ) wurcs_residue_code = "Aad21122h-2a_2-6"    ; // alpha
        else if ( name == "KDN" ) wurcs_residue_code = "Aad21122h-2b_2-6"    ; // beta
        else if ( name == "BDP" ) wurcs_residue_code = "a2122A-1b_1-5"   ; // beta
        else if ( name == "GCU" ) wurcs_residue_code = "a2122A-1a_1-5"   ; // alpha
        else if ( name == "MAV" ) wurcs_residue_code = "a1122A-1a_1-5"   ; // alpha
        else if ( name == "BEM" ) wurcs_residue_code = "a1122A-1b_1-5"   ; // beta
        else if ( name == "GTR" ) wurcs_residue_code = "a2112A-1b_1-5"   ; // beta
        else if ( name == "ADA" ) wurcs_residue_code = "a2112A-1a_1-5"   ; // alpha
        else if ( name == "LGU" ) wurcs_residue_code = "a1121A-1a_1-5" ; // alpha-L-Gulopyranuronic acid
        else if ( name == "GUP" ) wurcs_residue_code = "a1121h-1a_1-5" ; // alpha-l-Gulopyranoside
        else if ( name == "DAN" ) wurcs_residue_code = "Aad21122h-2x_2-6_5*N" ; // Undetermined. 

        // More unique residues
        else if ( name == "M6D" ) wurcs_residue_code = "a1122h-1b_1-5_6*OPO/3O/3=O" ; // beta-D-Mannose 6-phosphate 
        else if ( name == "NAA" ) wurcs_residue_code = "a2222h-1b_1-5_2*NCC/3=O" ; // 2-acetamido-2-deoxy-beta-D-Allopyranose
        else if ( name == "NGK" ) wurcs_residue_code = "a2112h-1a_1-5_2*NCC/3=O_4*OSO/3=O/3=O" ; // 2-acetamido-4-O-sulfono-2-deoxy-alpha-D-Galactopyranose
       
        else if ( name == "FRU" ) wurcs_residue_code = "ha122h-2b_2-5";

        // Rare bois, added 01/12/2021 onwards.
        else if ( name == "AC1" ) wurcs_residue_code = "a2122m-1a_1-5_4*NC^SC^SC^SC^RCCO/7=^ZC$3/6O/5O/4O";
        else if ( name == "RY7" ) wurcs_residue_code = "a2122m-1a_1-5_4*NC^SC^SC^SC^RC^RCO/7C$3/6O/5O/4O";
        else if ( name == "MAG" ) wurcs_residue_code = "a2122h-1b_1-5_1*OC_2*NCC/3=O";
        else if ( name == "SGN" ) wurcs_residue_code = "a2122h-1a_1-5_2*NSO/3=O/3=O_6*OSO/3=O/3=O";
        else if ( name == "UAP" ) wurcs_residue_code = "a21eEA-1a_1-5_2*OSO/3=O/3=O";
        else if ( name == "MGL" ) wurcs_residue_code = "a2122h-1b_1-5_1*OC";
        else if ( name == "SGC" ) wurcs_residue_code = "a2122h-1b_1-5_4*S";
        else if ( name == "NGS" ) wurcs_residue_code = "a2122h-1b_1-5_2*NCC/3=O_6*OSO/3=O/3=O";
        else if ( name == "MMA" ) wurcs_residue_code = "a1122h-1a_1-5_1*OC";
        else if ( name == "KDA" ) wurcs_residue_code = "Aad1122h-2a_2-6_2*OCC=C";
        else if ( name == "KDO" ) wurcs_residue_code = "Aad1122h-2a_2-6";
        else if ( name == "KDB" ) wurcs_residue_code = "Aazzd22h-2a_2-6";
        else if ( name == "DGS" ) wurcs_residue_code = "a2112h-1a_1-5_3-6_2*OSO/3=O/3=O";
        else if ( name == "G4S" ) wurcs_residue_code = "a2112h-1b_1-5_4*OSO/3=O/3=O";
        else if ( name == "GMH" ) wurcs_residue_code = "a11221h-1a_1-5";
        else if ( name == "TVS" ) wurcs_residue_code = "a2122h-1b_1-5_1*OCC=C_2*NCC/3=O";
        else if ( name == "TVV" ) wurcs_residue_code = "a2112h-1b_1-5_3*OCC#C";
        else if ( name == "RAM" ) wurcs_residue_code = "a2211m-1a_1-5";
        else if ( name == "GAD" ) wurcs_residue_code = "a21eEA-1a_1-5";
        else if ( name == "IDS" ) wurcs_residue_code = "a2121A-1a_1-5_2*OSO/3=O/3=O";
        else if ( name == "KO1" ) wurcs_residue_code = "Aa11122h-2a_2-6";
        else if ( name == "G2F" ) wurcs_residue_code = "a2122h-1a_1-5_2*F";
        else if ( name == "NFG" ) wurcs_residue_code = "a2122h-1b_1-5_1*O(C^ZC^EC^EC^EC^ZC^E$3)/6NO/9=O/4NO/12=O_2*F";
        else if ( name == "MAW" ) wurcs_residue_code = "a11eEA-1a_1-5";
        else if ( name == "ASG" ) wurcs_residue_code = "a2112h-1b_1-5_2*NCC/3=O_4*OSO/3=O/3=O";
        else if ( name == "GU4" ) wurcs_residue_code = "a2122h-1a_1-5_2*OSO/3=O/3=O_3*OSO/3=O/3=O_4*OSO/3=O/3=O_6*OSO/3=O/3=O";
        else if ( name == "YYJ" ) wurcs_residue_code = "ha122h-2b_2-5_1*OSO/3=O/3=O_3*OSO/3=O/3=O_4*OSO/3=O/3=O_6*OSO/3=O/3=O";
        else if ( name == "Z9N" ) wurcs_residue_code = "ha122h-2a_2-5";
        else if ( name == "G6P" ) wurcs_residue_code = "a2122h-1a_1-5_6*OPO/3O/3=O";
        else if ( name == "5N6" ) wurcs_residue_code = "Aad21122h-2a_2-6_5*NCC/3=O_9*OCC/3=O";
        else if ( name == "TYV" ) wurcs_residue_code = "a1d22m-1a_1-5";
        else if ( name == "PKM" ) wurcs_residue_code = "Aad21122h-2a_2-6_4*OCC/3=O_5*NCC/3=O";
        else if ( name == "ABE" ) wurcs_residue_code = "a2d12m-1a_1-5";
        else if ( name == "AAL" ) wurcs_residue_code = "a1221h-1a_1-5_3-6";
        else if ( name == "Z3Q" ) wurcs_residue_code = "a2122h-1b_1-5_1*OCCN=^ZN=N_2*NCC/3=O";
        else if ( name == "BHG" ) wurcs_residue_code = "a2112h-1b_1-5_1*OCCCCCC";
        else if ( name == "GNS" ) wurcs_residue_code = "a2122h-1a_1-5_2*NSO/3=O/3=O";
        else if ( name == "AMV" ) wurcs_residue_code = "a2122h-1b_1-5_1*OC_2*NCC/3=O_3*OC^RCO/4=O/3C";
        else if ( name == "AMU" ) wurcs_residue_code = "a2122h-1b_1-5_2*NCC/3=O_3*OC^RCO/4=O/3C";
        else if ( name == "1GL" ) wurcs_residue_code = "ad112m-1a_1-5_4*OC";
        else if ( name == "ARI" ) wurcs_residue_code = "ad212m-1b_1-5_4*OCC/3=O";
        else if ( name == "CDR" ) wurcs_residue_code = "ad222m-1b_1-5";
        else if ( name == "ERI" ) wurcs_residue_code = "ad611m-1a_1-5_3*C_4*OCC/3=O";
        else if ( name == "9WJ" ) wurcs_residue_code = "a1122m-1b_1-5_4*NCC/3=O";
        else if ( name == "IDY" ) wurcs_residue_code = "a2121A-1a_1-5_1*OC_2*OSO/3=O/3=O";
        else if ( name == "GNX" ) wurcs_residue_code = "a2122h-1a_1-5_2*NSO/3=O/3=O_3*OSO/3=O/3=O";
        else if ( name == "SUS" ) wurcs_residue_code = "a2122h-1a_1-5_2*NSO/3=O/3=O_3*OSO/3=O/3=O_6*OSO/3=O/3=O";
        else if ( name == "Z6W" ) wurcs_residue_code = "ha122d2ddddddm-2b_2-5";
        else if ( name == "83Y" ) wurcs_residue_code = "a2211m-1a_1-5_3*OSO/3=O/3=O";
        else if ( name == "GCD" ) wurcs_residue_code = "a21eEA-1a_1-5";
        else if ( name == "RR7" ) wurcs_residue_code = "ad122h-1b_1-5";
        else if ( name == "DDA" ) wurcs_residue_code = "ad122m-1b_1-5";
        else if ( name == "DDL" ) wurcs_residue_code = "ad112m-1b_1-5";
        else if ( name == "MDA" ) wurcs_residue_code = "ad622m-1b_1-5_3*C";
        else if ( name == "AQA" ) wurcs_residue_code = "a21eEA-1b_1-5";
        else if ( name == "9RN" ) wurcs_residue_code = "a2112h-1a_1-5_3-6";
        else if ( name == "YIO" ) wurcs_residue_code = "a2112h-1b_1-5_1*S";
        else if ( name == "GCO" ) wurcs_residue_code = "A2122h"; // linear form of sugar, not implemented in the above table. 
        else if ( name == "NBG" ) wurcs_residue_code = "a2122h-1b_1-5_1*NCC/3=O";
        else if ( name == "MBG" ) wurcs_residue_code = "a2112h-1b_1-5_1*OC";
        else if ( name == "AHR" ) wurcs_residue_code = "a211h-1a_1-4";
        else if ( name == "FUB" ) wurcs_residue_code = "a211h-1b_1-4";
        else if ( name == "GYP" ) wurcs_residue_code = "a2122h-1a_1-5_1*OC";
        else if ( name == "GZL" ) wurcs_residue_code = "a2112h-1b_1-4";
        else if ( name == "AIG" ) wurcs_residue_code = "a2112h-1b_1-5_1*OCCCCCC_3*N";
        else if ( name == "AOG" ) wurcs_residue_code = "a2112h-1b_1-5_1*OCCCCCCCC_3*N";
        else if ( name == "BXY" ) wurcs_residue_code = "a122h-1a_1-4";
        else if ( name == "YZ0" ) wurcs_residue_code = "a1122h-1b_1-5_1*OC";
        else if ( name == "Z4Y" ) wurcs_residue_code = "a1122h-1a_1-5_6*S";
        else if ( name == "DFX" ) wurcs_residue_code = "h212h_1-5_2*F";
        else if ( name == "Z4R" ) wurcs_residue_code = "a1122h-1a_1-5_1*OC_3*S";
        else if ( name == "GLD" ) wurcs_residue_code = "a21d2m-1a_1-5";
        else if ( name == "G6D" ) wurcs_residue_code = "a2122m-1a_1-5";
        else if ( name == "V3P" ) wurcs_residue_code = "a2122h-1b_1-5_1*S(C^EC^ZC^EC^EC^ZC^E$3)/6I";
        else if ( name == "TUP" ) wurcs_residue_code = "a2122h-1a_1-5_3*F";
        else if ( name == "M8C" ) wurcs_residue_code = "a2112A-1a_1-5_6*OC";
        else if ( name == "SHB" ) wurcs_residue_code = "a2112A-1b_1-5_6*OC";
        else if ( name == "MUB" ) wurcs_residue_code = "a2122h-1a_1-5_2*NCC/3=O_3*OC^RCO/4=O/3C";
        else if ( name == "GC4" ) wurcs_residue_code = "a21d2A-1b_1-5";
        else if ( name == "B6D" ) wurcs_residue_code = "a2122m-1b_1-5_2*NCC/3=O_4*NCC/3=O";
        else if ( name == "Z5L" ) wurcs_residue_code = "a1122h-1a_1-5_1*OC_2*S";
        else if ( name == "4NN" ) wurcs_residue_code = "AZz22h_1-5_2*NCC/3=O";
        else if ( name == "BG6" ) wurcs_residue_code = "a2122h-1b_1-5_6*OPO/3O/3=O";
        else if ( name == "K5B" ) wurcs_residue_code = "AOd2122h_4-7";
        else if ( name == "Z9L" ) wurcs_residue_code = "a2122h-1a_1-5_1*OC_2*OSO/3=O/3=O_3*OSO/3=O/3=O_6*OSO/3=O/3=O";
        else if ( name == "Z9K" ) wurcs_residue_code = "a2121A-1a_1-5_2*OSO/3=O/3=O_3*OC";
        else if ( name == "GU6" ) wurcs_residue_code = "a2122h-1a_1-5_2*OSO/3=O/3=O_3*OSO/3=O/3=O_6*OSO/3=O/3=O";
        else if ( name == "GU1" ) wurcs_residue_code = "a2122A-1b_1-5_2*OC_3*OC";
        else if ( name == "Z9H" ) wurcs_residue_code = "a2122h-1a_1-5_2*OSO/3=O/3=O_3*OC_4*OC_6*OSO/3=O/3=O";
        else if ( name == "NG6" ) wurcs_residue_code = "a2112h-1b_1-5_2*NCC/3=O_6*OSO/3=O/3=O";
        else if ( name == "VJ1" ) wurcs_residue_code = "a2122h-1a_1-5_1*OPO/3O/3=O_2*NCCC^RCCCCCCCCCCC/5O/3=O_3*OC^RCC^RCCCCCCCCCCC/5O/3O";
        else if ( name == "VJ4" ) wurcs_residue_code = "a2122h-1a_1-5_2*NC^SCC^ROCCCCC/7=O/5CCCCCCCC/3O_4*OPO/3O/3=O";
        else if ( name == "DGO" ) wurcs_residue_code = "zz122h_1-5";
        else if ( name == "Z61" ) wurcs_residue_code = "ad122h-1a_1-5";
        else if ( name == "SGA" ) wurcs_residue_code = "a2112h-1b_1-5_3*OSO/3=O/3=O";
        else if ( name == "U2A" ) wurcs_residue_code = "a2122h-1b_1-5_1*OC_2*S";
        else if ( name == "U1Y" ) wurcs_residue_code = "a2122h-1b_1-5_1*OC_6*S";
        else if ( name == "GCV" ) wurcs_residue_code = "a2122A-1a_1-5_4*OC";
        else if ( name == "MA3" ) wurcs_residue_code = "a2122h-1a_1-5_1*OC_4*S";
        else if ( name == "TWA" ) wurcs_residue_code = "a1222h-1b_1-5_2*OSO/3=O/3=O_3*OSO/3=O/3=O_4*OSO/3=O/3=O";
        else if ( name == "TWD" ) wurcs_residue_code = "a1211h-1a_1-5_2*OSO/3=O/3=O_3*OSO/3=O/3=O";
        else if ( name == "GTM" ) wurcs_residue_code = "a2122h-1b_1-5_1*OC_4*S";
        else if ( name == "GDA" ) wurcs_residue_code = "a2122h-1b_1-5_4*N";
        else if ( name == "SSG" ) wurcs_residue_code = "a2122h-1b_1-5_1*S_4*S";
        else if ( name == "RTV" ) wurcs_residue_code = "h1122h_1-5_2*NCC/3=O";
        else if ( name == "YYQ" ) wurcs_residue_code = "a1221h-1a_1-5_2*NCC/3=O";
        else if ( name == "SHG" ) wurcs_residue_code = "a2122h-1b_1-5_2*F";
        else if ( name == "NGC" ) wurcs_residue_code = "Aad21122h-2a_2-6_5*NCCO/3=O";
        else if ( name == "MN0" ) wurcs_residue_code = "Aad21122h-2a_2-6_2*OC_5*NCCO/3=O";
        else if ( name == "MGC" ) wurcs_residue_code = "a2112h-1a_1-5_1*OC_2*NCC/3=O";
        else if ( name == "BNG" ) wurcs_residue_code = "a2122h-1b_1-5_1*OCCCCCCCCC";
        else if ( name == "GPM" ) wurcs_residue_code = "h12122h_2-6_1*PO/2O/2=O";
        else if ( name == "GDL" ) wurcs_residue_code = "A2122h_1-5_2*NCC/3=O";
        else if ( name == "ARB" ) wurcs_residue_code = "a211h-1b_1-5";
        else if ( name == "M6P" ) wurcs_residue_code = "a1122h-1a_1-5_6*OPO/3O/3=O";
        else if ( name == "DT6" ) wurcs_residue_code = "a2122h-1b_1-5_2*NCC/3=O_4*NCC/3=O";
        else if ( name == "EGA" ) wurcs_residue_code = "a2112h-1b_1-5_1*OCC";
        else if ( name == "TQY" ) wurcs_residue_code = "a2122h-1a_1-5_6*OCCCCCCCC/3=O";
        else if ( name == "ZDO" ) wurcs_residue_code = "a2122h-1a_1-5_1*OC_2*NSO/3=O/3=O_6*OSO/3=O/3=O";
        else if ( name == "X2F" ) wurcs_residue_code = "a212h-1a_1-5_2*F";
        else if ( name == "AGL" ) wurcs_residue_code = "a2122m-1a_1-5_4*N"; // Current PDB version also doesn't have complete MOD for this UniqueRES
        else if ( name == "BOG" ) wurcs_residue_code = "a2122h-1b_1-5_1*OCCCCCCCC";
        else if ( name == "HSQ" ) wurcs_residue_code = "a2121h-1a_1-5_2*NCC/3=O";
        else if ( name == "AMG" ) wurcs_residue_code = "a2112h-1a_1-5_1*OC";
        else if ( name == "ZD0" ) wurcs_residue_code = "a1122m-1a_1-5_1*OC"; //incorrect, lacks long modification, but only 1 entry on PDB...
        else if ( name == "ZCZ" ) wurcs_residue_code = "a1122m-1a_1-5_2*OC"; //incorrect, lacks long modification, but only 1 entry on PDB...
        else if ( name == "RM4" ) wurcs_residue_code = "a2211m-1b_1-5";
        else if ( name == "7CV" ) wurcs_residue_code = "a2211m-1a_1-5_2*OC_3*OC";
        else if ( name == "XXR" ) wurcs_residue_code = "a1122m-1a_1-5";
        else if ( name == "KD5" ) wurcs_residue_code = "AOd1122h_4-7";
        else if ( name == "F6P" ) wurcs_residue_code = "ha122h-2b_2-5_6*OPO/3O/3=O";
        else if ( name == "BDR" ) wurcs_residue_code = "a222h-1b_1-4";
        else if ( name == "IDG" ) wurcs_residue_code = "a2121h-1b_1-5_2*N_6*N";
        else if ( name == "PA1" ) wurcs_residue_code = "a2122h-1a_1-5_2*N";
        else if ( name == "2GL" ) wurcs_residue_code = "ad112m-1b_1-5_4*OCC/3=O";
        else if ( name == "LDY" ) wurcs_residue_code = "a112h-1a_1-5";
        else if ( name == "GC1" ) wurcs_residue_code = "A1121h_2-6";
        else if ( name == "Z9M" ) wurcs_residue_code = "a2122h-1b_1-5_2*N_4*OPO/3O/3=O";
        else if ( name == "GP1" ) wurcs_residue_code = "a2122h-1a_1-5_1*OPO/3O/3=O_2*N";
        else if ( name == "RER" ) wurcs_residue_code = "ad621m-1a_1-5_3*C_3*N";
        else if ( name == "TMX" ) wurcs_residue_code = "a2122h-1b_1-5_2*NC/2C/2C";
        else if ( name == "BBV" ) wurcs_residue_code = "a2122h-1a_1-5_1*OC(C^EC^ZC^ZC^ZC^ZC^E$4)_2*NCC/3=O";
        else if ( name == "ASO" ) wurcs_residue_code = "h2122h_1-5";
        else if ( name == "G4D" ) wurcs_residue_code = "a2112h-1a_1-5";
        else if ( name == "BDF" ) wurcs_residue_code = "ha122h-2b_2-6";
        else if ( name == "BDG" ) wurcs_residue_code = "a2122h-1a_1-5_2*N_6*N";
        else if ( name == "RIB" ) wurcs_residue_code = "a222h-1a_1-4";
        else if ( name == "G6S" ) wurcs_residue_code = "a2112h-1b_1-5_6*OSO/3=O/3=O";
        else if ( name == "MXY" ) wurcs_residue_code = "a1221m-1b_1-5_2*OC";
        else if ( name == "WIA" ) wurcs_residue_code = "a2112h-1b_1-5_1*OC_6*S";
        else if ( name == "SHD" ) wurcs_residue_code = "a1222h-1a_1-5";
        else if ( name == "1GN" ) wurcs_residue_code = "a2112h-1b_1-5_2*N";
        else if ( name == "6PZ" ) wurcs_residue_code = "Aad22111m-2a_2-6_5*NCC/3=O_7*NCC/3=O";
        else if ( name == "6LW" ) wurcs_residue_code = "A211h_1-4_1*=NO";
        else if ( name == "ARA" ) wurcs_residue_code = "a211h-1a_1-5";
        else if ( name == "TT7" ) wurcs_residue_code = "ha122h-2b_2-5_4*OPO/3O/3=O";
        else if ( name == "4GL" ) wurcs_residue_code = "a2212h-1a_1-5";
        else if ( name == "MBF" ) wurcs_residue_code = "a1122h-1b_1-5_2*F";
        else if ( name == "DLG" ) wurcs_residue_code = "a2d12h-1b_1-5_1*OCCCCCC";
        else if ( name == "GS1" ) wurcs_residue_code = "a2122h-1b_1-5_1*S";
        else if ( name == "OPM" ) wurcs_residue_code = "a1122h-1a_1-5_1*OCCCCC";
        else if ( name == "DRI" ) wurcs_residue_code = "ad122m-1b_1-5_4*OC";
        else if ( name == "LXB" ) wurcs_residue_code = "a2212h-1b_1-5_2*NCC/3=O";
        else if ( name == "LXZ" ) wurcs_residue_code = "a1212h-1a_1-5_2*NCC/3=O";
        else if ( name == "NGZ" ) wurcs_residue_code = "a1211h-1a_1-5_2*NCC/3=O";
        else if ( name == "GL0" ) wurcs_residue_code = "a2212h-1b_1-5";
        else if ( name == "GXL" ) wurcs_residue_code = "a1221h-1a_1-5";
        else if ( name == "GM0" ) wurcs_residue_code = "a11221h-1a_1-5_4*OPO/3O/3=O";
        else if ( name == "GCN" ) wurcs_residue_code = "a2d22h-1a_1-5_2*N";
        else if ( name == "U2D" ) wurcs_residue_code = "a2122h-1a_1-5_6*OCCCCCCCCCC/3=O";
        else if ( name == "8EX" ) wurcs_residue_code = "a2112h-1b_1-5_2*NCC/3=O_4*OSO/3=O/3=O_6*OSO/3=O/3=O";
        else if ( name == "TXB" ) wurcs_residue_code = "a212h-1a_1-5_4*S";
        else if ( name == "TVG" ) wurcs_residue_code = "a2112h-1b_1-4_1*OCCC";
        else if ( name == "TRV" ) wurcs_residue_code = "ha122h-2b_2-5_6*OCCCCCCCC/3=O";
        else if ( name == "JHM" ) wurcs_residue_code = "ad122h-1a_1-5_6*OSO/3=O/3=O";

        else wurcs_residue_code = "ERROR: UNABLE TO FIND \'" + name + "\' RESIDUE CODE IN INTERNAL DATABASE";

        return wurcs_residue_code;
    } 
    
    // Consider getting rid of clipper::data::get_anomer when there is already clipper::MSugar::anomer()
    // Or maybe not, SNFG generation functions depend on this.
    std::string get_anomer( std::string name )
    {
        clipper::String anomer;

        // codes for hexoses

        if      ( name == "GLC" ) anomer = "alpha"   ; // alpha
        else if ( name == "BGC" ) anomer = "beta"   ; // beta
        else if ( name == "MAN" ) anomer = "alpha"   ; // alpha
        else if ( name == "BMA" ) anomer = "beta"   ; // beta
        else if ( name == "GLA" ) anomer = "alpha"   ; // alpha
        else if ( name == "GAL" ) anomer = "beta"   ; // beta
        else if ( name == "FCA" ) anomer = "alpha"   ; // alpha - d - fucose
        else if ( name == "FCB" ) anomer = "beta"   ; // beta - d - fucose
        else if ( name == "FUC" ) anomer = "alpha"   ; // alpha - l - fucose
        else if ( name == "FUL" ) anomer = "beta"   ; // beta - l - fucose
        else if ( name == "XYS" ) anomer = "alpha"   ; // alpha
        else if ( name == "XYP" ) anomer = "beta"   ; // beta

        // codes for hexosamines
        // couldn't find codes for: ManN (either), GalN (either)

        else if ( name == "GCS" ) anomer = "beta"  ; // beta
        else if ( name == "PA1" ) anomer = "alpha"  ; // alpha

        // codes for N-acetyl hexosamines
        // couldn't find codes for: ManNAc (beta)

        else if ( name == "NAG" ) anomer = "beta"; // beta
        else if ( name == "NDG" ) anomer = "alpha"; // alpha
        else if ( name == "NGA" ) anomer = "beta"; // beta
        else if ( name == "A2G" ) anomer = "alpha"; // alpha
        else if ( name == "BM3" ) anomer = "alpha"; // alpha
        else if ( name == "BM7" ) anomer = "beta"; // beta

        // codes for acidic sugars
        // couldn't find codes for: Neu5Gc (either)

        else if ( name == "SIA" ) anomer = "alpha" ; // alpha
        else if ( name == "SLB" ) anomer = "beta" ; // beta
        else if ( name == "IDR" ) anomer = "beta"   ; // beta
        else if ( name == "KDM" ) anomer = "alpha"    ; // alpha
        else if ( name == "KDN" ) anomer = "beta"    ; // beta
        else if ( name == "BDP" ) anomer = "beta"   ; // beta
        else if ( name == "GCU" ) anomer = "alpha"   ; // alpha
        else if ( name == "MAV" ) anomer = "alpha"   ; // alpha
        else if ( name == "BEM" ) anomer = "beta"   ; // beta
        else if ( name == "GTR" ) anomer = "beta"   ; // beta
        else if ( name == "ADA" ) anomer = "alpha"   ; // alpha
        else if ( name == "LGU" ) anomer = "alpha" ; // alpha-L-Gulopyranuronic acid
        else if ( name == "GUP" ) anomer = "alpha" ; // alpha-l-Gulopyranoside
        else if ( name == "DAN" ) anomer = "undetermined" ; // Undetermined. 

        // More unique residues
        else if ( name == "M6D" ) anomer = "beta" ; // beta-D-Mannose 6-phosphate 
        else if ( name == "NAA" ) anomer = "beta" ; // 2-acetamido-2-deoxy-beta-D-Allopyranose
        else if ( name == "NGK" ) anomer = "alpha" ; // 2-acetamido-4-O-sulfono-2-deoxy-alpha-D-Galactopyranose
        

        else anomer = "ERROR: UNABLE TO FIND \'" + name + "\' RESIDUE CODE IN INTERNAL DATABASE";

        return anomer;
    }


} // namespace data

} // namespace clipper
