//  $Id: ccp4srs_defs.h $
//  =================================================================
//
//   CCP4 SRS Library: Storage, Retrievak and Search support for
//   CCP4 ligand data.
//
//   Copyright (C) Eugene Krissinel 2010-2013.
//
//   This library is free software: you can redistribute it and/or
//   modify it under the terms of the GNU Lesser General Public
//   License version 3, modified in accordance with the provisions
//   of the license to address the requirements of UK law.
//
//   You should have received a copy of the modified GNU Lesser
//   General Public License along with this library. If not, copies
//   may be downloaded from http://www.ccp4.ac.uk/ccp4license.php
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU Lesser General Public License for more details.
//
//  =================================================================
//
//    18.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  ccp4srs_defs  <interface>
//       ~~~~~~~~~
//  **** Content :  CCP4 SRS definitions
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#ifndef CCP4SRS_DEFS_H
#define CCP4SRS_DEFS_H

// THIS FILE SHOULD BE NEVER UNCLUDED IN HEADERS. USE CCP4SRS_TYPES.H
// INSTEAD

//  mmCIF tags
#define MMCIF_COMP                   "comp_"
#define MMCIF_TAG_Name               "name"
#define MMCIF_TAG_Type               "type"
#define MMCIF_TAG_Formula            "formula"

#define MMCIF_STRUCT_CHEM_COMP       "_chem_comp"

#define MMCIF_LOOP_CHEM_COMP_ATOM    "_chem_comp_atom"
#define MMCIF_TAG_COMP_ID            "comp_id"
#define MMCIF_TAG_ATOM_ID            "atom_id"
#define MMCIF_TAG_ALT_ATOM_ID        "alt_atom_id"
#define MMCIF_TAG_TYPE_SYMBOL        "type_symbol"
#define MMCIF_TAG_TYPE_ENERGY        "type_energy"
#define MMCIF_TAG_PARTIAL_CHARGE     "partial_charge"
#define MMCIF_TAG_X                  "x"
#define MMCIF_TAG_Y                  "y"
#define MMCIF_TAG_Z                  "z"
#define MMCIF_TAG_MODEL_CARTN_X      "model_Cartn_x"
#define MMCIF_TAG_MODEL_CARTN_Y      "model_Cartn_y"
#define MMCIF_TAG_MODEL_CARTN_Z      "model_Cartn_z"
#define MMCIF_TAG_PDBX_MODEL_CARTN_X "pdbx_model_Cartn_x_ideal"
#define MMCIF_TAG_PDBX_MODEL_CARTN_Y "pdbx_model_Cartn_y_ideal"
#define MMCIF_TAG_PDBX_MODEL_CARTN_Z "pdbx_model_Cartn_z_ideal"

#define MMCIF_TAG_PDBX_STEREO_CONFIG     "pdbx_stereo_config"
#define MMCIF_TAG_PDBX_LEAVING_ATOM_FLAG "pdbx_leaving_atom_flag"
#define MMCIF_TAG_PDBX_AROMATIC_FLAG     "pdbx_aromatic_flag"
#define MMCIF_TAG_VALUE_ORDER            "value_order"

#define MMCIF_LOOP_CHEM_COMP_TREE   "_chem_comp_tree"
#define MMCIF_TAG_ATOM_BACK         "atom_back"
#define MMCIF_TAG_ATOM_FORWARD      "atom_forward"
#define MMCIF_TAG_CONNECT_TYPE      "connect_type"

#define MMCIF_LOOP_CHEM_COMP_BOND   "_chem_comp_bond"
#define MMCIF_TAG_ATOM_ID_1          "atom_id_1"
#define MMCIF_TAG_ATOM_ID_2          "atom_id_2"
#define MMCIF_TAG_TYPE               "type"
#define MMCIF_TAG_VALUE_DIST         "value_dist"
#define MMCIF_TAG_VALUE_DIST_ESD     "value_dist_esd"

#define MMCIF_LOOP_CHEM_COMP_ANGLE   "_chem_comp_angle"
#define MMCIF_TAG_ATOM_ID_3          "atom_id_3"
#define MMCIF_TAG_VALUE_ANGLE        "value_angle"
#define MMCIF_TAG_VALUE_ANGLE_ESD    "value_angle_esd"

#define MMCIF_LOOP_CHEM_COMP_CHIR    "_chem_comp_chir"
#define MMCIF_TAG_ID                 "id"
#define MMCIF_TAG_ATOM_ID_CENTRE     "atom_id_centre"
#define MMCIF_TAG_VOLUME_SIGN        "volume_sign"

#define MMCIF_LOOP_CHEM_COMP_TOR     "_chem_comp_tor"
#define MMCIF_TAG_ATOM_ID_4          "atom_id_4"
#define MMCIF_TAG_PERIOD             "period"

#define MMCIF_LOOP_CHEM_COMP_PLANE_ATOM  "_chem_comp_plane_atom"
#define MMCIF_TAG_PLANE_ID               "plane_id"
#define MMCIF_TAG_DIST_ESD               "dist_esd"

#define MMCIF_LOOP_CHEM_COMP_DESCRIPTOR  "_pdbx_chem_comp_descriptor"
#define MMCIF_TAG_PROGRAM                "program"
#define MMCIF_TAG_PROGRAM_VERSION        "program_version"
#define MMCIF_TAG_DESCRIPTOR             "descriptor"


//   Entry types
/*
#define CCP4SRS_Entry_None        0x00
#define CCP4SRS_Entry_Structure   0x01
#define CCP4SRS_Entry_Link        0x02
*/

#endif // CCP4SRS_DEFS_H
