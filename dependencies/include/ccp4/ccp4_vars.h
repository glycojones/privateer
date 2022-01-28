/*
     ccp4_vars.h: Standard strings for certain quantites
     Copyright (C) 2002  CCLRC, Martyn Winn

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
*/
/*
*/

/* Author: Martyn Winn */

/* Standard strings for certain quantites - for future use */

#ifndef __CCP4_VARS__
#define __CCP4_VARS__

#define MTZFILENAME "data::mtzfile::filename"
#define MTZTITLE "data::mtzfile::title"
#define MTZSPACEGROUP "data::mtzfile::spacegroup_num"        
#define MTZNUMREFLS "data::mtzfile::num_reflections"
#define MTZMNF "data::mtzfile::missing_number_flag"
#define MTZSORTORDER "data::mtzfile::sort_order"          

#define CRYSTALXTALNAME "data::crystal::crystal_name"
#define CRYSTALPNAME "data::crystal::project_name"
#define CRYSTALCELL "data::crystal::cell"

#define DATASETDNAME "data::crystal::dataset::dataset_name"
#define DATASETWAVELENGTH "data::crystal::dataset::wavelength"

#define COLUMNLABEL "data::crystal_i::dataset_i::column_i::label"
#define COLUMNTYPE "data::crystal_i::dataset_i::column_i::type"

#endif  /*!__CCP4_VARS__ */
