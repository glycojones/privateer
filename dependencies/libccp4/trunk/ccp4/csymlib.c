/*
     csymlib.c: C-level library for symmetry information.
     Copyright (C) 2001  CCLRC, Martyn Winn

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

/* DO NOT put Doxygen comments here - put in csymlib.h */
/* See csymlib.h for descriptions of functions */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "ccp4_parser.h"
#include "ccp4_types.h"
#include "ccp4_utils.h"
#include "csymlib.h"
#include "cvecmat.h"
#include "ccp4_errno.h"
#include "ccp4_unitcell.h"
#include "ccp4_general.h"
/* rcsid[] = "$Id$" */

/* stuff for error reporting */
#define CSYM_ERRNO(n) (CCP4_ERR_SYM | (n))

/* error defs */
#define  CSYMERR_Ok                  0
#define  CSYMERR_ParserFail          1
#define  CSYMERR_NoSyminfoFile       2
#define  CSYMERR_NullSpacegroup      3
#define  CSYMERR_NoAsuDefined        4
#define  CSYMERR_NoLaueCodeDefined   5
#define  CSYMERR_SyminfoTokensMissing 6

/* copy of constants from fortran/csymlib_f.c */
#define MAXSYMOPS 20
#define MAXLENSYMOPSTR 80

CCP4SPG *ccp4spg_load_by_standard_num(const int numspg) 
{ 
  return ccp4spg_load_spacegroup(numspg, 0, NULL, NULL, 0, NULL);
}

CCP4SPG *ccp4spg_load_by_ccp4_num(const int ccp4numspg) 
{ 
  return ccp4spg_load_spacegroup(0, ccp4numspg, NULL, NULL, 0, NULL);
}

CCP4SPG *ccp4spg_load_by_spgname(const char *spgname) 
{ 
  return ccp4spg_load_spacegroup(0, 0, spgname, NULL, 0, NULL);
}

CCP4SPG *ccp4spg_load_by_ccp4_spgname(const char *ccp4spgname) 
{ 
  return ccp4spg_load_spacegroup(0, 0, NULL, ccp4spgname, 0, NULL);
}

CCP4SPG *ccp4_spgrp_reverse_lookup(const int nsym1, const ccp4_symop *op1)
{
  return ccp4spg_load_spacegroup(0, 0, NULL, NULL, nsym1, op1);
}

CCP4SPG *ccp4spg_load_spacegroup(const int numspg, const int ccp4numspg,
         const char *spgname, const char *ccp4spgname, 
         const int nsym1, const ccp4_symop *op1) 

{ CCP4SPG *spacegroup;
  int i,j,k,l,debug=0,nsym2,symops_provided=0,ierr,ilaue;
  float sg_chb[4][4],limits[2],rot1[4][4],rot2[4][4],det;
  FILE *filein;
  char *symopfile, *ccp4dir, filerec[80];
  ccp4_symop *op2,*op3,opinv;

  /* spacegroup variables */
  int sg_num, sg_ccp4_num, sg_nsymp, sg_num_cent;
  float cent_ops[4][4];
  char sg_symbol_old[20],sg_symbol_Hall[40],sg_symbol_xHM[20],
       sg_point_group[20],sg_patt_group[40];
  char sg_basisop[80],sg_symop[192][80],sg_cenop[4][80];
  char sg_asu_descr[80], map_asu_x[12], map_asu_y[12], map_asu_z[12];    
  char map_asu_ccp4_x[12], map_asu_ccp4_y[12], map_asu_ccp4_z[12]; 

  /* For cparser */
  CCP4PARSERARRAY *parser;
  CCP4PARSERTOKEN *token=NULL;
  char *key;
  int iprint=0;

  /* initialisations */
  sg_nsymp=0;
  sg_num_cent=0;

  if (nsym1) symops_provided=1;

  spacegroup = (CCP4SPG *) ccp4_utils_malloc(sizeof(CCP4SPG));

  if (debug) {
    printf(" Entering ccp4spg_load_spacegroup, with arguments %d %d %d \n",
        numspg,ccp4numspg,nsym1);
    if (spgname) printf(" spgname = %s \n",spgname);
    if (ccp4spgname) printf(" ccp4spgname = %s \n",ccp4spgname);
    for (i = 0; i < nsym1; ++i) {
      printf(" %f %f %f \n",op1[i].rot[0][0],op1[i].rot[0][1],op1[i].rot[0][2]);
      printf(" %f %f %f \n",op1[i].rot[1][0],op1[i].rot[1][1],op1[i].rot[1][2]);
      printf(" %f %f %f \n",op1[i].rot[2][0],op1[i].rot[2][1],op1[i].rot[2][2]);
      printf(" %f %f %f \n\n",op1[i].trn[0],op1[i].trn[1],op1[i].trn[2]);
    }
  }

  /* if we are searching with symops, make sure translations are modulo 1 */
  if (symops_provided) {
    op3 = (ccp4_symop *) ccp4_utils_malloc(nsym1*sizeof(ccp4_symop));
    for (i = 0; i < nsym1; ++i) {
     for (k = 0; k < 3; ++k) {
      for (l = 0; l < 3; ++l) {
        op3[i].rot[k][l]=op1[i].rot[k][l];
      }
      op3[i].trn[k] = op1[i].trn[k];
     }
     ccp4spg_norm_trans(&op3[i]);
    }
  }

  /* Open the symop file: */
  if (!(symopfile = getenv("SYMINFO"))) {
    if (debug)
      printf("Environment variable SYMINFO not set ... guessing location of symmetry file. \n");
    if (!(ccp4dir = getenv("CLIBD"))) {
      printf("Environment variable CLIBD not set ... big trouble! \n");
      return NULL;
    }

    symopfile = ccp4_utils_malloc((strlen(ccp4dir)+22)*sizeof(char));
    strcpy(symopfile,ccp4_utils_joinfilenames(ccp4dir,"syminfo.lib"));
    symopfile[strlen(ccp4dir)+21] = '\0';
    ccp4printf(1," SYMINFO file set to %s \n",symopfile);
  } else {
    if (debug) {
      ccp4printf(1,"\n Spacegroup information obtained from library file: \n");
      ccp4printf(1," Logical Name: SYMINFO   Filename: %s\n\n",symopfile);
    }
  }

  filein = fopen(symopfile,"r");
  if (!filein) {
    ccp4_signal(CSYM_ERRNO(CSYMERR_NoSyminfoFile),"ccp4spg_load_spacegroup",NULL); 
    return NULL;
  }

  if (!(getenv("SYMINFO"))) free(symopfile);

  parser = ccp4_parse_start(20);
  if (parser == NULL) 
    ccp4_signal(CSYM_ERRNO(CSYMERR_ParserFail),"ccp4spg_load_spacegroup",NULL);
  /* "=" is used in map asu fields, so remove it as delimiter */
  ccp4_parse_delimiters(parser," \t,",",");
  /* Set some convenient pointers to members of the parser array */
  key   = parser->keyword;
  token = parser->token;

  if (debug) 
    printf(" parser initialised \n");

  while (fgets(filerec,80,filein)) {

    /* If syminfo.lib comes from a DOS platform, and we are on
       unix, need to strip spurious \r character. Note this is
       necessary because we have removed \r as parser delimiter. */
    if (strlen(filerec) > 1){
      if (filerec[strlen(filerec)-2]=='\r') {
        filerec[strlen(filerec)-2]='\n';
        filerec[strlen(filerec)-1]='\0';
      }
    }

    if (strlen(filerec) > 1) {

      ccp4_parser(filerec, 80, parser, iprint);

      if (ccp4_keymatch(key, "number")) {
        if (parser->ntokens < 2) {
           printf("Current SYMINFO line = %s\n",filerec);
           ccp4_signal(CCP4_ERRLEVEL(2) | CSYM_ERRNO(CSYMERR_SyminfoTokensMissing),"ccp4spg_load_spacegroup",NULL);
        } else {
          sg_num = (int) token[1].value;
        }
      }

      if (ccp4_keymatch(key, "basisop")) {
        strcpy(sg_basisop,filerec+8);
      }

      if (ccp4_keymatch(key, "symbol")) {
        if (parser->ntokens < 3) {
           printf("Current SYMINFO line = %s\n",filerec);
           ccp4_signal(CCP4_ERRLEVEL(2) | CSYM_ERRNO(CSYMERR_SyminfoTokensMissing),"ccp4spg_load_spacegroup",NULL);
        } else {
          if (strcmp(token[1].fullstring,"ccp4") == 0)
            sg_ccp4_num = (int) token[2].value;
          if (strcmp(token[1].fullstring,"Hall") == 0)
            strcpy(sg_symbol_Hall,token[2].fullstring);
          if (strcmp(token[1].fullstring,"xHM") == 0)
            strcpy(sg_symbol_xHM,token[2].fullstring);
          if (strcmp(token[1].fullstring,"old") == 0)
            strcpy(sg_symbol_old,token[2].fullstring);
          if (strcmp(token[1].fullstring,"patt") == 0)
            strcpy(sg_patt_group,token[3].fullstring);
          if (strcmp(token[1].fullstring,"pgrp") == 0)
            strcpy(sg_point_group,token[3].fullstring);
        }
      }

      if (ccp4_keymatch(key, "hklasu")) {
        if (parser->ntokens < 3) {
           printf("Current SYMINFO line = %s\n",filerec);
           ccp4_signal(CCP4_ERRLEVEL(2) | CSYM_ERRNO(CSYMERR_SyminfoTokensMissing),"ccp4spg_load_spacegroup",NULL);
        } else {
          if (strcmp(token[1].fullstring,"ccp4") == 0)
            strcpy(sg_asu_descr,token[2].fullstring);
        }
      }

      if (ccp4_keymatch(key, "mapasu")) {
        if (parser->ntokens < 5) {
           printf("Current SYMINFO line = %s\n",filerec);
           ccp4_signal(CCP4_ERRLEVEL(2) | CSYM_ERRNO(CSYMERR_SyminfoTokensMissing),"ccp4spg_load_spacegroup",NULL);
        } else {
          if (strcmp(token[1].fullstring,"zero") == 0) {
            strcpy(map_asu_x,token[2].fullstring);
            strcpy(map_asu_y,token[3].fullstring);
            strcpy(map_asu_z,token[4].fullstring);
          } else if (strcmp(token[1].fullstring,"ccp4") == 0) {
            strcpy(map_asu_ccp4_x,token[2].fullstring);
            strcpy(map_asu_ccp4_y,token[3].fullstring);
            strcpy(map_asu_ccp4_z,token[4].fullstring);
          }
        }
      }

      if (ccp4_keymatch(key, "symop")) {
        strcpy(sg_symop[sg_nsymp++],filerec+6);
      }

      if (ccp4_keymatch(key, "cenop")) {
        strcpy(sg_cenop[sg_num_cent++],filerec+6);
      }

      if (ccp4_keymatch(key, "end_spacegroup")) {
        /* end of spacegroup block, so check if right one */
        if (numspg) {
          if (sg_num == numspg)
            break;
        } else if (ccp4numspg) {
          if (sg_ccp4_num == ccp4numspg)
            break;
        } else if (spgname) {
          if (ccp4spg_name_equal_to_lib(sg_symbol_xHM,spgname))
            break;
        } else if (ccp4spgname) {
          if (ccp4spg_name_equal_to_lib(sg_symbol_old,ccp4spgname))
            break;
          if (ccp4spg_name_equal_to_lib(sg_symbol_xHM,ccp4spgname))
            break;
        } else if (symops_provided) {
          nsym2 = sg_nsymp*sg_num_cent;
          if (nsym2 == nsym1) {
            op2 = (ccp4_symop *) ccp4_utils_malloc(nsym2*sizeof(ccp4_symop));
            for (i = 0; i < sg_num_cent; ++i) {
             symop_to_mat4(sg_cenop[i],sg_cenop[i]+strlen(sg_cenop[i]),cent_ops[0]);
             for (j = 0; j < sg_nsymp; ++j) {
              symop_to_mat4(sg_symop[j],sg_symop[j]+strlen(sg_symop[j]),rot2[0]);
              ccp4_4matmul(rot1,(const float (*)[4])cent_ops,(const float (*)[4])rot2);
              op2[i*sg_nsymp+j] = mat4_to_rotandtrn((const float (*)[4])rot1);
	      /* combination of primitive and centering operators can 
                 produce translations greater than one. */
              ccp4spg_norm_trans(&op2[i*sg_nsymp+j]);
	     }
            }
	    /* op3 are requested operators and op2 are from SYMINFO file */
            if (ccp4_spgrp_equal(nsym1,op3,nsym2,op2)) {
              if (debug) printf(" ops match for sg %d ! \n",sg_num);
              free(op2);
              break;
            }
	    free(op2);
	  }
        }
        sg_nsymp = 0;
        sg_num_cent = 0;
        sg_symbol_xHM[0]='\0';
        sg_symbol_old[0]='\0';
      }
    }
  }
  if (symops_provided) free(op3);

  if (debug) 
    printf(" parser finished \n");

  /* Finished with the parser array */
  ccp4_parse_end(parser);
  fclose(filein);

  if (!sg_nsymp) {
    ccp4printf(0," Failed to find spacegroup in SYMINFO! \n");
    return NULL;
  } 

  /* extract various symbols for spacegroup */
  spacegroup->spg_num = sg_num;
  spacegroup->spg_ccp4_num = sg_ccp4_num;
  strcpy(spacegroup->symbol_Hall,sg_symbol_Hall);
  strcpy(spacegroup->symbol_xHM,sg_symbol_xHM);
  strcpy(spacegroup->symbol_old,sg_symbol_old);
  strcpy(spacegroup->point_group,"PG");
  strcpy(spacegroup->point_group+2,sg_point_group);
  if (sg_num <= 2) {
    strcpy(spacegroup->crystal,"TRICLINIC");
  } else if (sg_num >= 3 && sg_num <= 15) {
    strcpy(spacegroup->crystal,"MONOCLINIC");
  } else if (sg_num >= 16 && sg_num <= 74) {
    strcpy(spacegroup->crystal,"ORTHORHOMBIC");
  } else if (sg_num >= 75 && sg_num <= 142) {
    strcpy(spacegroup->crystal,"TETRAGONAL");
  } else if (sg_num >= 143 && sg_num <= 167) {
    strcpy(spacegroup->crystal,"TRIGONAL");
  } else if (sg_num >= 168 && sg_num <= 194) {
    strcpy(spacegroup->crystal,"HEXAGONAL");
  } else if (sg_num >= 195 && sg_num <= 230) {
    strcpy(spacegroup->crystal,"CUBIC");
  } else {
    strcpy(spacegroup->crystal," ");
  }

  if (debug) 
    printf(" Read in details of spacegroup %d %d \n",sg_num,sg_ccp4_num);

  /* change of basis */
  if (debug) 
    printf(" Change of basis %s \n",sg_basisop);
  symop_to_mat4(sg_basisop,sg_basisop+strlen(sg_basisop),sg_chb[0]);
  for (i = 0; i < 3; ++i) {
   for (j = 0; j < 3; ++j) {
    spacegroup->chb[i][j] = sg_chb[i][j];
   }
  }
  if (debug)
    for (k = 0; k < 3; ++k) 
      printf("chb: %f %f %f\n",spacegroup->chb[k][0],
           spacegroup->chb[k][1],spacegroup->chb[k][2]);

  /* symmetry operators */
  spacegroup->nsymop_prim = sg_nsymp;
  spacegroup->nsymop = sg_nsymp*sg_num_cent;
  spacegroup->symop = (ccp4_symop *) ccp4_utils_malloc(spacegroup->nsymop*sizeof(ccp4_symop));
  spacegroup->invsymop = (ccp4_symop *) ccp4_utils_malloc(spacegroup->nsymop*sizeof(ccp4_symop));
  if (symops_provided) {
   for (i = 0; i < nsym1; ++i) {
    opinv = ccp4_symop_invert(op1[i]);
    for (k = 0; k < 3; ++k) {
      for (l = 0; l < 3; ++l) {
        spacegroup->symop[i].rot[k][l]=op1[i].rot[k][l];
        spacegroup->invsymop[i].rot[k][l]=opinv.rot[k][l];
      }
      spacegroup->symop[i].trn[k] = op1[i].trn[k];
      spacegroup->invsymop[i].trn[k] = opinv.trn[k];
    }
   }
  } else {
   for (i = 0; i < sg_num_cent; ++i) {
    symop_to_mat4(sg_cenop[i],sg_cenop[i]+strlen(sg_cenop[i]),cent_ops[0]);
    for (j = 0; j < sg_nsymp; ++j) {
     strncpy(filerec,sg_symop[j],80);   /* symop_to_mat4 overwrites later sg_symop */
     symop_to_mat4(filerec,filerec+strlen(filerec),rot2[0]);
     ccp4_4matmul(rot1,(const float (*)[4])cent_ops,(const float (*)[4])rot2);
     det=invert4matrix((const float (*)[4])rot1,rot2);
     if (debug) printf("symop determinant: %f\n",det);
     for (k = 0; k < 3; ++k) {
      for (l = 0; l < 3; ++l) {
        spacegroup->symop[i*sg_nsymp+j].rot[k][l]=rot1[k][l];
        spacegroup->invsymop[i*sg_nsymp+j].rot[k][l]=rot2[k][l];
      }
      spacegroup->symop[i*sg_nsymp+j].trn[k] = rot1[k][3];
      spacegroup->invsymop[i*sg_nsymp+j].trn[k] = rot2[k][3];
     }
     /* unless symops have been provided, store normalised operators */
     ccp4spg_norm_trans(&spacegroup->symop[i*sg_nsymp+j]);
     ccp4spg_norm_trans(&spacegroup->invsymop[i*sg_nsymp+j]);
    }
   }
  }
  if (debug) 
   for (i = 0; i < sg_num_cent; ++i) 
    for (j = 0; j < sg_nsymp; ++j) {
     for (k = 0; k < 3; ++k) 
      printf("rot/trn: %f %f %f %f\n",spacegroup->symop[i*sg_nsymp+j].rot[k][0],
           spacegroup->symop[i*sg_nsymp+j].rot[k][1],
           spacegroup->symop[i*sg_nsymp+j].rot[k][2],
           spacegroup->symop[i*sg_nsymp+j].trn[k]);
     for (k = 0; k < 3; ++k) 
      printf("inv rot/trn: %f %f %f %f\n",spacegroup->invsymop[i*sg_nsymp+j].rot[k][0],
           spacegroup->invsymop[i*sg_nsymp+j].rot[k][1],
           spacegroup->invsymop[i*sg_nsymp+j].rot[k][2],
           spacegroup->invsymop[i*sg_nsymp+j].trn[k]);
    }

  /* reciprocal asymmetric unit */
  strcpy(spacegroup->asu_descr,sg_asu_descr);

  /* Select ASU function (referred to default basis) from asu desc */
  /* Also infer Laue and Patterson groups. This uses additional
     information from spacegroup name. In general, we use Hall symbol
     because syminfo.lib is missing a few xHM symbols. However, we
     need to use the latter for R vs. H settings. */
  ierr = 1;
  ilaue = 1;

  if ( strcmp( sg_asu_descr, "l>0 or (l==0 and (h>0 or (h==0 and k>=0)))" ) == 0 ) {
     spacegroup->asufn = &ASU_1b;
     ilaue = ccp4spg_load_laue(spacegroup,3);
     spacegroup->npatt = 2;
     strcpy(spacegroup->patt_name,"P-1");
     ierr = 0;
  }
  if ( strcmp( sg_asu_descr, "k>=0 and (l>0 or (l=0 and h>=0))" ) == 0 ) {
     spacegroup->asufn = &ASU_2_m;
     ilaue = ccp4spg_load_laue(spacegroup,4);
     if (strchr(spacegroup->symbol_Hall,'P')) {
       spacegroup->npatt = 10;
       strcpy(spacegroup->patt_name,"P2/m");
     } else if (strchr(spacegroup->symbol_Hall,'C')) {
       spacegroup->npatt = 12;
       strcpy(spacegroup->patt_name,"C2/m");
     }
     ierr = 0;
  }
  if ( strcmp( sg_asu_descr, "h>=0 and k>=0 and l>=0" ) == 0 ) {
     spacegroup->asufn = &ASU_mmm;
     ilaue = ccp4spg_load_laue(spacegroup,6);
     if (strchr(spacegroup->symbol_Hall,'P')) {
       spacegroup->npatt = 47;
       strcpy(spacegroup->patt_name,"Pmmm");
     } else if (strchr(spacegroup->symbol_Hall,'C')) {
       spacegroup->npatt = 65;
       strcpy(spacegroup->patt_name,"Cmmm");
     } else if (strchr(spacegroup->symbol_Hall,'I')) {
       spacegroup->npatt = 71;
       strcpy(spacegroup->patt_name,"Immm");
     } else if (strchr(spacegroup->symbol_Hall,'F')) {
       spacegroup->npatt = 69;
       strcpy(spacegroup->patt_name,"Fmmm");
     }
     ierr = 0;
  }
  if ( strcmp( sg_asu_descr, "l>=0 and ((h>=0 and k>0) or (h=0 and k=0))" ) == 0  &&
       strcmp( sg_patt_group, "4/m" ) == 0 ) {
     spacegroup->asufn = &ASU_4_m;
     ilaue = ccp4spg_load_laue(spacegroup,7);
     spacegroup->nlaue = 7;
     strcpy(spacegroup->laue_name,"4/m");
     spacegroup->laue_sampling[0] = 4;
     spacegroup->laue_sampling[1] = 4;
     spacegroup->laue_sampling[2] = 8;
     if (strchr(spacegroup->symbol_Hall,'P')) {
       spacegroup->npatt = 83;
       strcpy(spacegroup->patt_name,"P4/m");
     } else if (strchr(spacegroup->symbol_Hall,'I')) {
       spacegroup->npatt = 87;
       strcpy(spacegroup->patt_name,"I4/m");
     }
     ierr = 0;
  }
  if ( strcmp( sg_asu_descr, "h>=k and k>=0 and l>=0" ) == 0 &&
       strcmp( sg_patt_group, "4/mmm" ) == 0 ) {
     spacegroup->asufn = &ASU_4_mmm;
     ilaue = ccp4spg_load_laue(spacegroup,8);
     if (strchr(spacegroup->symbol_Hall,'P')) {
       spacegroup->npatt = 123;
       strcpy(spacegroup->patt_name,"P4/mmm");
     } else if (strchr(spacegroup->symbol_Hall,'I')) {
       spacegroup->npatt = 139;
       strcpy(spacegroup->patt_name,"I4/mmm");
     }
     ierr = 0;
  }
  if ( strcmp( sg_asu_descr, "(h>=0 and k>0) or (h=0 and k=0 and l>=0)" ) == 0 ) {
     spacegroup->asufn = &ASU_3b;
     ilaue = ccp4spg_load_laue(spacegroup,9);
     if (strchr(spacegroup->symbol_Hall,'P')) {
       spacegroup->npatt = 147;
       strcpy(spacegroup->patt_name,"P-3");
     } else if (strchr(spacegroup->symbol_xHM,'H')) {
       /* this is special case, as Hall doesn't specify H */
       spacegroup->npatt = 148;
       strcpy(spacegroup->patt_name,"H-3");
     } else if (strchr(spacegroup->symbol_Hall,'R')) {
       spacegroup->npatt = 1148;
       strcpy(spacegroup->patt_name,"R-3");
     }
     ierr = 0;
  }
  if ( strcmp( sg_asu_descr, "h>=k and k>=0 and (k>0 or l>=0)" ) == 0 ) {
     spacegroup->asufn = &ASU_3bm;
     ilaue = ccp4spg_load_laue(spacegroup,10);
     spacegroup->npatt = 162;
     strcpy(spacegroup->patt_name,"P-31m");
     ierr = 0;
  }
  if ( strcmp( sg_asu_descr, "h>=k and k>=0 and (h>k or l>=0)" ) == 0 ) {
     spacegroup->asufn = &ASU_3bmx;
     ilaue = ccp4spg_load_laue(spacegroup,11);
     if (strchr(spacegroup->symbol_Hall,'P')) {
       spacegroup->npatt = 164;
       strcpy(spacegroup->patt_name,"P-3m1");
     } else if (strchr(spacegroup->symbol_xHM,'H')) {
       /* this is special case, as Hall doesn't specify H */
       spacegroup->npatt = 166;
       strcpy(spacegroup->patt_name,"H-3m");
     } else if (strchr(spacegroup->symbol_Hall,'R')) {
       spacegroup->npatt = 1166;
       strcpy(spacegroup->patt_name,"R-3m");
     }
     ierr = 0;
  }
  if ( strcmp( sg_asu_descr, "l>=0 and ((h>=0 and k>0) or (h=0 and k=0))" ) == 0 &&
       strcmp( sg_patt_group, "6/m" ) == 0 ) {
     spacegroup->asufn = &ASU_6_m;
     ilaue = ccp4spg_load_laue(spacegroup,12);
     spacegroup->npatt = 175;
     strcpy(spacegroup->patt_name,"P6/m");
     ierr = 0;
  }
  if ( strcmp( sg_asu_descr, "h>=k and k>=0 and l>=0" ) == 0 &&
       strcmp( sg_patt_group, "6/mmm" ) == 0 ) {
     spacegroup->asufn = &ASU_6_mmm;
     ilaue = ccp4spg_load_laue(spacegroup,13);
     spacegroup->nlaue = 13;
     strcpy(spacegroup->laue_name,"6/mmm");
     spacegroup->laue_sampling[0] = 6;
     spacegroup->laue_sampling[1] = 6;
     spacegroup->laue_sampling[2] = 12;
     spacegroup->npatt = 191;
     strcpy(spacegroup->patt_name,"P6/mmm");
     ierr = 0;
  }
  if ( strcmp( sg_asu_descr, "h>=0 and ((l>=h and k>h) or (l=h and k=h))" ) == 0 ) {
     spacegroup->asufn = &ASU_m3b;
     ilaue = ccp4spg_load_laue(spacegroup,14);
     if (strchr(spacegroup->symbol_Hall,'P')) {
       spacegroup->npatt = 200;
       strcpy(spacegroup->patt_name,"Pm-3");
     } else if (strchr(spacegroup->symbol_Hall,'I')) {
       spacegroup->npatt = 204;
       strcpy(spacegroup->patt_name,"Im-3");
     } else if (strchr(spacegroup->symbol_Hall,'F')) {
       spacegroup->npatt = 202;
       strcpy(spacegroup->patt_name,"Fm-3");
     }
     ierr = 0;
  }
  if ( strcmp( sg_asu_descr, "k>=l and l>=h and h>=0" ) == 0 ) {
     spacegroup->asufn = &ASU_m3bm;
     ilaue = ccp4spg_load_laue(spacegroup,15);
     if (strchr(spacegroup->symbol_Hall,'P')) {
       spacegroup->npatt = 221;
       strcpy(spacegroup->patt_name,"Pm-3m");
     } else if (strchr(spacegroup->symbol_Hall,'I')) {
       spacegroup->npatt = 229;
       strcpy(spacegroup->patt_name,"Im-3m");
     } else if (strchr(spacegroup->symbol_Hall,'F')) {
       spacegroup->npatt = 225;
       strcpy(spacegroup->patt_name,"Fm-3m");
     }
     ierr = 0;
  }

  /* Raise an error if failed to match the ASU description */
  if (ierr) {
    ccp4_signal(CSYM_ERRNO(CSYMERR_NoAsuDefined),"ccp4spg_load_spacegroup",NULL);
    if (spacegroup) free(spacegroup);
    return NULL;
  }

  /* Raise an error if failed to match the Laue code */
  if (ilaue) {
    ccp4_signal(CSYM_ERRNO(CSYMERR_NoLaueCodeDefined),"ccp4spg_load_spacegroup",NULL);
    if (spacegroup) free(spacegroup);
    return NULL;
  }

  /* real asymmetric unit */
  /* origin-based choice */
  sprintf(spacegroup->mapasu_zero_descr,"%s %s %s",map_asu_x,map_asu_y,map_asu_z);
  range_to_limits(map_asu_x, limits);
  spacegroup->mapasu_zero[0] = limits[1];
  range_to_limits(map_asu_y, limits);
  spacegroup->mapasu_zero[1] = limits[1];
  range_to_limits(map_asu_z, limits);
  spacegroup->mapasu_zero[2] = limits[1];
  /* CCP4 choice a la SETLIM - defaults to origin-based choice */
  range_to_limits(map_asu_ccp4_x, limits);
  if (limits[1] > 0) {
    sprintf(spacegroup->mapasu_ccp4_descr,"%s %s %s",map_asu_ccp4_x,map_asu_ccp4_y,map_asu_ccp4_z);
    spacegroup->mapasu_ccp4[0] = limits[1];
    range_to_limits(map_asu_ccp4_y, limits);
    spacegroup->mapasu_ccp4[1] = limits[1];
    range_to_limits(map_asu_ccp4_z, limits);
    spacegroup->mapasu_ccp4[2] = limits[1];
  } else {
    strcpy(spacegroup->mapasu_ccp4_descr,spacegroup->mapasu_zero_descr);
    spacegroup->mapasu_ccp4[0] = spacegroup->mapasu_zero[0];
    spacegroup->mapasu_ccp4[1] = spacegroup->mapasu_zero[1];
    spacegroup->mapasu_ccp4[2] = spacegroup->mapasu_zero[2];
  }
  if (debug) {
    printf(" mapasu limits %f %f %f \n",spacegroup->mapasu_zero[0],
	   spacegroup->mapasu_zero[1],spacegroup->mapasu_zero[2]);
    printf(" CCP4 mapasu limits %f %f %f \n",spacegroup->mapasu_ccp4[0],
	   spacegroup->mapasu_ccp4[1],spacegroup->mapasu_ccp4[2]);
  }

  /* set up centric and epsilon zones for this spacegroup */
  ccp4spg_set_centric_zones(spacegroup);
  ccp4spg_set_epsilon_zones(spacegroup);

  if (debug) 
    printf(" Leaving ccp4spg_load_spacegroup \n");

  return spacegroup;
}


/* symfr_driver

   Convert one or more symop description strings into 4x4 matrix
   representations via multiple calls to symop_to_mat4.

   "line" is a string containing one or more symop description
   strings to be translated into matrix representations.
   Multiple symmetry operations can be specified in a single
   input line, and must be separated by * (with spaces either
   side).
   "rot" is an array of 4x4 matrices in which the symops are
   returned.

   On success, symfr returns the number of symops translated and
   stored; on failure -1 is returned.

   See comments for SYMFR2 for description of the symop formats.
*/
int symfr_driver (const char *line, float rot[MAXSYMOPS][4][4])
{
  CCP4PARSERARRAY *symops=NULL;
  int i,j,k,got_symop=0;
  int ns=0,nsym=0,maxsymops=MAXSYMOPS;
  char *symop=NULL,symopbuf[MAXLENSYMOPSTR];
  float tmp_rot[4][4];

  //CSYMLIB_DEBUG(puts("CSYMLIB: symfr_driver");)

  /* Set up a parser structure to break the line up into
     individual symop strings */
  if ((symops = ccp4_parse_start(maxsymops)) == NULL) {
    /* Couldn't set up a parser structure - abort */
    printf(" symfr_driver: failed to set up parser structure for reading symops.\n");
    return -1;
  }

  /* Tokenise the line, splitting on spaces */
  ccp4_parse_delimiters(symops," ","");
  if ((ns = ccp4_parse(line,symops)) > 0) {

    /* Initialise */
    got_symop = 0;
    symopbuf[0] = '\0';

    /* Loop over tokens and reconstruct symop strings */
    for (i=0; i<ns; ++i) {
      symop = symops->token[i].fullstring;

      /* If there are multiple symop strings then these
	 will be delimited by asterisks */
      if (strlen(symop) == 1 && symop[0] == '*') {
	/* End of symop */
	got_symop = 1;
      } else {
	/* Append token to symop */
	if (strlen(symopbuf)+strlen(symop)+1 <= MAXLENSYMOPSTR) {
	  strcat(symopbuf,symop);
	} else {
	  /* Error - symop string is too long */
	  printf("SYMFR: symmetry operator string is too long!\n");
	  if (symops) ccp4_parse_end(symops);
	  return -1;
	}
	/* Check if this is the last token, in which case
	   flag it to be processed */
	if ((i+1)==ns) got_symop = 1;
      }

      /* Process a complete symop */
      if (got_symop && strlen(symopbuf) > 0) {
	/* Translate */
	if (!symop_to_mat4(&(symopbuf[0]),&(symopbuf[0])+strlen(symopbuf),tmp_rot[0])) {
	  /* Error */
	  if (symops) ccp4_parse_end(symops);
	  return -1;
	}
	/* Load result into the appropriate array location */
	for (j = 0; j < 4; ++j) 
	  for (k = 0; k < 4; ++k) 
	    rot[nsym][j][k] = tmp_rot[j][k];
	nsym++;
	/* Reset for next symop */
	got_symop = 0;
	symopbuf[0] = '\0';
      }
    }
  }

  /* Tidy up and return the number of symops */
  if (symops) ccp4_parse_end(symops);
  return nsym;
}

void ccp4spg_free(CCP4SPG **sp) {
  free ((*sp)->symop);
  free ((*sp)->invsymop);
  free (*sp);
  *sp=NULL;
}

int ccp4_spg_get_centering(const char *symbol_Hall, float cent_ops[4][3])
{
  int debug=0;
  int i,j;

  for (i = 0; i < 4; ++i) 
   for (j = 0; j < 3; ++j) 
    cent_ops[i][j] = 0.0;

  if (strchr(symbol_Hall,'P')) {
    if (debug) printf("Primitive \n");
    return 1;
  } else if (strchr(symbol_Hall,'A')) {
    if (debug) printf("A centering \n");
    cent_ops[1][1] = 0.5;
    cent_ops[1][2] = 0.5;
    return 2;
  } else if (strchr(symbol_Hall,'B')) {
    if (debug) printf("B centering \n");
    cent_ops[1][0] = 0.5;
    cent_ops[1][2] = 0.5;
    return 2;
  } else if (strchr(symbol_Hall,'C')) {
    if (debug) printf("C centering \n");
    cent_ops[1][0] = 0.5;
    cent_ops[1][1] = 0.5;
    return 2;
  } else if (strchr(symbol_Hall,'F')) {
    if (debug) printf("F centering \n");
    cent_ops[1][1] = 0.5;
    cent_ops[1][2] = 0.5;
    cent_ops[2][0] = 0.5;
    cent_ops[2][2] = 0.5;
    cent_ops[3][0] = 0.5;
    cent_ops[3][1] = 0.5;
    return 4;
  } else if (strchr(symbol_Hall,'I')) {
    if (debug) printf("I centering \n");
    cent_ops[1][0] = 0.5;
    cent_ops[1][1] = 0.5;
    cent_ops[1][2] = 0.5;
    return 2;
  } else if (strchr(symbol_Hall,'H')) {
    /* fixme: Hall doesn't specify H, whereas xHM does */
    if (debug) printf("H centering \n");
    cent_ops[1][0] = 2.0/3.0;
    cent_ops[1][1] = 1.0/3.0;
    cent_ops[1][2] = 1.0/3.0;
    cent_ops[2][0] = 1.0/3.0;
    cent_ops[2][1] = 2.0/3.0;
    cent_ops[2][2] = 2.0/3.0;
    return 3;
  } else if (strchr(symbol_Hall,'R')) {
    if (debug) printf("R centering \n");
    return 1;
  }
  return 0;
}

/* standard asu tests for 11 Laue classes */

int ASU_1b   (const int h, const int k, const int l)
  { return (l>0 || (l==0 && (h>0 || (h==0 && k>=0)))); }
int ASU_2_m  (const int h, const int k, const int l)
  { return (k>=0 && (l>0 || (l==0 && h>=0))); }
int ASU_mmm  (const int h, const int k, const int l)
  { return (h>=0 && k>=0 && l>=0); }
int ASU_4_m  (const int h, const int k, const int l)
  { return (l>=0 && ((h>=0 && k>0) || (h==0 && k==0))); }
int ASU_4_mmm(const int h, const int k, const int l)
  { return (h>=k && k>=0 && l>=0); }
int ASU_3b   (const int h, const int k, const int l)
  { return (h>=0 && k>0) || (h==0 && k==0 && l >= 0); }
int ASU_3bm  (const int h, const int k, const int l)
  { return (h>=k && k>=0 && (k>0 || l>=0)); }
int ASU_3bmx (const int h, const int k, const int l)
  { return (h>=k && k>=0 && (h>k || l>=0)); }
int ASU_6_m  (const int h, const int k, const int l)
  { return (l>=0 && ((h>=0 && k>0) || (h==0 && k==0))); }
int ASU_6_mmm(const int h, const int k, const int l)
  { return (h>=k && k>=0 && l>=0); }
int ASU_m3b  (const int h, const int k, const int l)
  { return (h>=0 && ((l>=h && k>h) || (l==h && k==h))); }
int ASU_m3bm  (const int h, const int k, const int l)
  { return (h>=0 && k>=l && l>=h); }

char *ccp4spg_symbol_Hall(CCP4SPG* sp) {
  if (!sp) {  
    ccp4_signal(CSYM_ERRNO(CSYMERR_NullSpacegroup),"ccp4spg_symbol_Hall",NULL); 
    return NULL;
  }
  return sp->symbol_Hall; 
}

ccp4_symop ccp4_symop_invert( const ccp4_symop op1 )
{
  float rot1[4][4],rot2[4][4];

  rotandtrn_to_mat4(rot1,op1);
  invert4matrix((const float (*)[4])rot1,rot2);
  return (mat4_to_rotandtrn((const float (*)[4])rot2));
}

int ccp4spg_name_equal(const char *spgname1, const char *spgname2) {

  char *ch1, *ch2, *spgname1_upper, *spgname2_upper;

  /* create copies of input strings, and convert to upper case */
  spgname1_upper = strdup(spgname1);
  strtoupper(spgname1_upper,spgname1);
  spgname2_upper = strdup(spgname2);
  strtoupper(spgname2_upper,spgname2);

  ch1 = spgname1_upper;
  ch2 = spgname2_upper;
  while (*ch1 == *ch2) {
    if (*ch1 == '\0' && *ch2 == '\0') {
      free(spgname1_upper);
      free(spgname2_upper);
      return 1;
    }
    ++ch1;
    ++ch2;
  }
  free(spgname1_upper);
  free(spgname2_upper);
  return 0;
}

int ccp4spg_name_equal_to_lib(const char *spgname_lib, const char *spgname_match) {

  char *ch1, *ch2, *spgname1_upper, *spgname2_upper, *tmpstr;
  int have_one_1=0, have_one_2=0;

  /* create copies of input strings, convert to upper case, and
     deal with colons */
  spgname1_upper = strdup(spgname_lib);
  strtoupper(spgname1_upper,spgname_lib);
  ccp4spg_name_de_colon(spgname1_upper);
  spgname2_upper = strdup(spgname_match);
  strtoupper(spgname2_upper,spgname_match);
  ccp4spg_name_de_colon(spgname2_upper);

  /* see if strings are equal, except for spaces */
  ch1 = spgname1_upper;
  ch2 = spgname2_upper;
  while (*ch1 == *ch2) {
    if (*ch1 == '\0' && *ch2 == '\0') {
      free(spgname1_upper);
      free(spgname2_upper);
      return 1;
    }
    while (*(++ch1) == ' ') ;
    while (*(++ch2) == ' ') ;
  }

  /* if that didn't work, and spgname_match is a short name, try removing
     " 1 " from spgname_lib, and matching again. This would match P21 to
     'P 1 21 1' for instance. */

  /* try to identify if "short names" are being used. */
  if (strstr(spgname1_upper," 1 ")) have_one_1 = 1;
  if (strstr(spgname2_upper," 1 ")) have_one_2 = 1; 
  /* if spgname_lib has " 1 " and spgname_match doesn't, then strip
     out " 1" to do "short" comparison */
  if (have_one_1 && ! have_one_2) {
     tmpstr = strdup(spgname1_upper);
     ccp4spg_to_shortname(tmpstr,spgname1_upper);
     strcpy(spgname1_upper,tmpstr);
     free(tmpstr);
  }

  /* see if strings are equal, except for spaces */
  ch1 = spgname1_upper;
  ch2 = spgname2_upper;
  while (*ch1 == *ch2) {
    if (*ch1 == '\0' && *ch2 == '\0') {
      free(spgname1_upper);
      free(spgname2_upper);
      return 1;
    }
    while (*(++ch1) == ' ') ;
    while (*(++ch2) == ' ') ;
  }

  free(spgname1_upper);
  free(spgname2_upper);
  return 0;
}

char *ccp4spg_to_shortname(char *shortname, const char *longname) {

  const char *ch1;
  char *ch2;
  int trigonal=0;

  ch1 = longname;
  ch2 = shortname;

  /* "P 1" is an exception */
  if (!strcmp(ch1,"P 1")) {
    strcpy(ch2,"P1");
    return  ch2;
  }

  /* trigonal are another exception, don't want to lose significant " 1" */
  if (!strncmp(ch1,"P 3",3) || !strncmp(ch1,"P -3",4) || !strncmp(ch1,"R 3",3) || !strncmp(ch1,"R -3",4)) trigonal=1;

  while (*ch1 != '\0') {
    if (!trigonal && !strncmp(ch1," 1",2)) {
      ch1 += 2;
    } else {
      /* take out blanks - note check for " 1" takes precedence */
      while (*ch1 == ' ') ++ch1;
      if (*ch1 != '\0') { *ch2 = *ch1; ++ch2; ++ch1; }
    }
  }
  *ch2 = '\0';
  return ch2;
}

void ccp4spg_name_de_colon(char *name) {

  char *ch1;

  /* various spacegroup names have settings specified by colon. We'll
     deal with these on a case-by-case basis. */
  if ((ch1 = strstr(name,":R")) != NULL) {
  /* :R spacegroup should be R already so just replace with blanks */
    *ch1 = ' ';
    *(ch1+1) = ' ';
  } else if ((ch1 = strstr(name,":H")) != NULL) {
  /* :H spacegroup should be R so change to H */
    *ch1 = ' ';
    *(ch1+1) = ' ';
    ch1 = strstr(name,"R");
    if (ch1 != NULL) *ch1 = 'H';
  }
    
  return;
}

int ccp4spg_pgname_equal(const char *pgname1, const char *pgname2) {

  char *ch1, *ch2, *pgname1_upper, *pgname2_upper;

  pgname1_upper = strdup(pgname1);
  strtoupper(pgname1_upper,pgname1);
  pgname2_upper = strdup(pgname2);
  strtoupper(pgname2_upper,pgname2);

  ch1 = pgname1_upper;
  if (pgname1_upper[0] == 'P' && pgname1_upper[1] == 'G') ch1 += 2;
  ch2 = pgname2_upper;
  if (pgname2_upper[0] == 'P' && pgname2_upper[1] == 'G') ch2 += 2;
  while (*ch1 == *ch2) {
    if (*ch1 == '\0' && *ch2 == '\0') {
      free(pgname1_upper);
      free(pgname2_upper);
      return 1;
    }
    while (*(++ch1) == ' ') ;
    while (*(++ch2) == ' ') ;
  }
  free(pgname1_upper);
  free(pgname2_upper);
  return 0;
}

ccp4_symop *ccp4spg_norm_trans(ccp4_symop *op) {

  int i;

  for ( i = 0; i < 3; i++ ) {
    while (op->trn[i] < 0.0) op->trn[i] += 1.0;
    while (op->trn[i] >= 1.0) op->trn[i] -= 1.0;
  }

  return op;
}

int ccp4_spgrp_equal( int nsym1, const ccp4_symop *op1, int nsym2, const ccp4_symop *op2 )
{
  int i, n;
  int *symcode1, *symcode2;

  /* first check that we have equal number of symops */
  if ( nsym1 != nsym2 ) return 0;

  n = nsym1;

  /* now make the sym code arrays */
  symcode1 = ccp4_utils_malloc( n * sizeof(int) );
  symcode2 = ccp4_utils_malloc( n * sizeof(int) );
  for ( i = 0; i < n; i++ ) {
    symcode1[i] = ccp4_symop_code( op1[i] );
    symcode2[i] = ccp4_symop_code( op2[i] );
  }
  /* sort the symcodes */
  /* Kevin suggests maybe just compare all pairs rather than sort */
  qsort( symcode1, n, sizeof(int), &ccp4_int_compare );
  qsort( symcode2, n, sizeof(int), &ccp4_int_compare );
  /* compare the symcodes */
  for ( i = 0; i < n; i++ ) {
    if ( symcode1[i] != symcode2[i] ) break;
  }
  /* delete the symcodes */
  free(symcode1);
  free(symcode2);

  /* return true if they are equal */
  return ( i == n );
}

int ccp4_spgrp_equal_order( int nsym1, const ccp4_symop *op1, int nsym2, const ccp4_symop *op2 )
{
  int i;

  /* first check that we have equal number of symops */
  if ( nsym1 != nsym2 ) return 0;

  /* compare the symcodes */
  for ( i = 0; i < nsym1; i++ ) {
    if ( ccp4_symop_code( op1[i] ) != ccp4_symop_code( op2[i] ) ) break;
  }

  /* return true if they are equal */
  return ( i == nsym1 );
}

int ccp4_symop_code(ccp4_symop op)
{
  int i, j, code=0;
  for ( i=0; i<3; i++ )
    for ( j=0; j<3; j++ )
      code = ( code << 2 ) | ( (int) rint( op.rot[i][j] ) & 0x03 ) ;
  for ( i=0; i<3; i++ )
    code = ( code << 4 ) | ( (int) rint( op.trn[i]*12.0 ) & 0x0f ) ;
  return code;
}

int ccp4_int_compare( const void *p1, const void *p2 )
{
  return ( *((int*)p1) - *((int*)p2) );
}

int ccp4spg_is_in_pm_asu(const CCP4SPG* sp, const int h, const int k, const int l) {
  if (ccp4spg_is_in_asu(sp,h,k,l)) return (1);
  if (ccp4spg_is_in_asu(sp,-h,-k,-l)) return (-1);
  return 0;
}

int ccp4spg_is_in_asu(const CCP4SPG* sp, const int h, const int k, const int l) {
  if (!sp) {  
    ccp4_signal(CSYM_ERRNO(CSYMERR_NullSpacegroup),"ccp4spg_is_in_asu",NULL); 
    return 0;
  }
  if ( ccp4spg_do_chb(sp->chb) ) return sp->asufn(
      (int) rint( h*sp->chb[0][0] + k*sp->chb[1][0] + l*sp->chb[2][0] ),
      (int) rint( h*sp->chb[0][1] + k*sp->chb[1][1] + l*sp->chb[2][1] ),
      (int) rint( h*sp->chb[0][2] + k*sp->chb[1][2] + l*sp->chb[2][2] ) );
  else
    return sp->asufn( h, k, l );
}

int ccp4spg_put_in_asu(const CCP4SPG* sp, const int hin, const int kin, const int lin,
		       int *hout, int *kout, int *lout ) {

  int i, isign;

  if (!sp) {  
    ccp4_signal(CSYM_ERRNO(CSYMERR_NullSpacegroup),"ccp4spg_put_in_asu",NULL); 
    return 0;
  }

  /* cycle through all primitive symmetry operations until in asu */

  for (i = 0; i < sp->nsymop_prim; ++i) {
    *hout = (int) rint( hin*sp->symop[i].rot[0][0] + kin*sp->symop[i].rot[1][0] + 
                        lin*sp->symop[i].rot[2][0] ); 
    *kout = (int) rint( hin*sp->symop[i].rot[0][1] + kin*sp->symop[i].rot[1][1] + 
                        lin*sp->symop[i].rot[2][1] ); 
    *lout = (int) rint( hin*sp->symop[i].rot[0][2] + kin*sp->symop[i].rot[1][2] + 
                        lin*sp->symop[i].rot[2][2] ); 
    if ((isign = ccp4spg_is_in_pm_asu(sp,*hout,*kout,*lout)) != 0) {
      *hout = *hout * isign;
      *kout = *kout * isign;
      *lout = *lout * isign;
      return ( (isign > 0) ? 2*i+1 : 2*i+2 );
    }
  }

  printf ("Can't put in asu ! \n");
  return 0;
}

/* Generate indices according to symmetry operation isym */

void ccp4spg_generate_indices(const CCP4SPG* sp, const int isym,
                  const int hin, const int kin, const int lin,
		       int *hout, int *kout, int *lout ) {

  int jsym, isign;

  if (!sp) {  
    ccp4_signal(CSYM_ERRNO(CSYMERR_NullSpacegroup),"ccp4spg_generate_indices",NULL); 
    return;
  }

  jsym = (isym - 1) / 2;
  isign = (isym % 2) ? 1 : -1 ; 

  *hout = isign * (int) rint(hin*sp->invsymop[jsym].rot[0][0] + 
       kin*sp->invsymop[jsym].rot[1][0] + lin*sp->invsymop[jsym].rot[2][0]); 
  *kout = isign * (int) rint(hin*sp->invsymop[jsym].rot[0][1] + 
       kin*sp->invsymop[jsym].rot[1][1] + lin*sp->invsymop[jsym].rot[2][1]); 
  *lout = isign * (int) rint(hin*sp->invsymop[jsym].rot[0][2] + 
       kin*sp->invsymop[jsym].rot[1][2] + lin*sp->invsymop[jsym].rot[2][2]); 

}

/* shift phase value associated with hin,kin,lin according to translation 
and optional sign change. Return in range 0,360 */

float ccp4spg_phase_shift(const int hin, const int kin, const int lin,
                const float phasin, const float trans[3], const int isign)
{
  double phasout;

  phasout = (double) phasin;
  if (isign == -1) phasout = - phasout;

  phasout += (hin*trans[0] + kin*trans[1] + lin*trans[2]) * 360.0;

  phasout = fmod(phasout,360.0);
  if (phasout < 0.0) phasout += 360.0;

  return ((float) phasout);

}

int ccp4spg_do_chb(const float chb[3][3]) {

  return ( chb[0][0] != 1 || chb[1][1] != 1 || chb[2][2] != 1 ||
	   chb[0][1] != 0 || chb[0][2] != 0 || chb[1][2] != 0 ||
	   chb[1][0] != 0 || chb[2][0] != 0 || chb[2][1] != 0 );

}

/* functions to identify centrics - based on Randy's method */

void ccp4spg_set_centric_zones(CCP4SPG* sp) {

  int i,j,hnew,knew,lnew;
  int ihkl[12][3];

  if (!sp) {  
    ccp4_signal(CSYM_ERRNO(CSYMERR_NullSpacegroup),"ccp4spg_set_centric_zones",NULL); 
    return;
  }
  ihkl[0][0] = 0; ihkl[0][1] = 1; ihkl[0][2] = 2; 
  ihkl[1][0] = 1; ihkl[1][1] = 0; ihkl[1][2] = 2; 
  ihkl[2][0] = 1; ihkl[2][1] = 2; ihkl[2][2] = 0; 
  ihkl[3][0] = 1; ihkl[3][1] = 1; ihkl[3][2] = 10; 
  ihkl[4][0] = 1; ihkl[4][1] = 10; ihkl[4][2] = 1; 
  ihkl[5][0] = 10; ihkl[5][1] = 1; ihkl[5][2] = 1; 
  ihkl[6][0] = 1; ihkl[6][1] = -1; ihkl[6][2] = 10; 
  ihkl[7][0] = 1; ihkl[7][1] = 10; ihkl[7][2] = -1; 
  ihkl[8][0] = 10; ihkl[8][1] = 1; ihkl[8][2] = -1; 
  ihkl[9][0] = -1; ihkl[9][1] = 2; ihkl[9][2] = 10; 
  ihkl[10][0] = 2; ihkl[10][1] = -1; ihkl[10][2] = 10; 
  ihkl[11][0] = 1; ihkl[11][1] = 4; ihkl[11][2] = 8; 

  /* loop over all possible centric zones */
  for (i = 0; i < 12; ++i) {
   sp->centrics[i] = 0;
   for (j = 0; j < sp->nsymop; ++j) {
    hnew = (int) rint( ihkl[i][0]*sp->symop[j].rot[0][0] + 
      ihkl[i][1]*sp->symop[j].rot[1][0] + ihkl[i][2]*sp->symop[j].rot[2][0] );
    if (hnew == -ihkl[i][0]) {
     knew = (int) rint( ihkl[i][0]*sp->symop[j].rot[0][1] + 
       ihkl[i][1]*sp->symop[j].rot[1][1] + ihkl[i][2]*sp->symop[j].rot[2][1] );
     if (knew == -ihkl[i][1]) {
      lnew = (int) rint( ihkl[i][0]*sp->symop[j].rot[0][2] + 
        ihkl[i][1]*sp->symop[j].rot[1][2] + ihkl[i][2]*sp->symop[j].rot[2][2] );
      if (lnew == -ihkl[i][2]) {
        sp->centrics[i] = j+1;
        break;
      }
     }
    }
   }
  }
}

int ccp4spg_is_centric(const CCP4SPG* sp, const int h, const int k, const int l) {

  int i;

  if (!sp) {  
    ccp4_signal(CSYM_ERRNO(CSYMERR_NullSpacegroup),"ccp4spg_is_centric",NULL); 
    return -1;
  }
  /* loop over all possible centric zones */
  for (i = 0; i < 12; ++i) 
    if (sp->centrics[i]) 
      if (ccp4spg_check_centric_zone(i+1,h,k,l) == 0)
        return 1;

  return 0;
}

/* check indices against centric zones - return 0 if in zone "nzone" */

int ccp4spg_check_centric_zone(const int nzone, const int h, const int k, const int l) {

  switch (nzone) {
  case 1:
    return h;
  case 2:
    return k;
  case 3:
    return l;
  case 4:
    return h - k;
  case 5:
    return h - l;
  case 6:
    return k - l;
  case 7:
    return h + k;
  case 8:
    return h + l;
  case 9:
    return k + l;
  case 10:
    return 2*h + k;
  case 11:
    return h + 2*k;
  case 12:
    return 0;
  }
  printf ("Invalid nzone ! \n");
  return 0;
}

float ccp4spg_centric_phase(const CCP4SPG* sp, const int h, const int k, const int l) {

  int i,isym;
  float centric_phase;

  if (!sp) {  
    ccp4_signal(CSYM_ERRNO(CSYMERR_NullSpacegroup),"ccp4spg_centric_phase",NULL); 
    return 0.0;
  }
  /* loop over all possible centric zones */
  for (i = 0; i < 12; ++i) 
    if (sp->centrics[i]) 
      if (ccp4spg_check_centric_zone(i+1,h,k,l) == 0) {
        isym = sp->centrics[i];
        centric_phase = h*sp->symop[isym-1].trn[0] + 
	  k*sp->symop[isym-1].trn[1] + l*sp->symop[isym-1].trn[2];
        centric_phase = 180.0*(centric_phase - rint(centric_phase));
        if (centric_phase < 0.0) centric_phase = centric_phase + 180.0;
        return centric_phase;
      }

  return 0;
}

void ccp4spg_print_centric_zones(const CCP4SPG* sp) {

  int i,j=0;
  char centric_zone[8];

  if (!sp) {  
    ccp4_signal(CSYM_ERRNO(CSYMERR_NullSpacegroup),"ccp4spg_print_centric_zones",NULL); 
    return;
  }
  printf("\n  ******   CENTRIC ZONES  ****** \n\n");

  /* loop over all possible centric zones */
  for (i = 0; i < 12; ++i) 
    if (sp->centrics[i]) {
      printf("  CENTRIC Zone   %d\n",++j);
      printf("  Reflections of type  %s \n\n",
                ccp4spg_describe_centric_zone(i+1,centric_zone));
    }

  if (!j) printf("  (no centric zones for this spacegroup) \n\n");
}

char *ccp4spg_describe_centric_zone(const int nzone, char *centric_zone) {

  switch (nzone) {
  case 1:
    return ( strcpy(centric_zone,"0kl") );
  case 2:
    return ( strcpy(centric_zone,"h0l") );
  case 3:
    return ( strcpy(centric_zone,"hk0") );
  case 4:
    return ( strcpy(centric_zone,"hhl") );
  case 5:
    return ( strcpy(centric_zone,"hkh") );
  case 6:
    return ( strcpy(centric_zone,"hkk") );
  case 7:
    return ( strcpy(centric_zone,"h -hl") );
  case 8:
    return ( strcpy(centric_zone,"hk -h") );
  case 9:
    return ( strcpy(centric_zone,"hk -k") );
  case 10:
    return ( strcpy(centric_zone,"-h 2h l") );
  case 11:
    return ( strcpy(centric_zone,"2h -h l") );
  case 12:
    return ( strcpy(centric_zone,"hkl") );
  }
  printf ("Invalid nzone ! \n");
  return "null";
}

/* functions to identify epsilon zones - based on Randy's method */

void ccp4spg_set_epsilon_zones(CCP4SPG* sp) {

  int i,j,hnew,knew,lnew,neps;
  int ihkl[13][3];

  if (!sp) {  
    ccp4_signal(CSYM_ERRNO(CSYMERR_NullSpacegroup),"ccp4spg_set_epsilon_zones",NULL); 
    return;
  }
  ihkl[0][0] = 1; ihkl[0][1] = 0; ihkl[0][2] = 0; 
  ihkl[1][0] = 0; ihkl[1][1] = 2; ihkl[1][2] = 0; 
  ihkl[2][0] = 0; ihkl[2][1] = 0; ihkl[2][2] = 2; 
  ihkl[3][0] = 1; ihkl[3][1] = 1; ihkl[3][2] = 0; 
  ihkl[4][0] = 1; ihkl[4][1] = 0; ihkl[4][2] = 1; 
  ihkl[5][0] = 0; ihkl[5][1] = 1; ihkl[5][2] = 1; 
  ihkl[6][0] = 1; ihkl[6][1] = -1; ihkl[6][2] = 0; 
  ihkl[7][0] = 1; ihkl[7][1] = 0; ihkl[7][2] = -1; 
  ihkl[8][0] = 0; ihkl[8][1] = 1; ihkl[8][2] = -1; 
  ihkl[9][0] = -1; ihkl[9][1] = 2; ihkl[9][2] = 0; 
  ihkl[10][0] = 2; ihkl[10][1] = -1; ihkl[10][2] = 0; 
  ihkl[11][0] = 1; ihkl[11][1] = 1; ihkl[11][2] = 1; 
  ihkl[12][0] = 1; ihkl[12][1] = 2; ihkl[12][2] = 3; 

  /* Loop over all possible epsilon zones, except the catch-all 13th. For each 
     zone, "neps" counts the number of symmetry operators that map a representative
     reflection "ihkl" to itself. At least the identity will do this. If any
     more do, then this is a relevant epsilon zone. */
  for (i = 0; i < 12; ++i) {
   sp->epsilon[i] = 0;
   neps = 0;
   for (j = 0; j < sp->nsymop_prim; ++j) {
    hnew = (int) rint( ihkl[i][0]*sp->symop[j].rot[0][0] + 
      ihkl[i][1]*sp->symop[j].rot[1][0] + ihkl[i][2]*sp->symop[j].rot[2][0] );
    if (hnew == ihkl[i][0]) {
     knew = (int) rint( ihkl[i][0]*sp->symop[j].rot[0][1] + 
       ihkl[i][1]*sp->symop[j].rot[1][1] + ihkl[i][2]*sp->symop[j].rot[2][1] );
     if (knew == ihkl[i][1]) {
      lnew = (int) rint( ihkl[i][0]*sp->symop[j].rot[0][2] + 
        ihkl[i][1]*sp->symop[j].rot[1][2] + ihkl[i][2]*sp->symop[j].rot[2][2] );
      if (lnew == ihkl[i][2]) {
        ++neps;
      }
     }
    }
   }
   if (neps > 1) sp->epsilon[i] = neps * (sp->nsymop/sp->nsymop_prim);
  }
  /* hkl zone covers all with neps of 1 */
  sp->epsilon[12] = sp->nsymop/sp->nsymop_prim;
}

int ccp4spg_get_multiplicity(const CCP4SPG* sp, const int h, const int k, const int l) {

  int i;

  if (!sp) {  
    ccp4_signal(CSYM_ERRNO(CSYMERR_NullSpacegroup),"ccp4spg_get_multiplicity",NULL); 
    return 0;
  }
  /* loop over all possible epsilon zones */
  for (i = 0; i < 13; ++i) 
    if (sp->epsilon[i]) 
      if (ccp4spg_check_epsilon_zone(i+1,h,k,l) == 0)
        return sp->epsilon[i];

  return 0;
}

int ccp4spg_check_epsilon_zone(const int nzone, const int h, const int k, const int l) {

  int bigfac=1000;  /* this needs to be big enough to prevent accidental zeros */

  switch (nzone) {
  case 1:
    return bigfac*k + l;
  case 2:
    return h + bigfac*l;
  case 3:
    return h + bigfac*k;
  case 4:
    return h - k + bigfac*l;
  case 5:
    return h + bigfac*k - l;
  case 6:
    return bigfac*h + k - l;
  case 7:
    return h + k + bigfac*l;
  case 8:
    return h + bigfac*k + l;
  case 9:
    return bigfac*h + k + l;
  case 10:
    return 2*h + k + bigfac*l;
  case 11:
    return h + 2*k + bigfac*l;
  case 12:
    return h + bigfac*k - (bigfac+1)*l;
  case 13:
    return 0;
  }
  printf ("Invalid nzone ! \n");
  return 0;
}

void ccp4spg_print_epsilon_zones(const CCP4SPG* sp) {

  int i,j=0;
  char epsilon_zone[8];

  if (!sp) {  
    ccp4_signal(CSYM_ERRNO(CSYMERR_NullSpacegroup),"ccp4spg_print_epsilon_zones",NULL); 
    return;
  }
  printf("\n  ******   EPSILON ZONES -  Reflection Classes and their multiplicity ****** \n");

  /* loop over all possible epsilon zones */
  for (i = 0; i < 13; ++i) 
    if (sp->epsilon[i]) {
      printf("\n  EPSILON Zone   %d\n",++j);
      printf("  Reflections of type  %s \n",
                ccp4spg_describe_epsilon_zone(i+1,epsilon_zone));
      printf("  Multiplicity   %d\n",sp->epsilon[i]);
    }
}

char *ccp4spg_describe_epsilon_zone(const int nzone, char *epsilon_zone) {

  switch (nzone) {
  case 1:
    return ( strcpy(epsilon_zone,"h00") );
  case 2:
    return ( strcpy(epsilon_zone,"0k0") );
  case 3:
    return ( strcpy(epsilon_zone,"00l") );
  case 4:
    return ( strcpy(epsilon_zone,"hh0") );
  case 5:
    return ( strcpy(epsilon_zone,"h0h") );
  case 6:
    return ( strcpy(epsilon_zone,"0kk") );
  case 7:
    return ( strcpy(epsilon_zone,"h -h0") );
  case 8:
    return ( strcpy(epsilon_zone,"h0 -h") );
  case 9:
    return ( strcpy(epsilon_zone,"0k -k") );
  case 10:
    return ( strcpy(epsilon_zone,"-h 2h 0") );
  case 11:
    return ( strcpy(epsilon_zone,"2h -h 0") );
  case 12:
    return ( strcpy(epsilon_zone,"hhh") );
  case 13:
    return ( strcpy(epsilon_zone,"hkl") );
  }
  printf ("Invalid nzone ! \n");
  return "null";
}

int ccp4spg_is_sysabs(const CCP4SPG* sp, const int h, const int k, const int l)
{
  int j,hnew,knew,lnew;
  float del_phas;

  if (!sp) {  
    ccp4_signal(CSYM_ERRNO(CSYMERR_NullSpacegroup),"ccp4spg_is_sysabs",NULL); 
    return -1;
  }
  if (sp->nsymop > 1) {
   for (j = 1; j < sp->nsymop; ++j) {
    hnew = (int) rint( h*sp->invsymop[j].rot[0][0] + k*sp->invsymop[j].rot[1][0] + 
      l*sp->invsymop[j].rot[2][0] );
    if (hnew == h) {
     knew = (int) rint( h*sp->invsymop[j].rot[0][1] + k*sp->invsymop[j].rot[1][1] + 
       l*sp->invsymop[j].rot[2][1] );
     if (knew == k) {
      lnew = (int) rint( h*sp->invsymop[j].rot[0][2] + k*sp->invsymop[j].rot[1][2] + 
	l*sp->invsymop[j].rot[2][2] );
      if (lnew == l) {
	/* phase shift from translational component of sym op */
        del_phas = h*sp->symop[j].trn[0] + k*sp->symop[j].trn[1] + 
                   l*sp->symop[j].trn[2];
        if ( fabs(del_phas - rint( del_phas )) > 0.05 ) return (1);
      }
     }
    }
   }
  }
  return (0);

}

int ccp4spg_generate_origins(const char *namspg, const int nsym, const float rsym[][4][4],
			     float origins[][3], int *polarx, int *polary, int *polarz,
			     const int iprint)
{
  int i,j,k,norigins,k1,k2,k3,alt_orig,ichk;
  int id[6]={0,6,4,8,3,9},is[3];
  float xin=0.13,yin=0.17,zin=0.19,xout,yout,zout,rsymd[3][3];

  *polarx = *polary = *polarz = 1;

  for (i = 1; i < nsym; ++i) {
    xout = rsym[i][0][0]*xin + rsym[i][0][1]*yin + rsym[i][0][2]*zin;
    if (fabs(xout-xin) > 0.01) *polarx = 0;
    yout = rsym[i][1][0]*xin + rsym[i][1][1]*yin + rsym[i][1][2]*zin;
    if (fabs(yout-yin) > 0.01) *polary = 0;
    zout = rsym[i][2][0]*xin + rsym[i][2][1]*yin + rsym[i][2][2]*zin;
    if (fabs(zout-zin) > 0.01) *polarz = 0;
  }

  /*  First origin is 0,0,0 */
  norigins=1;
  origins[0][0]=0.0;
  origins[0][1]=0.0;
  origins[0][2]=0.0;

      /*  check which points can be an alternate origin.
	  only six possibilities which are 0 1/2 1/3 2/3 1/4 3/4 
          is/id expressed as twelfths */
  for (k1 = 0; k1 < 6; ++k1) {
    for (k2 = 0; k2 < 6; ++k2) {
      for (k3 = 0; k3 < 6; ++k3) {
        if (k1==0 && k2 == 0 && k3 ==0) continue;
	is[0]=id[k1];
	is[1]=id[k2];
	is[2]=id[k3];
        if ( *polarx && is[0] )  continue;
        if ( *polary && is[1] )  continue;
        if ( *polarz && is[2] )  continue;

/*  Let [Si] =[RSYMi] be (3x4) symmetry operator.
 Need to Check if the symmetry operator shifts of each alternate origin 
 [ORx,ORy,ORz)  are preserved for each symmetry operator.
 Origin (0,0,0) shifts to        Ti(1),     Ti(2)      Ti(3) 
                           == RSYMi(1,4),RSYMi(2,4),RSYMi(3,4) 

 [RSYMi] [OR]  =  [OR] + [Ti] + n[I]  = [1 0 0 RSYMi(1,4)] [OR1] +  n[I]
                                       [0 1 0 RSYMi(2,4)] [OR2]
                                       [0 0 1 RSYMi(3,4)] [OR3]

 Hence [RSYMi(1,1) -1   RSYMi(1,2)      RSYMi(1,3)      0] [OR1]   = n[I]
       [RSYMi(2,1)      RSYMi(2,2) -1   RSYMi(2,3)      0] [OR2] 
       [RSYMi(3,1)      RSYMi(3,2)      RSYMi(3,3) -1   0] [OR3] 
       [   0                0               0           1] [1  ]

 Use RSYM(..1) to respresent indentity.. Enough to use 3x3 matrix.. */

        alt_orig=1;
	for (i = 1; i < nsym && alt_orig; ++i) {
	  for (j = 0; j < 3; ++j) 
	    for (k = 0; k < 3; ++k) 
              rsymd[j][k] = rsym[i][j][k] - rsym[0][j][k];
	  for (j = 0; j < 3; ++j) {
            ichk = (int) rint( rsymd[j][0]*is[0]+rsymd[j][1]*is[1]+rsymd[j][2]*is[2] );
            if ( ichk % 12 ) {
              alt_orig=0;
              break;
	    }
          }
        }
        if (alt_orig) {
         norigins+=1;
         origins[norigins-1][0]=is[0]/12.0;
         origins[norigins-1][1]=is[1]/12.0;
         origins[norigins-1][2]=is[2]/12.0;
        }
      }
    }
  }

  if (iprint) {
    if( *polarx && *polary && *polarz) {
      printf(" this is p1: origin anywhere");
      printf("\n %s %s %s \n",
	     "Number of alternate origins for spacegroup:  ",namspg," is infinite.");
    } else if( *polarx && *polary) {
      printf(" this is a polar+ spacegroup: origin anywhere in a b plane");
      printf("\n %s %s %s %d \n",
     "Number of alternate origin containing planes for spacegroup:",
       namspg, " is:",norigins);
    } else if( *polarx && *polarz) {
      printf(" this is a polar+ spacegroup: origin anywhere in a c plane");
      printf("\n %s %s %s %d \n", 
     "Number of alternate origin containing planes for spacegroup:",
       namspg, " is:",norigins);
    } else if( *polary && *polarz) {
      printf(" this is a polar+ spacegroup: origin anywhere in b c plane");
      printf("\n %s %s %s %d \n", 
     "Number of alternate origin containing planes for spacegroup:",
       namspg, " is:",norigins);
    } else if( *polarx) {
      printf(" this is a polar spacegroup: origin is not fixed along a axis");
      printf("\n %s %s %s %d \n", 
     "Number of alternate origin containing lines for spacegroup: ",
       namspg, " is:",norigins);
    } else if( *polary) {
      printf(" this is a polar spacegroup: origin is not fixed along b axis");
      printf("\n %s %s %s %d \n", 
     "Number of alternate origin containing lines for spacegroup: ",
       namspg, " is:",norigins);
    } else if( *polarz) {
      printf(" this is a polar spacegroup: origin is not fixed along c axis");
      printf("\n %s %s %s %d \n", 
     "Number of alternate origin containing lines for spacegroup: ",
       namspg, " is:",norigins);
    } else {
      printf("\n %s %s %s %d \n",
     "Number of alternate origins for spacegroup:  ",namspg,
     " is:",norigins);
    } 

    printf("\n Norigin     Ox      Oy      Oz\n\n");
    for (i = 0; i < norigins; ++i) {
      if (*polary && *polarz && *polarx) {
        printf("%8d     ??      ??      ?? \n", i+1);
      } else if(*polarx && *polary) {
        printf("%8d     ??      ?? %8.4f\n", i+1,origins[i][2]);
      } else if(*polarx && *polarz) {
        printf("%8d     ?? %8.4f     ?? \n", i+1,origins[i][1]);
      } else if(*polary && *polarz) {
        printf("%8d%8.4f     ??      ?? \n", i+1,origins[i][0]);
      } else if( *polarx) {
        printf("%8d     ?? %8.4f%8.4f\n", i+1,origins[i][1],origins[i][2]);
      } else if(*polary) {
        printf("%8d%8.4f     ?? %8.4f\n", i+1,origins[i][0],origins[i][2]);
      } else if(*polarz) {
        printf("%8d%8.4f%8.4f     ?? \n", i+1,origins[i][0],origins[i][1]);
      } else  {
        printf("%8d%8.4f%8.4f%8.4f\n", i+1,origins[i][0],origins[i][1],origins[i][2]);
      }
    }
  }
  return norigins;
}

void ccp4spg_print_recip_spgrp(const CCP4SPG* sp)
{
  if (!sp) {  
    ccp4_signal(CSYM_ERRNO(CSYMERR_NullSpacegroup),"ccp4spg_print_recip_spgrp",NULL); 
    return;
  }

    printf("Reciprocal space symmetry: \n");
    printf("Space group: \"%s\" Point group: \"%s\" Laue group: \"%s\" \n",
       sp->symbol_xHM,sp->point_group,sp->laue_name); 
    printf("Reference asymmetric unit: \"%s\" \n",sp->asu_descr); 
    printf("  (change of basis may be applied) \n");
    ccp4spg_print_recip_ops(sp);
}

void ccp4spg_print_recip_ops(const CCP4SPG* sp)
{
  int i,j,k,l,nrow, n_in_last_row,rsymop_len=80;
  float tmp_symop[4][4];
  char rsymop[80];

  if (!sp) {  
    ccp4_signal(CSYM_ERRNO(CSYMERR_NullSpacegroup),"ccp4spg_print_recip_ops",NULL); 
    return;
  }

  nrow = (sp->nsymop_prim + 3)/ 4;
  n_in_last_row = sp->nsymop_prim % 4;
  if (n_in_last_row == 0) n_in_last_row = 4;

  printf("\n Spacegroup %d \"%s\" \n",sp->spg_ccp4_num,sp->symbol_xHM);
  printf(" Original indices for reflection hkl with symmetry number ISYM \n");
  printf("\n                              Bijvoet positive \n");
  printf("       %-18s%-18s%-18s%-18s\n","ISYM","ISYM","ISYM","ISYM");
  for (i = 0 ; i < nrow-1 ; ++i) {
    printf("  ISYM");
    for (j = 0 ; j < 4 ; ++j) {
      for (k = 0; k < 3; ++k) {
	/* note we use the transpose for reciprocal space operators */
        for (l = 0; l < 3; ++l) 
          tmp_symop[k][l] = sp->invsymop[4*i+j].rot[l][k];
        tmp_symop[k][3] = 0.0;
        tmp_symop[3][k] = 0.0;
      }
      tmp_symop[3][3] = 1.0;
      mat4_to_recip_symop(rsymop,rsymop+rsymop_len,(const float (*)[4])tmp_symop);
      rsymop[12] = '\0';
      printf(" %3d  %-12s",2*(4*i+j)+1,rsymop);
    }
    printf("\n");
  }
  printf("  ISYM");
  for (j = 0 ; j < n_in_last_row ; ++j) {
    for (k = 0; k < 3; ++k) {
      for (l = 0; l < 3; ++l) 
        tmp_symop[k][l] = sp->invsymop[4*i+j].rot[l][k];
      tmp_symop[k][3] = 0.0;
      tmp_symop[3][k] = 0.0;
    }
    tmp_symop[3][3] = 1.0;
    mat4_to_recip_symop(rsymop,rsymop+rsymop_len,(const float (*)[4])tmp_symop);
    rsymop[12] = '\0';
    printf(" %3d  %-12s",2*(4*(nrow-1)+j)+1,rsymop);
  }
  printf("\n");

  printf("\n                              Bijvoet negative \n");
  printf("       %-18s%-18s%-18s%-18s\n","ISYM","ISYM","ISYM","ISYM");
  for (i = 0 ; i < nrow-1 ; ++i) {
    printf("  ISYM");
    for (j = 0 ; j < 4 ; ++j) {
      for (k = 0; k < 3; ++k) {
        for (l = 0; l < 3; ++l) 
          tmp_symop[k][l] = - sp->invsymop[4*i+j].rot[l][k];
        tmp_symop[k][3] = 0.0;
        tmp_symop[3][k] = 0.0;
      }
      tmp_symop[3][3] = 1.0;
      mat4_to_recip_symop(rsymop,rsymop+rsymop_len,(const float (*)[4])tmp_symop);
      rsymop[12] = '\0';
      printf(" %3d  %-12s",2*(4*i+j)+2,rsymop);
    }
    printf("\n");
  }
  printf("  ISYM");
  for (j = 0 ; j < n_in_last_row ; ++j) {
    for (k = 0; k < 3; ++k) {
      for (l = 0; l < 3; ++l) 
        tmp_symop[k][l] = - sp->invsymop[4*i+j].rot[l][k];
      tmp_symop[k][3] = 0.0;
      tmp_symop[3][k] = 0.0;
    }
    tmp_symop[3][3] = 1.0;
    mat4_to_recip_symop(rsymop,rsymop+rsymop_len,(const float (*)[4])tmp_symop);
    rsymop[12] = '\0';
    printf(" %3d  %-12s",2*(4*(nrow-1)+j)+2,rsymop);
  }
  printf("\n");
}

int range_to_limits(const char *range, float limits[2])
{
  int i,in_value=1,neg=0,frac=0,equal=0;
  float value1,value2;
  float delta=0.00001;
  char ch;
  char buf[2];
  buf[1] = 0;

  for (i = 0 ; i < strlen(range) ; ++i) {
    ch = range[i];
    if (ch == '<') {
      if (in_value) {
	/* finishing lower value */
        limits[0] = value1;
        if (frac) limits[0] = value1/value2;
        if (neg) limits[0] = - limits[0];
        limits[0] += delta;
        neg = 0;
        frac = 0;
        in_value = 0;
      } else {
	/* starting upper value */

        in_value = 1;
      }
    } else if (ch == '-') {
      neg = 1;
    } else if (ch == '/') {
      frac = 1;
    } else if (ch == '=') {
      if (in_value) {
        equal = 1;
      } else {
        limits[0] -= 2.0*delta;        
      }
    } else if (ch == ';' || ch == ' ') {
      ;
    } else {
      if (in_value) {
	buf[0] = ch;
        if (frac) {
          value2 = (float) atoi(buf);
        } else {
          value1 = (float) atoi(buf);
        }
      }
    }
  } 
  /* finishing upper value */
  limits[1] = value1;
  if (frac) limits[1] = value1/value2;
  if (neg) limits[1] = - limits[1];
  limits[1] -= delta;
  if (equal) limits[1] += 2.0*delta;        

  return 0;
}

void set_fft_grid(CCP4SPG* sp, const int nxmin, const int nymin, const int nzmin, 
     const float sample, int *nx, int *ny, int *nz) 
{
  if (!sp) {  
    ccp4_signal(CSYM_ERRNO(CSYMERR_NullSpacegroup),"set_fft_grid",NULL); 
    return;
  }
  *nx = get_grid_sample(nxmin, sp->laue_sampling[0], sample);
  *ny = get_grid_sample(nymin, sp->laue_sampling[1], sample);
  *nz = get_grid_sample(nzmin, sp->laue_sampling[2], sample);
}

int all_factors_le_19(const int n)
{
  int i,ifact[8]={2,3,5,7,11,13,17,19};

  int nn = n;

  for (i = 0 ; i < 8 ; ++i) {
    while (nn % ifact[i] == 0) {
      /* factor found, divide & continue if required */
      nn = nn/ifact[i];
      /* success */
      if (nn == 1) 
	return 1;
    }
  }
  return 0;
}

int get_grid_sample(const int minsmp, const int nmul, const float sample)
{
  int n;
  float r1min=1.0, r1max=1.6, r2min=1.4, r2max=4.0;

  /*  check minsmp <= 0, if so set nsampl = nmul */
  if (minsmp <= 0) 
    return nmul;

  /* set search limits */
  if (sample >= 1.0) {
    r1max = sample;
    r2min = sample*0.95 < 1.0 ? 1.0 : sample*0.95;
  }

  /*  start with multiple of nmul */
  n = (int) rint((r1max*minsmp)/nmul)*nmul;

  while (n > (int) rint(r1min*minsmp)) {
    /* suitable sample interval found, accept it */
    if (all_factors_le_19(n)) 
      return n;
    /* decrement trial value & continue if still in range */
    n -= nmul;
  }

  /*  now try 2nd search if 1st unsuccesful */
  n = (int) rint((r2min*minsmp)/nmul)*nmul;

  while (n < (int) rint(r2max*minsmp)) {
    /* suitable sample interval found, accept it */
    if (all_factors_le_19(n)) 
      return n;
    /* increment trial value & continue if still in range */
    n += nmul;
  }

  /* failed */
  return -1;
}

int ccp4spg_load_laue(CCP4SPG *spacegroup, const int nlaue)
{
  int ierr = 1;

  if (!spacegroup) return ierr;

  if ( nlaue == 3 ) {
     spacegroup->asufn = &ASU_1b;
     spacegroup->nlaue = 3;
     strcpy(spacegroup->laue_name,"-1");
     spacegroup->laue_sampling[0] = 2;
     spacegroup->laue_sampling[1] = 2;
     spacegroup->laue_sampling[2] = 2;
     ierr = 0;
  }
  if ( nlaue == 4 ) {
     spacegroup->nlaue = 4;
     strcpy(spacegroup->laue_name,"2/m");
     spacegroup->laue_sampling[0] = 2;
     spacegroup->laue_sampling[1] = 4;
     spacegroup->laue_sampling[2] = 2;
     ierr = 0;
  }
  if ( nlaue == 5 ) {
     spacegroup->nlaue = 5;
     strcpy(spacegroup->laue_name,"2/m");
     spacegroup->laue_sampling[0] = 2;
     spacegroup->laue_sampling[1] = 8;
     spacegroup->laue_sampling[2] = 4;
     ierr = 0;
  }
  if ( nlaue == 6 ) {
     spacegroup->nlaue = 6;
     strcpy(spacegroup->laue_name,"mmm");
     spacegroup->laue_sampling[0] = 4;
     spacegroup->laue_sampling[1] = 4;
     spacegroup->laue_sampling[2] = 4;
     ierr = 0;
  }
  if ( nlaue == 7 ) {
     spacegroup->nlaue = 7;
     strcpy(spacegroup->laue_name,"4/m");
     spacegroup->laue_sampling[0] = 4;
     spacegroup->laue_sampling[1] = 4;
     spacegroup->laue_sampling[2] = 8;
     ierr = 0;
  }
  if ( nlaue == 8 ) {
     spacegroup->nlaue = 8;
     strcpy(spacegroup->laue_name,"4/mmm");
     spacegroup->laue_sampling[0] = 4;
     spacegroup->laue_sampling[1] = 4;
     spacegroup->laue_sampling[2] = 8;
     ierr = 0;
  }
  if ( nlaue == 9 ) {
     spacegroup->nlaue = 9;
     strcpy(spacegroup->laue_name,"-3");
     spacegroup->laue_sampling[0] = 6;
     spacegroup->laue_sampling[1] = 6;
     spacegroup->laue_sampling[2] = 6;
     ierr = 0;
  }
  if ( nlaue == 10 ) {
     spacegroup->nlaue = 10;
     strcpy(spacegroup->laue_name,"3bar1m");
     spacegroup->laue_sampling[0] = 6;
     spacegroup->laue_sampling[1] = 6;
     spacegroup->laue_sampling[2] = 6;
     ierr = 0;
  }
  if ( nlaue == 11 ) {
     spacegroup->nlaue = 11;
     strcpy(spacegroup->laue_name,"3barm");
     spacegroup->laue_sampling[0] = 6;
     spacegroup->laue_sampling[1] = 6;
     spacegroup->laue_sampling[2] = 6;
     ierr = 0;
  }
  if ( nlaue == 12 ) {
     spacegroup->nlaue = 12;
     strcpy(spacegroup->laue_name,"6/m");
     spacegroup->laue_sampling[0] = 6;
     spacegroup->laue_sampling[1] = 6;
     spacegroup->laue_sampling[2] = 12;
     ierr = 0;
  }
  if ( nlaue == 13 ) {
     spacegroup->nlaue = 13;
     strcpy(spacegroup->laue_name,"6/mmm");
     spacegroup->laue_sampling[0] = 6;
     spacegroup->laue_sampling[1] = 6;
     spacegroup->laue_sampling[2] = 12;
     ierr = 0;
  }
  if ( nlaue == 14 ) {
     spacegroup->nlaue = 14;
     strcpy(spacegroup->laue_name,"m3bar");
     spacegroup->laue_sampling[0] = 4;
     spacegroup->laue_sampling[1] = 4;
     spacegroup->laue_sampling[2] = 4;
     ierr = 0;
  }
  if ( nlaue == 15 ) {
     spacegroup->nlaue = 15;
     strcpy(spacegroup->laue_name,"m3barm");
     spacegroup->laue_sampling[0] = 8;
     spacegroup->laue_sampling[1] = 8;
     spacegroup->laue_sampling[2] = 8;
     ierr = 0;
  }
  return ierr;
}

int ccp4spg_check_symm_cell(int nsym, float rsym[][4][4], float cell[6]) {

  CCP4SPG *spacegroup;
  int i,k,l,status=1;
  ccp4_symop *op1;

  if (nsym <= 0) return 0;

  /* identify spacegroup from supplied symops */
  op1 = (ccp4_symop *) ccp4_utils_malloc(nsym*sizeof(ccp4_symop));
  for (i = 0; i < nsym; ++i) {
    for (k = 0; k < 3; ++k) {
      for (l = 0; l < 3; ++l) {
	op1[i].rot[k][l] = rsym[i][k][l];
      }
      op1[i].trn[k] = rsym[i][k][3];
    }
  }
  spacegroup = ccp4_spgrp_reverse_lookup(nsym,op1);

  /* test cell against symmetry on case-by-case basis */
  if (strstr(spacegroup->symbol_xHM,":R")) {
    status = ccp4uc_is_rhombohedral(cell,0.01F);
  } else if (strstr(spacegroup->symbol_xHM,":H")) {
    status = ccp4uc_is_hexagonal(cell,0.01F);
  } else if (spacegroup->spg_num >= 168 && spacegroup->spg_num <= 194) {
    status = ccp4uc_is_hexagonal(cell,0.01F);
  }

  free(op1);

  return status;
}
