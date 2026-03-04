/*
     cmtzlib.c: functions for handling MTZ files
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

/* DO NOT put Doxygen comments here - put in cmtzlib.h */
/* See cmtzlib.h for descriptions of functions */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cmtzlib.h"
#include "ccp4_types.h"
#include "ccp4_array.h"
#include "ccp4_parser.h"
#include "ccp4_vars.h"
#include "ccp4_errno.h"
#include "ccp4_unitcell.h"
/* "$Id$" */

/* stuff for error reporting */
#define CMTZ_ERRNO(n) (CCP4_ERR_MTZ | (n))

/* error defs */
#define  CMTZERR_Ok                  0
#define  CMTZERR_NoChannel           1
#define  CMTZERR_NoFile              2
#define  CMTZERR_NoLogicalName       3
#define  CMTZERR_CantOpenFile        4
#define  CMTZERR_NoHeader            5
#define  CMTZERR_ReadFail            6
#define  CMTZERR_WriteFail           7
#define  CMTZERR_ParamError          8
#define  CMTZERR_Cellerr             9
#define  CMTZERR_FileStamp           10
#define  CMTZERR_SymErr              11
#define  CMTZERR_AllocFail           12
#define  CMTZERR_MaxFile             13
#define  CMTZERR_ParserFail          14
#define  CMTZERR_NotMTZ              15
#define  CMTZERR_DatasetIncomplete   16
#define  CMTZERR_NoArch              17
#define  CMTZERR_NullDataset         18
#define  CMTZERR_BadVersion          19
#define  CMTZERR_SYMINFIncomplete    20
#define  CMTZERR_COLUMNIncomplete    21
#define  CMTZERR_BadBatchHeader    22
#define  CMTZERR_DifferentVersion  23
#define  CMTZERR_ColTypeMismatch   24
#define  CMTZERR_ColGroupError     25
#define  CMTZERR_ColSourceError    26


MTZ *MtzGet(const char *logname, int read_refs)

{ return MtzGetUserCellTolerance(logname, read_refs, 0.002);
}

MTZ *MtzGetUserCellTolerance(const char *logname, int read_refs, const double cell_tolerance)

{ MTZ *mtz;
  CCP4File *filein;
  int istat, newproj, cset_warn=0, length;
  MTZCOL *colin[MCOLUMNS], *newcol;
  char *filename;
  char crysin[MXTALS][65],projin[MXTALS][65],crystal[65],project[65];
  double cellin[MXTALS][6],cell[6];
  int jxtalin[MSETS];
  char mkey[4], keyarg[76], hdrrec[MTZRECORDLENGTH+1], label[31], type[3];
  int i, j, hdrst, ntotcol, nref, ntotset=0, nbat, nhist=0, icolin;
  int ixtal, jxtal, iset, iiset, icset, nxtal=0, nset[MXTALS]={0}, isym=0;
  int indhigh[3],indlow[3],isort[5],ind_xtal,ind_set,ind_col[3],debug=0;
  float min,max,totcell[6],minres,maxres;
  float *refldata;
  double coefhkl[6];
  int k; long xmllen;

  /* For cparser */
  CCP4PARSERARRAY *parser;
  CCP4PARSERTOKEN *token=NULL;
  char *key;
  int ntok,iprint=0;

  /* For batches */
  int ibat,nintegers,nreals;
  float buf[NBATCHWORDS];
  int *intbuf = (int *) buf;
  float *fltbuf = buf + NBATCHINTEGERS;
  MTZBAT *batch;

  /* known headers */
  char known_headers[][5] =
    { "PROJ","DATA","DCEL","DRES","DWAV","VERS","TITL","CELL",
      "SORT","SYMI","SYMM","COLU","VALM","RESO","COLS","COLG",
      "NCOL","NDIF","CRYS","MTZH","MTZB","BH" };
  int n_known_headers = sizeof(known_headers)/sizeof(known_headers[0]);

  if (debug) 
    printf(" Entering MtzGet \n");

  /* check input */
  if (!logname) return NULL;

  /* Open the mtz file: */
  if (getenv(logname) != NULL) {
    filename = strdup(getenv(logname));
  } else {
    filename = strdup(logname);
  }

  if (debug) 
    printf(" Opening file %s \n",filename);

  filein = ccp4_file_open(filename,O_RDONLY);
  if (! filein ) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_CantOpenFile),"MtzGet",NULL);
    free(filename);
    return NULL;
  }

  if (debug) 
    printf(" File opened successfully \n");

  /* specify location of stamp as 2*sizeof(float), where float is default mode */
  ccp4_file_setstamp(filein, 2);
  /* Read architecture */
  istat = ccp4_file_rarch (filein);
  if (!istat) {
   ccp4_signal(CCP4_ERRLEVEL(2) | CMTZ_ERRNO(CMTZERR_NoArch),
                        "MtzGet", NULL);
  }

  parser = ccp4_parse_start(20);
  if (parser == NULL) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_ParserFail),"MtzGet",NULL);
    free(filename);
    ccp4_file_close(filein);
    return NULL;
  }
  /* Set some convenient pointers to members of the parser array */
  key   = parser->keyword;
  token = parser->token;

  /* return to beginning of the file */
  ccp4_file_seek (filein, 0, SEEK_SET);
  /* set reading characters */
  ccp4_file_setmode(filein,0);
  istat = ccp4_file_readchar(filein, (uint8 *) hdrrec, 4);
  /* We don't test all reads, but this one should trap for e.g. truncated files */
  if (istat == EOF || hdrrec[0] == '\0') {
    ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_ReadFail),"MtzGet",NULL);
    ccp4_parse_end(parser);
    ccp4_file_close(filein);
    free(filename);
    return NULL;
  }
  hdrrec[4] = '\0';
  ntok = ccp4_parser(hdrrec, MTZRECORDLENGTH, parser, iprint);

  if (!ccp4_keymatch(key,"MTZ")) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_NotMTZ),"MtzGet",NULL);
    ccp4_parse_end(parser);
    ccp4_file_close(filein);
    free(filename);
    return(NULL);
  }

  if (debug) 
    printf(" MTZ file confirmed \n");

  /* set reading integers */
  ccp4_file_setmode(filein,6);
  istat = ccp4_file_read(filein, (uint8 *) &hdrst, 1);
  if (debug) printf(" hdrst read as %d \n",hdrst);

  /* 1st Pass: Read ntotcol, nref, nbat and dataset info.  
     nxtal and nset are used to assign memory for MTZ structure.
     Position at top of header */
  /* We don't test all seeks, but this one might trap duff files */
  if ( ccp4_file_seek(filein, hdrst-1, SEEK_SET) ) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_ReadFail),"MtzGet",NULL);
    ccp4_parse_end(parser);
    ccp4_file_close(filein);
    free(filename);
    return NULL;
  }

  /* set up base dataset in case it is in not in file */
  iiset = 0;
  nxtal = 1;
  jxtalin[0]=0;
  strcpy(projin[0],"HKL_base");
  strcpy(crysin[0],"HKL_base");
  nset[0] = 1;

  if (debug) 
    printf(" Start first pass \n");

  strcpy(project,"dummy");
  strcpy(crystal,"dummy");
  ccp4_file_setmode(filein,0);
  istat = ccp4_file_readchar(filein, (uint8 *) hdrrec, MTZRECORDLENGTH);
  if (debug) 
    printf(" Read first header record with istat = %d \n",istat);
  /* We don't test all reads, but this one should trap for e.g. truncated files */
  if (istat == EOF || hdrrec[0] == '\0') {
    ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_ReadFail),"MtzGet",NULL);
    ccp4_parse_end(parser);
    ccp4_file_close(filein);
    free(filename);
    return NULL;
  }

  hdrrec[MTZRECORDLENGTH] = '\0';
  ntok = ccp4_parser(hdrrec, MTZRECORDLENGTH, parser, iprint);
  while (!ccp4_keymatch(key,"END")) {

    /* read total number of columns, reflections and batches */
    if (ccp4_keymatch(key, "NCOL")) {
      ntotcol = (int) token[1].value;
      nref = (int) token[2].value;
      nbat = (int) token[3].value;
    }

    /* read total number of datasets over all projects/crystals */
    else if (ccp4_keymatch(key, "NDIF")) {
      ntotset = (int) token[1].value;
      if (debug) printf(" MtzGet: NDIF is %d\n",ntotset);
    }

    /* PROJECT line. Projects are not part of data structure, but
       may imply new crystal in hierarchy. */
    else if (ccp4_keymatch(key, "PROJ")) {
      ++iiset;
      if (iiset >= MSETS) {
        if (ccp4_liberr_verbosity(-1))
          printf("MtzGet: Maximum number of datasets exceeded! \n");
        ccp4_parse_end(parser);
        ccp4_file_close(filein);
        free(filename);
        return NULL;
      }
      strcpy(project,"dummy");
      if (ntok > 2) strcpy(project,token[2].fullstring);
      strcpy(crystal,project);
      jxtal = -1;
      for (ixtal = 0; ixtal < nxtal; ++ixtal) {
        if (strcmp(projin[ixtal],project) == 0) {
          jxtal = ixtal;
          jxtalin[iiset] = jxtal;
        }
      }
      /* New project implies new crystal */
      newproj=0;
      if (jxtal == -1) {
        ++nxtal;
        if (nxtal > MXTALS) {
          if (ccp4_liberr_verbosity(-1))
            printf("MtzGet: Maximum number of crystals exceeded! \n");
          ccp4_parse_end(parser);
          ccp4_file_close(filein);
          free(filename);
          return NULL;
        }
        jxtalin[iiset]=nxtal-1;
        strcpy(projin[nxtal-1],project);
        strcpy(crysin[nxtal-1],crystal);
        newproj=1;
      }
    }

    /* CRYSTAL line. This will be missing in old files! 
       If the line is present but incomplete we treat it as missing. */
    else if (ccp4_keymatch(key, "CRYS")) {
      if (ntok >= 3) {
       strcpy(crystal,token[2].fullstring);
       if (newproj == 1) {
        strcpy(crysin[nxtal-1],crystal);
       } else { 
        jxtal = -1;
        for (ixtal = 0; ixtal < nxtal; ++ixtal) {
          if (strcmp(crysin[ixtal],crystal) == 0) {
            jxtal = ixtal;
            jxtalin[iiset] = jxtal;
          }
        }
        if (jxtal == -1) {
          ++nxtal;
          if (nxtal > MXTALS) {
            if (ccp4_liberr_verbosity(-1))
              printf("MtzGet: Maximum number of crystals exceeded! \n");
            ccp4_parse_end(parser);
            ccp4_file_close(filein);
            free(filename);
            return NULL;
          }
          jxtalin[iiset]=nxtal-1;
          strcpy(projin[nxtal-1],project);
          strcpy(crysin[nxtal-1],crystal);
        }
       }
      }
    }

    /* DATASET line. This should present for every dataset so use to
       increment dataset count. However, if this is the base dataset,
       don't increment as we already have it. */
    else if (ccp4_keymatch(key, "DATA")) {
      if ( ntok <= 2 || (ntok > 2 && strcmp(token[2].fullstring,"HKL_base")) )
        ++nset[jxtalin[iiset]];
    }

    /* DCELL line. */
    else if (ccp4_keymatch(key, "DCEL")) {
      for (i = 0; i < 6; ++i) 
        cell[i] = token[i+2].value;
      /* If old crystal but cell dimensions differ, make new crystal.
         This is primarily for old files with no CRYSTAL cards. 
         This test doesn't apply to base dataset. 
         Chosen tolerance is arbitrary - there is no single correct value! */
      if (jxtal > 0 && iiset > 0 && 
	  ccp4uc_cells_differ(cellin[jxtal], cell, cell_tolerance)) {
        if (debug) {
          printf(" MtzGet: Old crystal %d but new cell dimensions. \n",jxtal);
          for (i = 0; i < 6; ++i) 
            printf(" %lf %lf \n",cellin[jxtal][i],cell[i]);
	}
        ++nxtal;
        if (nxtal > MXTALS) {
          if (ccp4_liberr_verbosity(-1))
            printf("MtzGet: Maximum number of crystals exceeded! \n");
          ccp4_parse_end(parser);
          ccp4_file_close(filein);
          free(filename);
          return NULL;
        }
        strcpy(projin[nxtal-1],project);
        strcpy(crysin[nxtal-1],crystal);
	/* Try to make crystal name unique */
        sprintf(crysin[nxtal-1]+strlen(crystal),"%d",nxtal);
	/* correct DATASET increment */
        --nset[jxtalin[iiset]];
        jxtalin[iiset]=nxtal-1;
        ++nset[jxtalin[iiset]];
      }
      for (i = 0; i < 6; ++i) 
        cellin[jxtalin[iiset]][i] = cell[i];
    }

    istat = ccp4_file_readchar(filein, (uint8 *) hdrrec, MTZRECORDLENGTH);
    hdrrec[MTZRECORDLENGTH] = '\0';
    ntok = ccp4_parser(hdrrec, MTZRECORDLENGTH, parser, iprint);
  }

  if (debug) {
    printf(" MtzGet: Found %d crystals \n",nxtal);
    for (i = 0; i < nxtal; ++i) 
      printf(" MtzGet: Crystal %d has %d datasets \n",i+1,nset[i]);
  }

  if (debug) 
    printf(" MtzGet: end of 1st pass \n");

  /* Allocate memory for input MTZ file */
  if (! (mtz = MtzMalloc(nxtal, nset))) {
    ccp4_parse_end(parser);
    ccp4_file_close(filein);
    free(filename);
    return NULL;
  }
  if (debug) 
    printf(" MtzGet: created mtz \n");
  mtz->filein = filein;
  mtz->nref_filein = nref;
  mtz->nref = nref;
  mtz->ncol_read = ntotcol;
  mtz->n_orig_bat = nbat;
  mtz->batch = NULL;
  mtz->refs_in_memory = read_refs;

  /* set up base dataset in case it is in not in file */
  nset[0] = 0;
  for (i = 1; i < nxtal; ++i) 
    nset[i] = -1;
  iiset = 0;
  strcpy(mtz->xtal[0]->pname,"HKL_base");
  strcpy(mtz->xtal[0]->xname,"HKL_base");
  mtz->xtal[0]->xtalid = 0;
  mtz->xtal[0]->cell[0] = 0.0;
  mtz->xtal[0]->set[0]->setid = 0;
  strcpy(mtz->xtal[0]->set[0]->dname,"HKL_base");
  mtz->xtal[0]->set[0]->wavelength = 0.0;

  if (debug) 
    printf(" MtzGet: starting 2nd pass \n");

  /* 2nd Pass: Copy dataset information to MTZ structure.
     Position at top of header */
  ccp4_file_setmode(filein,6);
  ccp4_file_seek(filein, hdrst-1, SEEK_SET);

  /* Read dataset information */
  ccp4_file_setmode(filein,0);
  istat = ccp4_file_readchar(filein, (uint8 *) hdrrec, MTZRECORDLENGTH);
  hdrrec[MTZRECORDLENGTH] = '\0';
  ntok = ccp4_parser(hdrrec, MTZRECORDLENGTH, parser, iprint);
  while (strncmp((strncpy(mkey,hdrrec,4)),"END",3) != 0) {

    if (strncmp (mkey, "PROJ",4) == 0) {
      ++iiset;
      strcpy(mtz->xtal[jxtalin[iiset]]->pname,projin[jxtalin[iiset]]);
      strcpy(mtz->xtal[jxtalin[iiset]]->xname,crysin[jxtalin[iiset]]);
      mtz->xtal[jxtalin[iiset]]->xtalid = jxtalin[iiset] + 1;
    }

    else if (strncmp (mkey, "DATA",4) == 0) {
      if ( ntok <= 2 || (ntok > 2 && strcmp(token[2].fullstring,"HKL_base")) ) {
        iset = (int) token[1].value;
        ++nset[jxtalin[iiset]];
        /* Test that dataset exists (i.e. pointer is non-NULL) */
        if (!mtz->xtal[jxtalin[iiset]]->set[nset[jxtalin[iiset]]]) {
          ccp4_signal(CCP4_ERRLEVEL(3) | 
                      CMTZ_ERRNO(CMTZERR_NullDataset),"MtzGet",NULL);
          ccp4_parse_end(parser);
          ccp4_file_close(filein);
          free(filename);
          return NULL;
	}
        mtz->xtal[jxtalin[iiset]]->set[nset[jxtalin[iiset]]]->setid = iset;
        strcpy(mtz->xtal[jxtalin[iiset]]->set[nset[jxtalin[iiset]]]->dname,"dummy");
        if (ntok > 2) strcpy(mtz->xtal[jxtalin[iiset]]->set[nset[jxtalin[iiset]]]->dname,
                              token[2].fullstring);
      }
    }

    else if (strncmp (mkey, "DCEL",4) == 0) {
      for (i = 0; i < 6; ++i) 
        mtz->xtal[jxtalin[iiset]]->cell[i] = (float) token[i+2].value;
    }

    /* this keyword not in use yet */
    else if (strncmp (mkey, "DRES",4) == 0) {
      for (i = 0; i < 3; ++i) {
        indhigh[i] = (int) token[i+2].value;
        indlow[i] = (int) token[i+5].value;
      }
      MtzHklcoeffs(mtz->xtal[jxtalin[iiset]]->cell, coefhkl);
      mtz->xtal[jxtalin[iiset]]->resmax = MtzInd2reso(indhigh, coefhkl);
      mtz->xtal[jxtalin[iiset]]->resmin = MtzInd2reso(indlow, coefhkl);
    }

    else if (strncmp (mkey, "DWAV",4) == 0) {
      mtz->xtal[jxtalin[iiset]]->set[nset[jxtalin[iiset]]]->wavelength = (float) token[2].value;
    }

    istat = ccp4_file_readchar(filein, (uint8 *) hdrrec, MTZRECORDLENGTH);
    if (istat == EOF) {
      /* Unexpected end-of-file */
      ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_ReadFail),"MtzGet",NULL);
      ccp4_parse_end(parser);
      ccp4_file_close(filein);
      free(filename);
      return NULL;
    }
    hdrrec[MTZRECORDLENGTH] = '\0';
    ntok = ccp4_parser(hdrrec, MTZRECORDLENGTH, parser, iprint);
  }

  if (debug) 
    printf(" MtzGet: end of 2nd pass \n");

  /* 3rd Pass: Position at top of header */
  ccp4_file_setmode(filein,6);
  ccp4_file_seek(filein, hdrst-1, SEEK_SET);

  icolin = -1;
  ccp4_file_setmode(filein,0);
  istat = ccp4_file_readchar(filein, (uint8 *) hdrrec, MTZRECORDLENGTH);
  hdrrec[MTZRECORDLENGTH] = '\0';
  ntok = ccp4_parser(hdrrec, MTZRECORDLENGTH, parser, iprint);
  while (strncmp((strncpy(mkey,hdrrec,4)),"END",3) != 0) {

    if (debug) 
      printf(" MtzGet: header line %s \n",hdrrec);

    if (strncmp (mkey, "VERS",4) == 0) {
      if (atoi(hdrrec+10) != MTZ_MAJOR_VERSN) {
         if (ccp4_liberr_verbosity(-1))
           printf("Input MTZ file has major version %d and minor version %d \n",
	       atoi(hdrrec+10),atoi(hdrrec+12));
         ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_BadVersion),"MtzGet",NULL);
         ccp4_parse_end(parser);
         ccp4_file_close(filein);
         free(filename);
         return(NULL);
         }  
      if (atoi(hdrrec+12) != MTZ_MINOR_VERSN) {
         if (ccp4_liberr_verbosity(-1))
           printf("Input MTZ file has major version %d and minor version %d \n",
	       atoi(hdrrec+10),atoi(hdrrec+12));
         ccp4_signal(CCP4_ERRLEVEL(2) | CMTZ_ERRNO(CMTZERR_DifferentVersion),"MtzGet",NULL);
         }  
       }
    else if (strncmp (mkey, "TITL",4) == 0) {
       strncpy(mtz->title,hdrrec+6,70); 
       length = 70;
       while ((--length >= 0) && (mtz->title[length] == ' '));
       mtz->title[length+1] = '\0';
       }

    else if (strncmp (mkey, "CELL",4) == 0) {
      for (i = 0; i < 6; ++i) 
        totcell[i] = (float) token[i+1].value;
      for (i = 0; i < mtz->nxtal; ++i) {
        if (mtz->xtal[i]->cell[0] < 0.01) {
          mtz->xtal[i]->cell[0] = totcell[0];
          mtz->xtal[i]->cell[1] = totcell[1];
          mtz->xtal[i]->cell[2] = totcell[2];
          mtz->xtal[i]->cell[3] = totcell[3];
          mtz->xtal[i]->cell[4] = totcell[4];
          mtz->xtal[i]->cell[5] = totcell[5];
        }
      } 
    }
    else if (strncmp (mkey, "SORT",4) == 0) {
      for (i = 0; i < 5; ++i) 
        isort[i] = (int) token[i+1].value;
    }
    else if (strncmp (mkey, "SYMI",4) == 0) {
      /* Check that there are enough tokens in the header record */
      if (ntok < 7) {
	ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_SYMINFIncomplete),
                        "MtzGet", NULL);
        ccp4_parse_end(parser);
        ccp4_file_close(filein);
        free(filename);
	return(NULL);
      }
      mtz->mtzsymm.nsym = (int) token[1].value;
      mtz->mtzsymm.nsymp = (int) token[2].value;
      mtz->mtzsymm.symtyp = token[3].fullstring[0];
      mtz->mtzsymm.spcgrp = (int) token[4].value;
      strcpy(mtz->mtzsymm.spcgrpname,token[5].fullstring);
      strcpy(mtz->mtzsymm.pgname,token[6].fullstring);
      if (ntok > 7) {
        mtz->mtzsymm.spg_confidence = token[7].fullstring[0];
      } else {
        mtz->mtzsymm.spg_confidence = 'X';
      }
    }
    else if (strncmp (mkey, "SYMM",4) == 0) {
      symop_to_mat4(hdrrec+4,hdrrec+MTZRECORDLENGTH,mtz->mtzsymm.sym[isym++][0]);
       }

    else if (strncmp (mkey, "COLU",4) == 0) {
      /* Check that there are enough tokens in the header record */
      if (ntok < 5) {
	ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_COLUMNIncomplete),
                        "MtzGet", NULL);
        ccp4_parse_end(parser);
        ccp4_file_close(filein);
        free(filename);
	return(NULL);
      }
      ++icolin;
      if (icolin >= MCOLUMNS) {
        if (ccp4_liberr_verbosity(-1))
          printf("MtzGet: Maximum number of columns exceeded! \n");
        ccp4_parse_end(parser);
        ccp4_file_close(filein);
        free(filename);
        return NULL;
      }
      strcpy(label,token[1].fullstring);
      strcpy(type,token[2].fullstring);
      min = (float) token[3].value;
      max = (float) token[4].value;
      /* Dataset id for this column
	 Very old MTZ files may not have this value */
      if (ntok < 6) {
        if (!cset_warn) {
          if (ccp4_liberr_verbosity(-1)) {
            printf("\nWARNING: Dataset id missing from COLUMN records in MTZ header. \n");
            printf("WARNING: Making default dataset assignments. \n");
	  }
          ccp4_signal(CCP4_ERRLEVEL(2) | CMTZ_ERRNO(CMTZERR_DatasetIncomplete),
                        "MtzGet", NULL);
          cset_warn = 1;
	}
	icset = 0;
      } else {
	icset = (int) token[5].value;
      }
      /* Special trap for M/ISYM */
      if (type[0] == 'Y' && strncmp (label,"M/ISYM",6) == 0)
        strcpy(label,"M_ISYM");
      /* Find dataset corresponding to this column */
      ixtal = 0; iset = 0;
      for (i = 0; i < mtz->nxtal; ++i) {
       for (j = 0; j < mtz->xtal[i]->nset; ++j) {
        if (mtz->xtal[i]->set[j]->setid == icset) {
          ixtal = i;
          iset = j;
          break;
        }
       }
      }

      /* Create column. */
      newcol = MtzAddColumn(mtz, mtz->xtal[ixtal]->set[iset], label, type);
      newcol->source = icolin + 1;
      newcol->min = min;
      newcol->max = max;
      colin[icolin] = newcol;
    }

    else if (strncmp (mkey, "VALM",4) == 0) {
      strcpy(keyarg,token[1].fullstring);
      if (strncmp (keyarg, "NAN",3) == 0) {
        sprintf(mtz->mnf.amnf,"NAN");
      } else {
	mtz->mnf.fmnf = (float) token[1].value;
      }
    }

    else if (strncmp (mkey, "RESO",4) == 0) {
      minres = (float) token[1].value;
      maxres = (float) token[2].value;
      for (i = 0; i < mtz->nxtal; ++i) {
        if (mtz->xtal[i]->resmax == 0.0)
          mtz->xtal[i]->resmax = maxres;
        if (mtz->xtal[i]->resmin == 100.0)
          mtz->xtal[i]->resmin = minres;
      }
    }

    istat = ccp4_file_readchar(filein, (uint8 *) hdrrec, MTZRECORDLENGTH);
    hdrrec[MTZRECORDLENGTH] = '\0';
    ntok = ccp4_parser(hdrrec, MTZRECORDLENGTH, parser, iprint);
  }

  /* 4th Pass: Column group and source extensions and unknown keywords */
  /* 4th Pass: Position at top of header */
  ccp4_file_setmode(filein,6);
  ccp4_file_seek(filein, hdrst-1, SEEK_SET);
  ccp4_file_setmode(filein,0);
  istat = ccp4_file_readchar(filein, (uint8 *) hdrrec, MTZRECORDLENGTH);
  hdrrec[MTZRECORDLENGTH] = '\0';
  ntok = ccp4_parser(hdrrec, MTZRECORDLENGTH, parser, iprint);
  while (strncmp((strncpy(mkey,hdrrec,4)),"END",3) != 0) {
    if (strncmp (mkey, "COLS",4) == 0 ) {
      strcpy(label,token[1].fullstring);
      /* Special trap for M/ISYM */
      if (strncmp (label,"M/ISYM",6) == 0)
        strcpy(label,"M_ISYM");
      icset = (int) token[3].value;
      newcol = NULL;
      for (i = 0; i < mtz->nxtal; ++i) {
	for (j = 0; j < mtz->xtal[i]->nset; ++j) {
	  if (mtz->xtal[i]->set[j]->setid == icset) {
	    for ( k = 0; k < mtz->xtal[i]->set[j]->ncol; k++ ) {
	      if (strcmp(mtz->xtal[i]->set[j]->col[k]->label,label) == 0) {
		newcol = mtz->xtal[i]->set[j]->col[k];
		break;
	      }
	    }
	  }
	}
      }
      if ( newcol == NULL ) {
 	ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_ColSourceError),
		    "MtzGet", NULL);
	ccp4_parse_end(parser);
	ccp4_file_close(filein);
	free(filename);
	return(NULL);
      }
      strncpy( newcol->colsource, token[2].fullstring, 36 );
      newcol->colsource[36] = '\0';
    } else if (strncmp (mkey, "COLG",4) == 0 ) {
      strcpy(label,token[1].fullstring);
      /* Special trap for M/ISYM */
      if (strncmp (label,"M/ISYM",6) == 0)
        strcpy(label,"M_ISYM");
      icset = (int) token[5].value;
      newcol = NULL;
      for (i = 0; i < mtz->nxtal; ++i) {
	for (j = 0; j < mtz->xtal[i]->nset; ++j) {
	  if (mtz->xtal[i]->set[j]->setid == icset) {
	    for ( k = 0; k < mtz->xtal[i]->set[j]->ncol; k++ ) {
	      if (strcmp(mtz->xtal[i]->set[j]->col[k]->label,label) == 0) {
		newcol = mtz->xtal[i]->set[j]->col[k];
		break;
	      }
	    }
	  }
	}
      }
      if ( newcol == NULL ) {
 	ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_ColGroupError),
		    "MtzGet", NULL);
	ccp4_parse_end(parser);
	ccp4_file_close(filein);
	free(filename);
	return(NULL);
      }
      strncpy( newcol->grpname, token[2].fullstring, 30 );
      newcol->grpname[30] = '\0';
      strncpy( newcol->grptype, token[3].fullstring, 4 );
      newcol->grptype[4] = '\0';
      newcol->grpposn = (int) token[4].value;
    }
    istat = ccp4_file_readchar(filein, (uint8 *) hdrrec, MTZRECORDLENGTH);
    hdrrec[MTZRECORDLENGTH] = '\0';
    ntok = ccp4_parser(hdrrec, MTZRECORDLENGTH, parser, iprint);
  }

  /* 5th Pass: Deal with unknown headers */
  /* 5th Pass: Position at top of header */
  ccp4_file_setmode(filein,6);
  ccp4_file_seek(filein, hdrst-1, SEEK_SET);
  ccp4_file_setmode(filein,0);
  istat = ccp4_file_readchar(filein, (uint8 *) hdrrec, MTZRECORDLENGTH);
  hdrrec[MTZRECORDLENGTH] = '\0';
  ntok = ccp4_parser(hdrrec, MTZRECORDLENGTH, parser, iprint);
  while (strncmp((strncpy(mkey,hdrrec,4)),"END",3) != 0) {
    for ( i = 0; i < n_known_headers; ++i )
      if (strncmp (mkey,known_headers[i],4) == 0 )
	break;
    if ( i == n_known_headers ) {
      mtz->unknown_headers = ccp4_utils_realloc( mtz->unknown_headers, mtz->n_unknown_headers*MTZRECORDLENGTH+MTZRECORDLENGTH );  // if null, malloc
      memcpy( mtz->unknown_headers+mtz->n_unknown_headers*MTZRECORDLENGTH, hdrrec, MTZRECORDLENGTH );
      mtz->n_unknown_headers++;
    }
    istat = ccp4_file_readchar(filein, (uint8 *) hdrrec, MTZRECORDLENGTH);
    hdrrec[MTZRECORDLENGTH] = '\0';
    ntok = ccp4_parser(hdrrec, MTZRECORDLENGTH, parser, iprint);
  }

  /* copy sort order */
  for (i = 0; i < 5; ++i) {
    if (isort[i] > 0) 
      mtz->order[i] = colin[isort[i]-1];
  }

  if (debug) 
    printf(" MtzGet: end of 3rd pass \n");

  /* Now read history if any */
  istat = ccp4_file_readchar(filein, (uint8 *) hdrrec, MTZRECORDLENGTH);
  hdrrec[MTZRECORDLENGTH] = '\0';
  ntok = ccp4_parser(hdrrec, MTZRECORDLENGTH, parser, iprint);
  while (!ccp4_keymatch(key,"MTZE")) {

    if (ccp4_keymatch(key, "MTZH")) {
      nhist = (int) token[1].value;
      /* allocate memory for nhist lines */
      mtz->hist = MtzCallocHist(nhist);
      mtz->histlines = nhist;

      for (i = 0; i < nhist; ++i) {
        istat = ccp4_file_readchar(filein, (uint8 *) hdrrec, MTZRECORDLENGTH);
        strncpy(mtz->hist + MTZRECORDLENGTH*i,hdrrec,MTZRECORDLENGTH);
      }

    } else if (ccp4_keymatch(key, "MTZB")) {
      for (ibat = 0; ibat < nbat; ++ibat) {

        istat = ccp4_file_readchar(filein, (uint8 *) hdrrec, MTZRECORDLENGTH);
        hdrrec[MTZRECORDLENGTH] = '\0';
        ntok = ccp4_parser(hdrrec, MTZRECORDLENGTH, parser, iprint);
        if (!ccp4_keymatch(key, "BH")) {
          ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_BadBatchHeader),
                        "MtzGet", NULL);
          ccp4_parse_end(parser);
          ccp4_file_close(filein);
          free(filename);
          return(NULL);
        }

	/* allocate memory for this batch */
        if (ibat == 0) {
          mtz->batch = MtzMallocBatch();
          batch = mtz->batch;
        } else {
          batch->next = MtzMallocBatch();
          batch = batch->next;
        }
        batch->next = NULL;
        batch->num = (int) token[1].value;
        /* nwords = (int) token[2].value; */
        nintegers = (int) token[3].value;
        nreals = (int) token[4].value;
	/* read batch title */
        istat = ccp4_file_readchar(filein, (uint8 *) hdrrec, MTZRECORDLENGTH);
        strncpy(batch->title,hdrrec+6,70); 
        batch->title[70]='\0';

        ccp4_file_setmode(filein,6);
        istat = ccp4_file_read(filein, (uint8 *) intbuf, nintegers);
        ccp4_file_setmode(filein,2);
        istat = ccp4_file_read(filein, (uint8 *) fltbuf, nreals);

        MtzArrayToBatch(intbuf, fltbuf, batch);

        ccp4_file_setmode(filein,0);
	istat = ccp4_file_readchar(filein, (uint8 *) hdrrec, MTZRECORDLENGTH); 
        hdrrec[MTZRECORDLENGTH] = '\0';
        ntok = ccp4_parser(hdrrec, MTZRECORDLENGTH, parser, iprint);
        if (ntok == 4) {
          strcpy(batch->gonlab[0],token[1].fullstring); 
          strcpy(batch->gonlab[1],token[2].fullstring); 
          strcpy(batch->gonlab[2],token[3].fullstring); 
          batch->gonlab[0][8] = batch->gonlab[1][8] = batch->gonlab[2][8] = '\0';
        } else if (ntok == 2) {
          strcpy(batch->gonlab[0],token[1].fullstring);
          batch->gonlab[0][8] = '\0';
          batch->gonlab[1][0] = batch->gonlab[2][0] = '\0';
        } else {
          batch->gonlab[0][0] = batch->gonlab[1][0] = batch->gonlab[2][0] = '\0';
	}
      }
    }

    istat = ccp4_file_readchar(filein, (uint8 *) hdrrec, MTZRECORDLENGTH);
    hdrrec[MTZRECORDLENGTH] = '\0';
    ntok = ccp4_parser(hdrrec, MTZRECORDLENGTH, parser, iprint);
  }

  /* Finished with the parser array */
  ccp4_parse_end(parser);

  if (debug) 
    printf(" MtzGet: end of batch pass \n");

  /* Read XML datablock */
  xmllen = ccp4_file_length(filein) - ccp4_file_tell(filein);
  if ( xmllen > 0 ) {
    mtz->xml = (char *)ccp4_utils_malloc( xmllen+1 );
    if ( mtz->xml != NULL ) {
      istat = ccp4_file_readchar(filein, (uint8 *) mtz->xml, xmllen);
      mtz->xml[xmllen] = '\0';
    }
  }

  /* Position at start of reflections */
  ccp4_file_setmode(filein,6);
  ccp4_file_seek(filein, SIZE1, SEEK_SET);

  if (read_refs) {

    refldata = (float *) ccp4_utils_malloc(ntotcol*sizeof(float));
    /* Read all reflections into memory - make this optional? */
    for (i = 0; i < mtz->nref_filein; ++i) {
      MtzRrefl(filein, ntotcol, refldata);
      for (j = 0; j < ntotcol; ++j)
        colin[j]->ref[i] = refldata[j];
    }
    free(refldata);

    /* Recalculate resolution limits */

    /* Find dataset of indices */
    MtzFindInd(mtz,&ind_xtal,&ind_set,ind_col);

    for (i = 0; i < mtz->nxtal; ++i) {
      MtzHklcoeffs(mtz->xtal[i]->cell, coefhkl);
      for (j = 0; j < mtz->nref; ++j) {
        indhigh[0] = (int) mtz->xtal[ind_xtal]->set[ind_set]->col[ind_col[0]]->ref[j];
        indhigh[1] = (int) mtz->xtal[ind_xtal]->set[ind_set]->col[ind_col[1]]->ref[j];
        indhigh[2] = (int) mtz->xtal[ind_xtal]->set[ind_set]->col[ind_col[2]]->ref[j];
        maxres = MtzInd2reso(indhigh, coefhkl);
        if (maxres > mtz->xtal[i]->resmax) mtz->xtal[i]->resmax = maxres;
        if (maxres < mtz->xtal[i]->resmin) mtz->xtal[i]->resmin = maxres;
      }
    }

    /* And close the mtz file: */
    ccp4_file_close(filein);
    mtz->filein = NULL;
  }

  free(filename); 

  return(mtz);}

int MtzArrayToBatch(const int *intbuf, const float *fltbuf, MTZBAT *batch)

{  int i;

  batch->iortyp = intbuf[3];
  for (i = 0; i < 6; ++i)
    batch->lbcell[i] = intbuf[4 + i];
  batch->misflg = intbuf[10];
  batch->jumpax = intbuf[11];
  batch->ncryst = intbuf[12];
  batch->lcrflg = intbuf[13];
  batch->ldtype = intbuf[14];
  batch->jsaxs = intbuf[15];
  batch->nbscal = intbuf[16];
  batch->ngonax = intbuf[17];
  batch->lbmflg = intbuf[18];
  batch->ndet = intbuf[19];
  batch->nbsetid = intbuf[20];

  for (i = 0; i < 6; ++i)
    batch->cell[i] = fltbuf[i];
  for (i = 0; i < 9; ++i)
    batch->umat[i] = fltbuf[6 + i];
  for (i = 0; i < 3; ++i) 
    batch->phixyz[0][i] = fltbuf[15 + i];
  for (i = 0; i < 3; ++i) 
    batch->phixyz[1][i] = fltbuf[18 + i];
  for (i = 0; i < 12; ++i) 
    batch->crydat[i] = fltbuf[21 + i];
  for (i = 0; i < 3; ++i) 
    batch->datum[i] = fltbuf[33 + i];
  batch->phistt = fltbuf[36];
  batch->phiend = fltbuf[37];
  for (i = 0; i < 3; ++i) 
    batch->scanax[i] = fltbuf[38 + i];
  batch->time1 = fltbuf[41];
  batch->time2 = fltbuf[42];
  batch->bscale = fltbuf[43];
  batch->bbfac = fltbuf[44];
  batch->sdbscale = fltbuf[45];
  batch->sdbfac = fltbuf[46];
  batch->phirange = fltbuf[47];
  for (i = 0; i < 3; ++i)
    batch->e1[i] = fltbuf[59 + i];
  for (i = 0; i < 3; ++i)
    batch->e2[i] = fltbuf[62 + i];
  for (i = 0; i < 3; ++i)
    batch->e3[i] = fltbuf[65 + i];
  for (i = 0; i < 3; ++i)
    batch->source[i] = fltbuf[80 + i];
  for (i = 0; i < 3; ++i)
    batch->so[i] = fltbuf[83 + i];
  batch->alambd = fltbuf[86];
  batch->delamb = fltbuf[87];
  batch->delcor = fltbuf[88];
  batch->divhd = fltbuf[89];
  batch->divvd = fltbuf[90];
  for (i = 0; i < 2; ++i)
  { batch->dx[i] = fltbuf[111 + (i * 6)];
    batch->theta[i] = fltbuf[112 + (i * 6)];
    batch->detlm[i][0][0] = fltbuf[113 + (i * 6)];
    batch->detlm[i][0][1] = fltbuf[114 + (i * 6)];
    batch->detlm[i][1][0] = fltbuf[115 + (i * 6)];
    batch->detlm[i][1][1] = fltbuf[116 + (i * 6)];}

  return 1;
}

int MtzRrefl(CCP4File *filein, int ncol, float *refldata) {

  int istat;

  ccp4_file_setmode(filein,2);
  istat = ccp4_file_read(filein, (uint8 *) refldata, ncol);

  /* This will return EOF if end-of-file is reached. But by then
     you will have read the MTZ file header, so not so useful. */
  return istat;
}

int MtzFindInd(const MTZ *mtz, int *ind_xtal, int *ind_set, int ind_col[3]) {

  int i,j,k;

  /* default to first 3 columns of 1st datset */
  *ind_xtal = 0;
  *ind_set = 0;
  ind_col[0] = 0;
  ind_col[1] = 1;
  ind_col[2] = 2;  

  for (i = 0; i < mtz->nxtal; ++i) 
   for (j = 0; j < mtz->xtal[i]->nset; ++j) 
    for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) {
     if (mtz->xtal[i]->set[j]->col[k]->label[0] == 'H' &&
         mtz->xtal[i]->set[j]->col[k]->type[0] == 'H') {
       *ind_xtal = i;
       *ind_set = j;
       ind_col[0] = k;
     }
     if (mtz->xtal[i]->set[j]->col[k]->label[0] == 'K' &&
         mtz->xtal[i]->set[j]->col[k]->type[0] == 'H') 
       ind_col[1] = k;
     if (mtz->xtal[i]->set[j]->col[k]->label[0] == 'L' &&
         mtz->xtal[i]->set[j]->col[k]->type[0] == 'H') 
       ind_col[2] = k;
    }

  return 1;
}

float MtzInd2reso(const int in[3], const double coefhkl[6]) {

  int ih,ik,il;
  float reso;

  ih = in[0];
  ik = in[1];
  il = in[2];

  reso = (float) 4.0 * (ih*ih*coefhkl[0] + ih*ik*coefhkl[1] + 
             ih*il*coefhkl[2] + ik*ik*coefhkl[3] + 
             ik*il*coefhkl[4] + il*il*coefhkl[5]);

  return reso;

}

int MtzHklcoeffs(const float cell[6], double coefhkl[6]) {

  /* generate coefhkl coefficients from given cell parameters */

  int i;
  double alpha,beta,gamma,degtorad,denom;
  double ax,bx,by,cx,cy,cz;
  double axst,ayst,azst,byst,bzst,czst;

  /* sanity clause (but there ain't no sanity clause!) */
  for (i = 0; i < 6; ++i) 
    coefhkl[i] = 0.0;
  for (i = 0; i < 6; ++i)
    if (cell[i] < 0.001) {
      ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_Cellerr),"MtzHklcoeffs",NULL);
      return 0;
    }

  degtorad = atan(1.0)/45.0;

  alpha = degtorad*cell[3];
  beta = degtorad*cell[4];
  gamma = degtorad*cell[5];

  /* orthogonal frame for calculation 
   a along x, b in x-y plane */

  ax = cell[0];
  bx = cell[1] * cos(gamma);
  by = cell[1] * sin(gamma);
  cx = cell[2] * cos(beta);
  cy = (cell[1]*cell[2]*cos(alpha) - bx*cx)/by;
  cz = sqrt(cell[2]*cell[2] - cx*cx - cy*cy);

  /* find reciprocal vectors in orthogonal frame */

  denom = ax*by*cz;
  axst = 1.0/ax;
  ayst = -bx*cz/denom;
  azst = (bx*cy - by*cx)/denom;
  byst = 1.0/by;
  bzst = -ax*cy/denom;
  czst = 1.0/cz;

  coefhkl[0] = 0.25*(axst*axst + ayst*ayst + azst*azst);
  coefhkl[1] = 0.5*(ayst*byst + azst*bzst);
  coefhkl[2] = 0.5*(azst*czst);
  coefhkl[3] = 0.25*(byst*byst + bzst*bzst);
  coefhkl[4] = 0.5*(bzst*czst);
  coefhkl[5] = 0.25*(czst*czst);

  return 1;
}

int ccp4_lrtitl(const MTZ *mtz, char *title) {

  int length;

  length = (int) strlen(strcpy(title, mtz->title));
  if (length > 0) {
    while ((--length >= 0) && (title[length] == ' '));
    ++length;
  }
  return(length);
}

int ccp4_lrhist(const MTZ *mtz, char history[][MTZRECORDLENGTH], int nlines) {

  int i,nhist;

  if (nlines < mtz->histlines)
    nhist = nlines;
  else
    nhist = mtz->histlines;

  for (i = 0; i < nhist; ++i) {
    strncpy(history[i],mtz->hist + MTZRECORDLENGTH*i,MTZRECORDLENGTH);
  }

  return nhist;
}

int ccp4_lrsort(const MTZ *mtz, int isort[5]) {

  int i,j,k,l,icol;

  icol = 0;
  for (i = 0; i < 5; ++i) 
    isort[i] = 0;
  /* Loop over crystals */
  for (i = 0; i < mtz->nxtal; ++i) {
    /* Loop over datasets for each crystal */
    for (j = 0; j < mtz->xtal[i]->nset; ++j) {
      /* Loop over columns for each dataset */
      for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) {
        ++icol;
        for (l = 0; l < 5; ++l) {
          if (mtz->order[l] == mtz->xtal[i]->set[j]->col[k])
            isort[l] = icol;
        }
      }
    }
  }

  return 1;
}

int ccp4_lrbats(const MTZ *mtz, int *nbatx, int batchx[]) {

  int i=0;
  MTZBAT *batch;

  *nbatx = mtz->n_orig_bat;
  batch = mtz->batch;
  while (batch != NULL) {
    batchx[i++] = batch->num;
    batch = batch->next;
  }

  return i;
}

void MtzDebugHierarchy(const MTZ *mtz) {

  int i,j,k;

  if (mtz->filein)
    printf("MtzDebugHierarchy: input file = %s \n",mtz->filein->name);
  if (mtz->fileout)
    printf("MtzDebugHierarchy: output file = %s \n",mtz->fileout->name);

  printf("MtzDebugHierarchy: nxtal = %d \n",mtz->nxtal);
  for (i = 0; i < mtz->nxtal; ++i) {
   printf("MtzDebugHierarchy: xtal = %s, cell = %f %f %f %f %f %f \n",
          mtz->xtal[i]->xname,
	  mtz->xtal[i]->cell[0],mtz->xtal[i]->cell[1],mtz->xtal[i]->cell[2],
	  mtz->xtal[i]->cell[3],mtz->xtal[i]->cell[4],mtz->xtal[i]->cell[5]);
   printf("MtzDebugHierarchy: xtal = %s, nset = %d \n",mtz->xtal[i]->xname,
              mtz->xtal[i]->nset);
   for (j = 0; j < mtz->xtal[i]->nset; ++j) {
    printf("MtzDebugHierarchy: xtal = %s, set = %s, setid = %d, ncol = %d \n",
              mtz->xtal[i]->xname,mtz->xtal[i]->set[j]->dname,
              mtz->xtal[i]->set[j]->setid,mtz->xtal[i]->set[j]->ncol);
     for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) {
      printf("MtzDebugHierarchy: col = %s (in: %d) (out: %d) \n",
              mtz->xtal[i]->set[j]->col[k]->label,
              mtz->xtal[i]->set[j]->col[k]->source,
              mtz->xtal[i]->set[j]->col[k]->active);
     }
   }
  }

}

/* List of column information: label, type, dataset.
   Returns number of columns in current structure. */
int MtzListColumn(const MTZ *mtz, char clabs[][31], char ctyps[][3], int csetid[]) {

 int i,j,k,icol=0;

 /* Loop over crystals */
   for (i = 0; i < mtz->nxtal; ++i) {
 /* Loop over datasets for each crystal */
    for (j = 0; j < mtz->xtal[i]->nset; ++j) {
 /* Loop over columns for each dataset */
     for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) {
       if (strcmp(mtz->xtal[i]->set[j]->col[k]->type,"Y") == 0 && 
           strcmp(mtz->xtal[i]->set[j]->col[k]->label,"M_ISYM") == 0) {
         strcpy(clabs[icol],"M/ISYM");
       } else {
         strcpy(clabs[icol],mtz->xtal[i]->set[j]->col[k]->label);
       }
       strcpy(ctyps[icol],mtz->xtal[i]->set[j]->col[k]->type);
       csetid[icol] = mtz->xtal[i]->set[j]->setid;
       ++icol;
     }
    }
   }
   return icol;
}

/* List of column information from input file: label, type, dataset.
   Returns number of columns in input file. */
int MtzListInputColumn(const MTZ *mtz, char clabs[][31], char ctyps[][3], int csetid[]) {

 int i,j,k,colin,icol=0;

 /* Loop over crystals */
   for (i = 0; i < mtz->nxtal; ++i) {
 /* Loop over datasets for each crystal */
    for (j = 0; j < mtz->xtal[i]->nset; ++j) {
 /* Loop over columns for each dataset */
     for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) {
      if ((colin = mtz->xtal[i]->set[j]->col[k]->source) != 0) {
       if (strcmp(mtz->xtal[i]->set[j]->col[k]->type,"Y") == 0 && 
           strcmp(mtz->xtal[i]->set[j]->col[k]->label,"M_ISYM") == 0) {
         strcpy(clabs[colin - 1],"M/ISYM");
       } else {
         strcpy(clabs[colin - 1],mtz->xtal[i]->set[j]->col[k]->label);
       }
       strcpy(ctyps[colin - 1],mtz->xtal[i]->set[j]->col[k]->type);
       csetid[colin - 1] = mtz->xtal[i]->set[j]->setid;
       ++icol;
      }
     }
    }
   }
   return icol;
}

int ccp4_lrcell(const MTZXTAL *xtl, float cell[]) {

  int i;

  for (i = 0; i < 6; ++i) {
    cell[i] = xtl->cell[i];
  }
 
  return 1;
}

int MtzResLimits(const MTZ *mtz, float *minres, float *maxres) {

  int i;

  *maxres = 0.0;
  *minres = 100.0;
 /* Loop over crystals */
  for (i = 0; i < mtz->nxtal; ++i) {
    if (mtz->xtal[i]->resmax > *maxres) *maxres = mtz->xtal[i]->resmax;
    if (mtz->xtal[i]->resmin < *minres) *minres = mtz->xtal[i]->resmin;
  }

  return 1;
}

int ccp4_lrsymi(const MTZ *mtz, int *nsympx, char *ltypex, int *nspgrx, 
       char *spgrnx, char *pgnamx) {
  char spgconf_temp[2];

  return ccp4_lrsymi_c(mtz,nsympx,ltypex,nspgrx,spgrnx,pgnamx,spgconf_temp);
}

int ccp4_lrsymi_c(const MTZ *mtz, int *nsympx, char *ltypex, int *nspgrx, 
       char *spgrnx, char *pgnamx, char *spgconf) {

  *nsympx = mtz->mtzsymm.nsymp;
  *nspgrx = mtz->mtzsymm.spcgrp;
  ltypex[0] = mtz->mtzsymm.symtyp;
  ltypex[1] = '\0';
  strcpy(spgrnx,mtz->mtzsymm.spcgrpname);
  strcpy(pgnamx,mtz->mtzsymm.pgname);
  spgconf[0] = mtz->mtzsymm.spg_confidence;
  spgconf[1] = '\0';

  return *nspgrx;
}

int MtzSpacegroupNumber(const MTZ *mtz)
/* get the spacegroup number (likely CCP4 convention) */
{
  if (!mtz) return 0;
  return mtz->mtzsymm.spcgrp;
}

int ccp4_lrsymm(const MTZ *mtz, int *nsymx, float rsymx[192][4][4]) {

  int i,j,k;

  *nsymx = mtz->mtzsymm.nsym;
  for (i = 0; i < *nsymx; ++i) {
    for (j = 0; j < 4; ++j) {
      for (k = 0; k < 4; ++k) {
        rsymx[i][j][k] = mtz->mtzsymm.sym[i][j][k];
      }
    }
  }

  return *nsymx;
}

int MtzParseLabin(char *labin_line, const char prog_labels[][31], 
            const int nlprgi, char user_labels[][2][31]) 
 
{ int i,j,imatch,nlabels=0,err=0;
  char label1[31],label2[31];

  /* For cparser */
  CCP4PARSERARRAY *parser=NULL;
  CCP4PARSERTOKEN *token=NULL;
  char *key;
  int ntok,iprint=0;

  parser = ccp4_parse_start(strlen(labin_line));
  if (parser == NULL) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_ParserFail),"MtzParseLabin",NULL);
    return -1;
  }
  /* Set some convenient pointers to members of the parser array */
  key   = parser->keyword;
  token = parser->token;

  ntok = ccp4_parser(labin_line, strlen(labin_line), parser, iprint);

  if (ccp4_keymatch(key,"LABI")) {
    if (iprint) printf("Interpreting LABIN line.\n");
  } else if (ccp4_keymatch(key,"LABO")) {
    if (iprint) printf("Interpreting LABOUT line.\n");
  } else if (ccp4_keymatch(key,"COMP")) {
    if (iprint) printf("Interpreting COMPLETE line (freerflag).\n");
  } else {
    printf("Warning in MtzParseLabin: Input is not LABIN or LABOUT line !!\n");
  }

  /* initialise user labels */
  for (j = 0; j < nlprgi; ++j) {
    strcpy(user_labels[j][0],"");
    strcpy(user_labels[j][1],"");
  }

  for (i = 1; i < ntok; i += 2) {
    strcpy(label1,token[i].fullstring);

    if (strlen(label1)>30) {
      printf("MtzParseLabin: labels cannot be longer than 30 characters: \"%s\"\n",label1);
      err++;
      break;
    }

    /* Trap against trying to access tokens that don't exist */
    if (i+1 < ntok) {
      strcpy(label2,token[i+1].fullstring);

      if (strlen(label2)>30) {
        printf("MtzParseLabin: labels cannot be longer than 30 characters: \"%s\"\n",label2);
        err++;
        break;
      }

      /* check first label against program labels */
      imatch = 0;
      for (j = 0; j < nlprgi; ++j) {
	if (strcmp(label1,prog_labels[j]) == 0) {
	  strcpy(user_labels[j][0],label1);
	  strcpy(user_labels[j][1],label2);
	  imatch = 1;
	  ++nlabels;
	  break;
	}
      }
 
      if (imatch == 0) {
	/* check second label against program labels */
	for (j = 0; j < nlprgi; ++j) {
	  if (strcmp(label2,prog_labels[j]) == 0) {
	    strcpy(user_labels[j][0],label2);
	    strcpy(user_labels[j][1],label1);
	    imatch = 1;
	    ++nlabels;
	    break;
	  }
	}
      }
    } else {
      printf("MtzParseLabin: run out of labels trying to match \"%s\"\n",label1);
      /* Stop here - there are no more labels to process */
      err++;
      break;
    }

    if (imatch == 0) {
      /* no match */
      printf("MtzParseLabin: neither label recognised: %s %s \n",label1,label2);
      err++;
    }
  }

  /* Finished with the parser array */
  ccp4_parse_end(parser);

  return err ? -1 : nlabels;
}

MTZCOL **ccp4_lrassn(const MTZ *mtz, const char labels[][31], const int nlabels, 
             char types[][3]) 
{
  int ilab;
  char label[31];
  MTZCOL *col, **lookup;

  lookup = (MTZCOL **) ccp4_utils_malloc(nlabels*sizeof(MTZCOL *));

 /* Loop over labels */
   for (ilab = 0; ilab < nlabels; ++ilab) {

     strcpy(label,labels[ilab]);
     /* column not assigned */
     if (label[0] == '\0') {
       lookup[ilab] = NULL;
     } else {
       /* Special trap for M/ISYM */
       if (strcmp(types[ilab],"Y") == 0 && strcmp(label,"M/ISYM") == 0)
         strcpy(label,"M_ISYM");
       col = MtzColLookup(mtz,label);
       if (col != NULL) {

	 /* if requested type is blank, return actual type */
         if (!strcmp(types[ilab],"")) {
           if (strcmp(col->type,"")) {
             strcpy(types[ilab],col->type);
	   } else {
             strcpy(types[ilab],"R");
	   }

	 /* check requested column type against file type. */
	 } else if (strncmp(col->type,types[ilab],1)) {
           ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_ColTypeMismatch),"ccp4_lrassn",NULL);
           printf("   From ccp4_lrassn: expected type %s does not match file type %s for column %s\n", 
             types[ilab],col->type,col->label);
           if (!strcmp(types[ilab],"R") || !strcmp(types[ilab],"I"))
             printf("   (This may be intended for generic types R/I.) \n");
	 }
       }

       lookup[ilab] = col;
     }
   }

   return lookup;
}

int ccp4_lridx(const MTZ *mtz, const MTZSET *set, char crystal_name[64], 
            char dataset_name[64], char project_name[64], int *isets, 
            float datcell[6], float *datwave)
{
  int i;
  MTZXTAL *xtl;

  /* find which crystal this dataset  belongs to */
  xtl = MtzSetXtal(mtz, set);

  /* copy crystal and dataset information */
  strncpy(crystal_name,xtl->xname,63);
  crystal_name[63] = '\0';
  strncpy(dataset_name,set->dname,63);
  dataset_name[63] = '\0';
  strncpy(project_name,xtl->pname,63);
  project_name[63] = '\0';
  *isets = set->setid;
  for (i = 0; i < 6; ++i)
     datcell[i] = xtl->cell[i];
  *datwave = set->wavelength;

  return 1;
}

/* Return MTZ record in file order */
int ccp4_lrrefl(const MTZ *mtz, float *resol, float adata[], int logmss[], int iref) {

  int i,j,k;
  int ind[3],ixtal;
  unsigned int colin;
  float *refldata;
  double coefhkl[6];

  /* If we are past the last reflection, indicate this with return value. */
  if (iref > mtz->nref_filein) 
    return 1;

  /* If reflections not in memory, read next record from file. */
  if (!mtz->refs_in_memory) {
    refldata = (float *) ccp4_utils_malloc(mtz->ncol_read*sizeof(float));
    if (MtzRrefl( mtz->filein, mtz->ncol_read, refldata) == EOF) {
      free(refldata);
      return 1;
    }
  }

 /* Loop over all columns in the MTZ struct, and select those which
    derive from the input file. */
  for (i = 0; i < mtz->nxtal; ++i) {
    for (j = 0; j < mtz->xtal[i]->nset; ++j) {
     for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) {
       if ((colin = mtz->xtal[i]->set[j]->col[k]->source) != 0) {
         if (mtz->refs_in_memory) {
           adata[colin - 1] = mtz->xtal[i]->set[j]->col[k]->ref[iref-1];
         } else {
           adata[colin - 1] = refldata[colin - 1];
	 }
         logmss[colin - 1] = ccp4_ismnf(mtz, adata[colin - 1]);
         if (mtz->xtal[i]->set[j]->col[k]->type[0] == 'H') {
           if (strcmp(mtz->xtal[i]->set[j]->col[k]->label,"H") == 0)
             ind[0] = (int) adata[colin - 1];
           if (strcmp(mtz->xtal[i]->set[j]->col[k]->label,"K") == 0)
             ind[1] = (int) adata[colin - 1];
           if (strcmp(mtz->xtal[i]->set[j]->col[k]->label,"L") == 0)
             ind[2] = (int) adata[colin - 1];
	 }
       }
     }
    }
  }

  /* calculate resolution of this reflection, based on cell of 
     first crystal with non-zero cell dimensions */
  for (ixtal = 0; ixtal < mtz->nxtal; ++ixtal)
   if (mtz->xtal[ixtal]->cell[0] > 0.001) {
     MtzHklcoeffs(mtz->xtal[ixtal]->cell, coefhkl);
     break;
   }
  *resol = MtzInd2reso(ind, coefhkl);
  /* kludge taken from mtzlib.f */
  if (*resol > mtz->xtal[ixtal]->resmax) *resol = mtz->xtal[ixtal]->resmax;
  if (*resol < mtz->xtal[ixtal]->resmin) *resol = mtz->xtal[ixtal]->resmin;

  free(refldata);
  return 0;
}

/* Return MTZ record in lookup order */
int ccp4_lrreff(const MTZ *mtz, float *resol, float adata[], int logmss[],
   const MTZCOL *lookup[], const int ncols, const int iref) {

  int icol,l;
  int ind[3],ixtal,ind_xtal,ind_set,ind_col[3];
  unsigned int colin;
  float *refldata;
  double coefhkl[6];
  union float_uint_uchar uf;

  /* If we are past the last reflection, indicate this with return value. */
  if (iref > mtz->nref_filein) 
    return 1;

  /* If reflections not in memory, read next record from file. */
  if (!mtz->refs_in_memory) {
    refldata = (float *) ccp4_utils_malloc(mtz->ncol_read*sizeof(float));
    if (MtzRrefl( mtz->filein, mtz->ncol_read, refldata) == EOF) {
      free(refldata);
      return 1;
    }
  }

  if (strncmp (mtz->mnf.amnf,"NAN",3) == 0) {
    uf = ccp4_nan();
  } else {
    uf.f = mtz->mnf.fmnf;
  }

  /* loop over columns requested in lookup array. */
  for (icol=0; icol < ncols; icol++) {
    logmss[icol] = 1;
    if (lookup[icol]) {
      if (mtz->refs_in_memory) {
        adata[icol] = lookup[icol]->ref[iref-1];
        logmss[icol] = ccp4_ismnf(mtz, adata[icol]);
      } else {
         if ((colin = lookup[icol]->source) != 0) {
           adata[icol] = refldata[colin - 1];
           logmss[icol] = ccp4_ismnf(mtz, adata[icol]);
	 } else {
           adata[icol] = uf.f;
           logmss[icol] = 1;
	 }
      }
    }
  }

  /* Check if HKL are first 3 columns */
  if (lookup[0]->type[0] == 'H' && lookup[1]->type[0] == 'H' && 
      lookup[2]->type[0] == 'H') {
    ind[0] = (int) adata[0];
    ind[1] = (int) adata[1];
    ind[2] = (int) adata[2];
  } else {
    MtzFindInd(mtz,&ind_xtal,&ind_set,ind_col);
    for (l = 0; l < ncols; ++l) {
      if (lookup[l] == mtz->xtal[ind_xtal]->set[ind_set]->col[ind_col[0]])
        ind[0] = (int) adata[l];
      if (lookup[l] == mtz->xtal[ind_xtal]->set[ind_set]->col[ind_col[1]])
        ind[1] = (int) adata[l];
      if (lookup[l] == mtz->xtal[ind_xtal]->set[ind_set]->col[ind_col[2]])
        ind[2] = (int) adata[l];
    }
  }

  /* calculate resolution of this reflection, based on cell of 
     first crystal with non-zero cell dimensions */
  for (ixtal = 0; ixtal < mtz->nxtal; ++ixtal)
   if (mtz->xtal[ixtal]->cell[0] > 0.001) {
     MtzHklcoeffs(mtz->xtal[ixtal]->cell, coefhkl);
     break;
   }
  *resol = MtzInd2reso(ind, coefhkl);
  /* kludge taken from mtzlib.f */
  if (*resol > mtz->xtal[ixtal]->resmax) *resol = mtz->xtal[ixtal]->resmax;
  if (*resol < mtz->xtal[ixtal]->resmin) *resol = mtz->xtal[ixtal]->resmin;

  free(refldata);
  return 0;
}

void MtzRewdInput(MTZ *mtz) {
  if (mtz->filein) {
    ccp4_file_seek(mtz->filein, SIZE1, SEEK_SET);
  } else {
    printf("MtzRewdInput: No associated file. Was MtzGet called with read_refs option?\n");
  }
}

int ccp4_ismnf(const MTZ *mtz, const float datum) {

  if (strncmp (mtz->mnf.amnf,"NAN",3) == 0) {
    return ccp4_utils_isnan((union float_uint_uchar *) &datum);
  } else {
    if (datum == mtz->mnf.fmnf)
        return 1;
  }
  return 0;
}

int ccp4_lhprt(const MTZ *mtz, int iprint) {

  int i,j,k,numbat,isort[5],base_set_exists=0;
  float maxres=0.0,minres=100.0;
  char buffer[MTZRECORDLENGTH+1],symline[81];
  MTZSET *baseset=NULL;

  if (iprint <= 0) return 2;

  printf(" * Title:\n\n");
  printf(" %s\n\n",mtz->title);

  if ((baseset = MtzSetLookup(mtz,"HKL_base/HKL_base")) != NULL) {
    if ( MtzNumActiveColsInSet(baseset) ||
         MtzNbatchesInSet(mtz,baseset) ) {
      printf(" * Base dataset:\n\n");
      printf(" %8d %s\n",baseset->setid,"HKL_base");
      printf("          %s\n","HKL_base");
      printf("          %s\n","HKL_base");
      base_set_exists=1;
    }
  }

  printf("\n * Number of Datasets = %d\n\n",MtzNumActiveSet(mtz)-base_set_exists);
  printf(" * Dataset ID, project/crystal/dataset names, cell dimensions, wavelength:\n\n");
 /* Loop over crystals */
  for (i = 0; i < mtz->nxtal; ++i) {
 /* Loop over datasets for each crystal */
    for (j = 0; j < mtz->xtal[i]->nset; ++j) {
      /* is this the base dataset? */
      if (mtz->xtal[i]->set[j] == baseset) continue;
      /* check if dataset contains any active columns */
      if ( (MtzNumActiveColsInSet(mtz->xtal[i]->set[j]) == 0) &&
           (MtzNbatchesInSet(mtz,mtz->xtal[i]->set[j]) == 0) ) continue;
      printf(" %8d %s\n",mtz->xtal[i]->set[j]->setid,mtz->xtal[i]->pname);
      printf("          %s\n",mtz->xtal[i]->xname);
      printf("          %s\n",mtz->xtal[i]->set[j]->dname);
      printf("          %10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n",
        mtz->xtal[i]->cell[0],mtz->xtal[i]->cell[1],mtz->xtal[i]->cell[2],
        mtz->xtal[i]->cell[3],mtz->xtal[i]->cell[4],mtz->xtal[i]->cell[5]);
      printf("          %10.5f\n",mtz->xtal[i]->set[j]->wavelength);
    }
  }
  printf("\n * Number of Columns = %d\n\n",MtzNumActiveCol(mtz));
  printf(" * Number of Reflections = %d\n\n",mtz->nref);
  if (strncmp (mtz->mnf.amnf,"NAN",3) == 0) {
   printf(" * Missing value set to NaN in input mtz file\n\n");
  } else {
   printf(" * Missing value set to %f in input mtz file\n\n",mtz->mnf.fmnf);
  }

  /* if new batch headers have been written, lose the old ones */
  if (MtzNbat(mtz) > mtz->n_orig_bat) {
    numbat = MtzNbat(mtz) - mtz->n_orig_bat;
  } else {
    numbat = mtz->n_orig_bat;
  }
  if (numbat > 0)
    printf(" * Number of Batches = %d\n\n",numbat);

  if (iprint == 2 || iprint == 3) {
    printf(" * HISTORY for current MTZ file :\n\n");
    for (i = 0; i < mtz->histlines; ++i) {
      strncpy(buffer,mtz->hist + MTZRECORDLENGTH*i,MTZRECORDLENGTH);
      buffer[MTZRECORDLENGTH] = '\0';
      printf(" %s\n",buffer);
    }
    printf("\n");
  }

  if (iprint == 1 || iprint == 2 || iprint >=4 ) {

  printf(" * Column Labels :\n\n");
 /* Loop over crystals */
  for (i = 0; i < mtz->nxtal; ++i) {
 /* Loop over datasets for each crystal */
   for (j = 0; j < mtz->xtal[i]->nset; ++j) {
 /* Loop over columns for each dataset */
    for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) {
     if (mtz->xtal[i]->set[j]->col[k]->active) {
      if (strcmp(mtz->xtal[i]->set[j]->col[k]->type,"Y") == 0 && 
         strcmp(mtz->xtal[i]->set[j]->col[k]->label,"M_ISYM") == 0) {
       printf(" M/ISYM");
      } else {
       printf(" %s",mtz->xtal[i]->set[j]->col[k]->label);
      }
     }
    }
   }
  }
  printf("\n\n * Column Types :\n\n");
 /* Loop over crystals */
  for (i = 0; i < mtz->nxtal; ++i) {
 /* Loop over datasets for each crystal */
   for (j = 0; j < mtz->xtal[i]->nset; ++j) {
 /* Loop over columns for each dataset */
    for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) {
     if (mtz->xtal[i]->set[j]->col[k]->active) 
      printf(" %s",mtz->xtal[i]->set[j]->col[k]->type);
    }
   }
  }
  printf("\n\n * Associated datasets :\n\n");
 /* Loop over crystals */
  for (i = 0; i < mtz->nxtal; ++i) {
 /* Loop over datasets for each crystal */
   for (j = 0; j < mtz->xtal[i]->nset; ++j) {
 /* Loop over columns for each dataset */
    for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) {
     if (mtz->xtal[i]->set[j]->col[k]->active) 
      printf(" %d",mtz->xtal[i]->set[j]->setid);
    }
   }
  }

  } else if ( iprint == 3 ) {

  printf(" * Column Labels, Types, Ranges [and Dataset IDs] :\n\n");

  /* Loop over crystals/datasets/columns */
  for (i = 0; i < mtz->nxtal; ++i) {
   for (j = 0; j < mtz->xtal[i]->nset; ++j) {
    for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) {
     if (mtz->xtal[i]->set[j]->col[k]->active) {
      if (strcmp(mtz->xtal[i]->set[j]->col[k]->type,"Y") == 0 && 
         strcmp(mtz->xtal[i]->set[j]->col[k]->label,"M_ISYM") == 0) {
       printf(" M/ISYM                         %2s %19.4f %19.4f %8d \n",
         mtz->xtal[i]->set[j]->col[k]->type,
         mtz->xtal[i]->set[j]->col[k]->min,mtz->xtal[i]->set[j]->col[k]->max,
         mtz->xtal[i]->set[j]->setid);
      } else {
       printf(" %-30s %2s %19.4f %19.4f %8d \n",
         mtz->xtal[i]->set[j]->col[k]->label,mtz->xtal[i]->set[j]->col[k]->type,
         mtz->xtal[i]->set[j]->col[k]->min,mtz->xtal[i]->set[j]->col[k]->max,
         mtz->xtal[i]->set[j]->setid);
      }
     }
    }
   }
  }

  }

  /* write overall cell - just for scripts which grep for this */
  printf("\n\n * Cell Dimensions : (obsolete - refer to dataset cell dimensions above)\n\n");
  for (i = 0; i < mtz->nxtal; ++i) 
    if (mtz->xtal[i]->cell[0] > 0.001) {
      printf(" %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f \n\n",
        mtz->xtal[i]->cell[0],mtz->xtal[i]->cell[1],mtz->xtal[i]->cell[2],
        mtz->xtal[i]->cell[3],mtz->xtal[i]->cell[4],mtz->xtal[i]->cell[5]);
      break;
    }

  /* Calculate overall resolution limits. Two cases: If we have written some reflections
     to file, we probably want to know the resolution limits for these. In this case,
     mtz->resmax_out and mtz->resmin_out have been set and we use those. Otherwise, 
     we use the resolution limits of the crystals in memory. */
  if (mtz->resmax_out > 0.0001) {
    maxres = mtz->resmax_out;
    minres = mtz->resmin_out;
  } else {
    for (i = 0; i < mtz->nxtal; ++i) {
      if (mtz->xtal[i]->resmax > maxres) maxres = mtz->xtal[i]->resmax;
      if (mtz->xtal[i]->resmin < minres) minres = mtz->xtal[i]->resmin;
    }
  }
  printf(" *  Resolution Range :\n\n");
  if (maxres > 0.0 && minres > 0.0) {
    printf(" %10.5f %10.5f     ( %10.3f - %10.3f A )\n\n", 
            minres,maxres,1.0/sqrt(minres),1.0/sqrt(maxres));
  } else if (maxres > 0.0)  {
    printf(" %10.5f %10.5f     ( inf  - %10.3f A )\n\n", 
            minres,maxres,1.0/sqrt(maxres));
  } else {
    printf("   Not set - no crystals or reflections? \n\n");
  }
  ccp4_lrsort(mtz, isort);
  printf(" * Sort Order :\n\n  %5d %5d %5d %5d %5d\n\n",isort[0],isort[1],isort[2],
       isort[3],isort[4]);

  if (iprint == 3 || iprint == 4 ) {

    printf(" * Number of Symmetry Operations = %d \n",mtz->mtzsymm.nsym);
    printf(" * Number of Primitive Operations = %d \n",mtz->mtzsymm.nsymp);
    printf(" * Space Group = %d \'%s\' \n",mtz->mtzsymm.spcgrp,mtz->mtzsymm.spcgrpname);
    printf(" * Lattice Type = %c \n",mtz->mtzsymm.symtyp);
    printf(" * Point Group Name = %s \n",mtz->mtzsymm.pgname);

    printf("\n * Symmetry Operations : \n\n");
    for (i = 0; i < mtz->mtzsymm.nsym; ++i) {
      mat4_to_symop(symline,symline+80,(const float (*)[4])mtz->mtzsymm.sym[i]);
      symline[60] = '\0';
      printf(" Symmetry %d %s\n",i+1,symline);
      for (j = 0; j < 4; ++j) 
        printf(" %5.2f %5.2f %5.2f %5.2f \n",mtz->mtzsymm.sym[i][j][0],
           mtz->mtzsymm.sym[i][j][1],mtz->mtzsymm.sym[i][j][2],
	       mtz->mtzsymm.sym[i][j][3]);
    }
    printf("\n");

  } else {
    printf(" * Space group = \'%s\' (number     %d)\n\n",mtz->mtzsymm.spcgrpname,
       mtz->mtzsymm.spcgrp);
  }

  if (mtz->mtzsymm.spg_confidence == 'L') {
    printf("  (only Bravais lattice is fixed so far)\n\n");
  } else if (mtz->mtzsymm.spg_confidence == 'P') {
    printf("  (only pointgroup is fixed so far)\n\n");
  } else if (mtz->mtzsymm.spg_confidence == 'E') {
    printf("  (one of pair of enantiomorphic spacegroups)\n\n");
  } else if (mtz->mtzsymm.spg_confidence == 'S') {
    printf("  (spacegroup is known)\n\n");
  }

  return 1;
}

int ccp4_lhprt_adv(const MTZ *mtz, int iprint) {

  int i,j,k;
  char buffer[MTZRECORDLENGTH+1];

  printf(" HEADER INFORMATION FROM MTZ FILE \n\n");

  printf(" * File information :\n\n");

  printf("%s       %s\n",MTZTITLE,mtz->title);
  printf("%s       %d\n",MTZSPACEGROUP,mtz->mtzsymm.spcgrp);
  printf("%s       %d\n",MTZNUMREFLS,mtz->nref);
  if (strncmp (mtz->mnf.amnf,"NAN",3) == 0) {
   printf("%s       %s\n",MTZMNF,"NaN");
  } else {
   printf("%s       %f\n",MTZMNF,mtz->mnf.fmnf);
  }
  printf("%s       %s\n",MTZSORTORDER,"(not implemented)");

  printf("\n * Crystals, datasets :\n");
 /* Loop over crystals */
  for (i = 0; i < mtz->nxtal; ++i) {

    printf("\n%s       %s\n",CRYSTALXTALNAME,mtz->xtal[i]->xname);
    printf("%s       %s\n",CRYSTALPNAME,mtz->xtal[i]->pname);
    printf("%s       %10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n",CRYSTALCELL,
        mtz->xtal[i]->cell[0],mtz->xtal[i]->cell[1],mtz->xtal[i]->cell[2],
        mtz->xtal[i]->cell[3],mtz->xtal[i]->cell[4],mtz->xtal[i]->cell[5]);

 /* Loop over datasets for each crystal */
    for (j = 0; j < mtz->xtal[i]->nset; ++j) {
      printf("\n    %s       %s\n",DATASETDNAME,mtz->xtal[i]->set[j]->dname);
      printf("    %s       %10.5f\n",DATASETWAVELENGTH,mtz->xtal[i]->set[j]->wavelength);
      if (mtz->xtal[i]->set[j]->ncol > 0) {
        printf("\n        %s %s\n",COLUMNLABEL,COLUMNTYPE);
 /* Loop over columns for each dataset */
        for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) {
          printf("                     %-31s %-3s\n",mtz->xtal[i]->set[j]->col[k]->label,
                                       mtz->xtal[i]->set[j]->col[k]->type);
        }
      }
    }
  }

  printf("\n * HISTORY for current MTZ file :\n\n");
  for (i = 0; i < mtz->histlines; ++i) {
    strncpy(buffer,mtz->hist + MTZRECORDLENGTH*i,MTZRECORDLENGTH);
    buffer[MTZRECORDLENGTH] = '\0';
    printf(" %s\n",buffer);
  }

  return 1;
}

int ccp4_lrbat(MTZBAT *batch, float *buf, char *charbuf, int iprint)

{ int nwords=NBATCHWORDS,nintegers=NBATCHINTEGERS,nreals=NBATCHREALS;
  int *intbuf = (int *) buf;
  float *fltbuf = buf + NBATCHINTEGERS;

  if (!batch) return 0;

  MtzBatchToArray(batch,intbuf,fltbuf);
  intbuf[0] = nwords;
  intbuf[1] = nintegers;
  intbuf[2] = nreals;

  strncpy(charbuf,batch->title,70); 
  strncpy(charbuf+70,batch->gonlab[0],8); 
  strncpy(charbuf+78,batch->gonlab[1],8); 
  strncpy(charbuf+86,batch->gonlab[2],8); 

  if (iprint == 1) {
    printf(" Batch number: \n %6d    %s\n",batch->num,batch->title);
  } else if (iprint > 1) {
    MtzPrintBatchHeader(batch);
  }

  return 1;
}

int MtzPrintBatchHeader(const MTZBAT *batch) {

  int i;
  char labtype[26],axes[5],string1[40],string2[40];

  switch (batch->ldtype) {
  case 1:
    strcpy(labtype,"oscillation data");
    break;
  case 2:
    strcpy(labtype,"area detector data");
    break;
  case 3:
    strcpy(labtype,"Laue data");
    break;
  default:
    strcpy(labtype,"*** unknown data type ***");
  }

  switch (batch->jumpax) {
  case 1:
    strcpy(axes,"a*");
    break;
  case 2:
    strcpy(axes,"b*");
    break;
  case 3:
    strcpy(axes,"c*");
    break;
  default:
    strcpy(axes,"none");
  }

  printf(" Batch number: \n");
  printf(" %6d    %s\n",batch->num,batch->title);
  printf("\n %s \n\n %s %7d     %s  \n\n %s %7d\n %s %7d\n %s %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n %s %7d %7d %7d %7d %7d %7d \n",
	 "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++",
	 "Orientation data for batch",batch->num,labtype,
	 "  Crystal number ...................",batch->ncryst,
         "  Associated dataset ID ............",batch->nbsetid,
         "  Cell dimensions ..................",
         batch->cell[0],batch->cell[1],batch->cell[2],
         batch->cell[3],batch->cell[4],batch->cell[5],
         "  Cell fix flags ...................",
         batch->lbcell[0],batch->lbcell[1],batch->lbcell[2],
         batch->lbcell[3],batch->lbcell[4],batch->lbcell[5]);
  if (!batch->misflg) {
    strcpy(string1,"Orientation matrix U .............");
    strcpy(string2,"    (including setting angles)    ");
  } else {
    strcpy(string1,"Standard orientation matrix U ....");
    strcpy(string2,"                                  ");
  }    
  printf("   %s %9.4f %9.4f %9.4f \n   %s %9.4f %9.4f %9.4f \n   %s %9.4f %9.4f %9.4f \n",
         string1,batch->umat[0],batch->umat[3],batch->umat[6],
         string2,batch->umat[1],batch->umat[4],batch->umat[7],
         "                                  ",batch->umat[2],batch->umat[5],batch->umat[8]);
  if (batch->misflg == 1) {
    printf("   %s %6.2f %6.2f %6.2f\n",
         "Missetting angles PhiX PhiY PhiZ..",
         batch->phixyz[0][0],batch->phixyz[0][1],batch->phixyz[0][2]);
  } else if (batch->misflg > 1) {
    printf("   %s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
         "Missetting angles PhiX PhiY PhiZ..",
         batch->phixyz[0][0],batch->phixyz[0][1],batch->phixyz[0][2],
         batch->phixyz[1][0],batch->phixyz[1][1],batch->phixyz[1][2]);
  }
  printf("   %s%s%s   %s\n",
         "Reciprocal axis nearest ",batch->gonlab[batch->ngonax-1],"..",axes);
  if (!batch->lcrflg) {
    printf("   %s %6.3f \n",
	   "Mosaicity ........................",batch->crydat[0]);
  } else {
    printf("   %s %6.3f %6.3f \n",
           "Mosaicity (horizontal, vertical)..",batch->crydat[0],batch->crydat[1]);
  }
  printf("   Datum goniostat angles (degrees)..");
  for (i = 0; i < batch->ngonax; ++i) 
    printf(" %8.3f",batch->datum[i]);
  printf("\n");

  if (batch->jsaxs > 0 && batch->jsaxs <= batch->ngonax) 
    printf("   %s  %s \n",
	 "Scan axis ........................",batch->gonlab[batch->jsaxs-1]);
  printf("   %s %8.3f %8.3f \n   %s %8.3f \n   %s %8.2f %8.2f \n",
	 "Start & stop Phi angles (degrees).",batch->phistt,batch->phiend,
	 "Range of Phi angles (degrees).....",batch->phirange,
         "Start & stop time (minutes).......",batch->time1,batch->time2);

  if (batch->nbscal == 4) {
    printf("   %s %9.4f %9.4f \n   %s %9.4f %9.4f \n",
	   "   Batch scale & SD .................",batch->bscale,batch->sdbscale,
	   "   Batch B-factor & SD ..............",batch->bbfac,batch->sdbfac);
  }

  printf("   %s  \n   %s %7d \n   %s %s %s %9.4f %9.4f %9.4f \n   %s %s %s %9.4f %9.4f %9.4f \n   %s %s %s %9.4f %9.4f %9.4f \n",
         " Crystal goniostat information :-",
	 "   Number of goniostat axes..........",batch->ngonax,
	 "   Goniostat vectors.....",batch->gonlab[0],"....",batch->e1[0],batch->e1[1],batch->e1[2],
	 "                    .....",batch->gonlab[1],"....",batch->e2[0],batch->e2[1],batch->e2[2],
	 "                    .....",batch->gonlab[2],"....",batch->e3[0],batch->e3[1],batch->e3[2]);

  printf("   %s \n   %s  %9.4f %9.4f %9.4f \n   %s %9.4f %9.4f %9.4f \n",
	 " Beam information :-",
	 "   Idealized X-ray beam vector.......",batch->source[0],batch->source[1],batch->source[2],
	 "   X-ray beam vector with tilts......",batch->so[0],batch->so[1],batch->so[2]);

  if (batch->lbmflg == 0) {
    printf("   %s %9.5f %9.5f \n",
           "   Wavelength and dispersion ........",batch->alambd,batch->delamb);
  } else if (batch->lbmflg == 1) {
    printf("   %s %9.5f %9.5f %9.5f \n   %s %7.3f %7.3f \n",
	   "   Wavelength and dispersion ........",batch->alambd,batch->delamb,batch->delcor,
	   "   Divergence .......................",batch->divhd,batch->divvd);
  }

  printf(" Detector information :-\n   Number of detectors...............%7d \n",batch->ndet);
  printf("   %s%9.3f\n%s%9.3f\n%s%7.1f%7.1f%7.1f%7.1f\n",
         "   Crystal to Detector distance (mm).",batch->dx[0],
         "   Detector swing angle..............",batch->theta[0],
         "   Pixel limits on detector..........",batch->detlm[0][0][0],batch->detlm[0][0][1],batch->detlm[0][1][0],batch->detlm[0][1][1]);
  printf(" ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n");

  return 1;
}

int ccp4_lwtitl(MTZ *mtz, const char *ftitle, int flag) {

  int length;

  if (flag == 0) {

    strncpy(mtz->title,ftitle,70);

  } else {

    /* Append ftitle to existing title.
       BEWARE this has been fixed a few times for special
       cases. There is often a reaons behind numbers
       such as 69, so don't change it lightly */

    length = (int) strlen(mtz->title);
    /* this shouldn't happen if title is NULL terminated */
    if (length > 70) length = 70;
    while ((--length >= 0) && mtz->title[length] == ' ');
    /* if there is an existing title and it doesn't take
       up all 70 chars, then add a space before appending
       new title */
    if (length >= 0 && length < 69)
      mtz->title[++length] = ' ';
    strncpy(mtz->title+length+1,ftitle,69-length);

  }
  mtz->title[70] = '\0';

  return 1;
}

int MtzSetSortOrder(MTZ *mtz, MTZCOL *colsort[5]) {

  int i;

  for (i = 0; i < 5; ++i) 
    mtz->order[i] = colsort[i];

  return 1;
}

int MtzAddHistory(MTZ *mtz, const char history[][MTZRECORDLENGTH], const int nlines) {

  int i,j,numlines=0;
  char *newhist;

  newhist = MtzCallocHist(mtz->histlines + nlines);
  /* write new history lines */
  for (i = 0; i < nlines; ++i) {
   for (j = 0; j < MTZRECORDLENGTH; ++j) {
    /* remove leading blanks and blank lines */
    if ( *(history[i]+j) != ' ') {
     strncpy(newhist + MTZRECORDLENGTH*i,history[i]+j,MTZRECORDLENGTH-j);
     ++numlines;
     break;
    }
   }
  }
  /* copy old history lines */
  for (i = 0; i < mtz->histlines; ++i) {
    strncpy(newhist + MTZRECORDLENGTH*numlines + MTZRECORDLENGTH*i,
        mtz->hist + MTZRECORDLENGTH*i,MTZRECORDLENGTH);
  }
  MtzFreeHist(mtz->hist);
  mtz->hist = newhist;
  mtz->histlines += numlines;

  return mtz->histlines;
}

int ccp4_lwidx(MTZ *mtz, const char crystal_name[],  const char dataset_name[],
       const char project_name[], const float datcell[6], const float *datwave) {

  MTZXTAL *xtl;
  MTZSET *set;
  int i;
  char path1[200];

  /* Is it a new crystal? */
  if ((xtl = MtzXtalLookup(mtz,crystal_name)) == NULL) {
    xtl = MtzAddXtal(mtz,crystal_name,project_name,datcell);
    MtzAddDataset(mtz,xtl,dataset_name,*datwave);
  } else {
    /* Existing crystal - update parameters */
    if (project_name && strlen(project_name) > 0) {
      strncpy(xtl->pname,project_name,64);
      xtl->pname[64] = '\0';
    }
    if (datcell[0] > 0.0)
      for (i = 0; i < 6; ++i) 
        xtl->cell[i] = datcell[i];
    strcpy( path1, "/" );
    strcat( path1, xtl->xname );
    strcat( path1, "/" );
    strcat( path1, dataset_name );
    /* Is it a new dataset? */
    if ((set = MtzSetLookup(mtz,path1)) == NULL) {
      MtzAddDataset(mtz,xtl,dataset_name,*datwave);
    } else {
      if (*datwave > 0.0)
        set->wavelength = *datwave;
    }
  }
  return 1;
}

int MtzAssignHKLtoBase(MTZ *mtz)
{
  int i,j,k,l=0;
  MTZSET *baseset=NULL;
  MTZCOL *colarray[3];

  /* get base dataset if it exists */
  baseset = MtzSetLookup(mtz,"HKL_base/HKL_base");

  if (baseset) {

   for (i = 0; i < mtz->nxtal; ++i)
    for (j = 0; j < mtz->xtal[i]->nset; ++j)
     for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) 
       if ( strcmp(mtz->xtal[i]->set[j]->col[k]->type,"H") == 0 ) {
        colarray[l++] = mtz->xtal[i]->set[j]->col[k];
        if (l == 3) goto assign;
       }

    assign:
    for (l = 0; l < 3; ++l)
      if (colarray[l]) MtzAssignColumn(mtz, colarray[l], "HKL_base","HKL_base");

  }
  return 1;
}

int MtzAssignColumn(MTZ *mtz, MTZCOL *col, const char crystal_name[],  
     const char dataset_name[]) 
{

  MTZXTAL *xtl;
  MTZSET *set, *oldset;
  int i,j;
  float datcell[6] = {0.0}, datwave = 0.0;
  char path1[200], *path2;

  if ( !mtz || !col || !crystal_name || !dataset_name || 
       !strcmp(crystal_name,"") || !strcmp(dataset_name,"") )
      ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_ParamError),"MtzAssignColumn",NULL);

  /* if column already belongs in this dataset, do nothing and return */
  oldset = MtzColSet(mtz, col);
  path2 = MtzSetPath(mtz, oldset);
  strcpy( path1, "/" );
  strcat( path1, crystal_name );
  strcat( path1, "/" );
  strcat( path1, dataset_name );
  if ( MtzPathMatch( path1, path2 ) ) {
    free (path2);
    return 1;
  }
  free (path2);

  /* remove column from existing set */
  for (i = 0; i < oldset->ncol; ++i) {
    if ( oldset->col[i] == col ) {
      for (j = i; j < oldset->ncol - 1; ++j) 
	oldset->col[j] = oldset->col[j+1];
      oldset->col[--oldset->ncol] = NULL;
      break;
    }
  }

  /* Does the requested new dataset exist? If not, create it. */
  if ( !(set = MtzSetLookup(mtz,path1)) ) {
    if ( !(xtl = MtzXtalLookup(mtz,crystal_name)) ) 
      xtl = MtzAddXtal(mtz,crystal_name,crystal_name,datcell);
    set = MtzAddDataset(mtz,xtl,dataset_name,datwave);
  }

  /* Add column to new dataset */
  if ( ++set->ncol > ccp4array_size(set->col))
    ccp4array_resize(set->col, set->ncol + 9);
  set->col[set->ncol - 1] = col;

  return 1;
}

int ccp4_lwsymconf(MTZ *mtz, char spgconf[])
{
  if (spgconf[0] != ' ' && spgconf[0] != '\0') mtz->mtzsymm.spg_confidence = spgconf[0];

  return 1;
}

int ccp4_lwsymm(MTZ *mtz, int nsymx, int nsympx, float rsymx[192][4][4], 
   char ltypex[], int nspgrx, char spgrnx[], char pgnamx[])
{
  /* Could set this to "X" but beware of legacy programs where lwsymm
     still used. Don't want to overwrite flag in newer file. */
  char spgconf_temp[2]="";

  return ccp4_lwsymm_c(mtz, nsymx, nsympx, rsymx, ltypex, nspgrx, spgrnx,
		     pgnamx, spgconf_temp);
}

int ccp4_lwsymm_c(MTZ *mtz, int nsymx, int nsympx, float rsymx[192][4][4], 
		  char ltypex[], int nspgrx, char spgrnx[], char pgnamx[], 
                  char spgconf[])
{
  int i,j,k,length;

  if (nsymx > 0) {
    mtz->mtzsymm.nsym = nsymx;
    mtz->mtzsymm.nsymp = nsympx;
    for (i = 0; i < nsymx; ++i) {
      for (j = 0; j < 4; ++j) {
        for (k = 0; k < 4; ++k) {
          mtz->mtzsymm.sym[i][j][k] = rsymx[i][j][k];
        }
      }
    }
  }
  if (ltypex[0] != ' ' && ltypex[0] != '\0') mtz->mtzsymm.symtyp = ltypex[0];
  if (nspgrx != 0) mtz->mtzsymm.spcgrp = nspgrx;
  if (spgconf[0] != ' ' && spgconf[0] != '\0') mtz->mtzsymm.spg_confidence = spgconf[0];

  if (strcmp(spgrnx,"")) {
    length = ( strlen(spgrnx) < MAXSPGNAMELENGTH ) ? strlen(spgrnx) : MAXSPGNAMELENGTH;
    strncpy(mtz->mtzsymm.spcgrpname,spgrnx,length);
    mtz->mtzsymm.spcgrpname[length] = '\0';
  }
  if (strcmp(pgnamx,"")) {
    length = ( strlen(pgnamx) < MAXPGNAMELENGTH ) ? strlen(pgnamx) : MAXPGNAMELENGTH;
    strncpy(mtz->mtzsymm.pgname,pgnamx,length);
    mtz->mtzsymm.pgname[length] = '\0';
  }

  return 1;
}

MTZCOL **ccp4_lwassn(MTZ *mtz, const char labels[][31], const int nlabels, 
             const char types[][3], const int iappnd) 
{
  int i,j,k,ilab;
  MTZCOL *col, **lookup;
  MTZSET *defaultset;

  lookup = (MTZCOL **) ccp4_utils_malloc(nlabels*sizeof(MTZCOL *));

  /* if iappnd = 0, deactivate existing columns */
  if (iappnd == 0) {
 /* Loop over crystals */
   for (i = 0; i < mtz->nxtal; ++i) {
 /* Loop over datasets for each crystal */
    for (j = 0; j < mtz->xtal[i]->nset; ++j) {
 /* Loop over columns for each dataset */
     for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) {
       mtz->xtal[i]->set[j]->col[k]->active = 0;
     }
    }
   }
  } 

  /* new columns need to be assigned to a dataset. Set this
     as the base dataset if it exists, else the first dataset. */
  if ( !(defaultset = MtzSetLookup(mtz,"HKL_base/HKL_base")) )
    defaultset = mtz->xtal[0]->set[0];

  /* Loop over labels */
  for (ilab = 0; ilab < nlabels; ++ilab) {
    if (strcmp(types[ilab],"Y") == 0 && strcmp(labels[ilab],"M/ISYM") == 0) {
      col = MtzColLookup(mtz,"M_ISYM");
    } else {
      col = MtzColLookup(mtz,labels[ilab]);
    }
    if (col) {
      col->active = 1;
      lookup[ilab] = col;
    } else {
      /* add new column to first dataset - MtzAssignColumn corrects this */
      if (strcmp(types[ilab],"Y") == 0 && strcmp(labels[ilab],"M/ISYM") == 0) {
        lookup[ilab] = MtzAddColumn(mtz, defaultset, 
                      "M/ISYM", types[ilab]);
      } else {
        lookup[ilab] = MtzAddColumn(mtz, defaultset, 
                      labels[ilab], types[ilab]);
      }
    }
  }

  return lookup;
}

int ccp4_lwbat(MTZ *mtz, MTZBAT *batch, const int batno, const float *buf, const char *charbuf)

{  
  int *intbuf = (int *) buf;
  const float *fltbuf = buf + NBATCHINTEGERS;
  char cbatch[95]=" ";
  int i,cbatch_len;
  MTZBAT *otherbat;

  if (batch == NULL) {
    /* add new batch at end of list */
    batch = mtz->batch;
    /* is this the first ever batch? */
    if (batch == NULL) {
      mtz->batch = MtzMallocBatch();
      batch = mtz->batch;
      batch->num = batno;
      batch->next = NULL;
    } else {
      /* first, skip over n_orig_bat batches if some were read in */
      for (i=0; i < mtz->n_orig_bat - 1; ++i)
        batch = batch->next;
      if (mtz->n_orig_bat == 0 && batch->num == batno) {
        printf("From ccp4_lwbat: warning: attempt to add new batch with existing batch number %d!\n",batno);
        return 0;
      }
      while (batch->next != NULL) {
        batch = batch->next;
        if (batch->num == batno) {
          printf("From ccp4_lwbat: warning: attempt to add new batch with existing batch number %d!\n",batno);
          return 0;
        }
      }
      batch->next = MtzMallocBatch();
      batch = batch->next;
      batch->num = batno;
      batch->next = NULL;
    }
  } else {
    if (batch->num != batno) {
      /* renumbering - check unique */
      otherbat = mtz->batch;
      while (otherbat != NULL) {
        if (otherbat->num == batno && otherbat != batch) {
          printf("From ccp4_lwbat: warning: attempt to change batch number to existing batch number %d!\n",batno);
          return 0;
        }
        otherbat = otherbat->next;
      }
      batch->num = batno;
    }
  }

  MtzArrayToBatch(intbuf,fltbuf,batch);

  cbatch_len = ( strlen(charbuf) < 94 ) ? strlen(charbuf) : 94;
  strncpy(cbatch,charbuf,cbatch_len);

  strncpy(batch->title,cbatch,70); 
  strncpy(batch->gonlab[0],cbatch+70,8); 
  strncpy(batch->gonlab[1],cbatch+78,8); 
  strncpy(batch->gonlab[2],cbatch+86,8); 
  batch->gonlab[0][8] = batch->gonlab[1][8] = batch->gonlab[2][8] = '\0';

  return 1;
}

int ccp4_lwbsetid(MTZ *mtz, MTZBAT *batch, const char xname[], const char dname[])

{
  MTZXTAL *xtl;
  MTZSET *set;
  char path1[200];

  if ((xtl = MtzXtalLookup(mtz,xname)) != NULL) {
    strcpy( path1, "/" );
    strcat( path1, xtl->xname );
    strcat( path1, "/" );
    strcat( path1, dname );
    if ((set = MtzSetLookup(mtz,path1)) != NULL) {
      batch->nbsetid = set->setid;
      return 1;
    }
  }

  printf("From ccp4_lwbsetid: warning: dataset id not found!\n");
  return 0;
}

int MtzDeleteRefl(MTZ *mtz, int iref)

{
  int i,j,k;

  /* only possible if reflections in memory */
  if (mtz->refs_in_memory) {
    for (i = 0; i < mtz->nxtal; ++i) 
     for (j = 0; j < mtz->xtal[i]->nset; ++j) 
      for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) 
        ccp4array_delete_ordered(mtz->xtal[i]->set[j]->col[k]->ref,iref);
    --mtz->nref;
  }

  return 1;
}

int ccp4_lwrefl(MTZ *mtz, const float adata[], MTZCOL *lookup[], 
           const int ncol, const int iref)

{ int i,j,k,l,icol,ind[3],ind_xtal,ind_set,ind_col[3];
  float refldata[MCOLUMNS],res;
  double coefhkl[6];

  /* if this is extra reflection, check memory for in-memory mode */
  if (mtz->refs_in_memory && iref > mtz->nref) {
    if (iref > ccp4array_size(lookup[0]->ref)) {
     /* Loop over crystals */
      for (i = 0; i < mtz->nxtal; ++i) {
     /* Loop over datasets for each crystal */
       for (j = 0; j < mtz->xtal[i]->nset; ++j) {
      /* Loop over columns for each dataset */
        for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) {
         ccp4array_resize(mtz->xtal[i]->set[j]->col[k]->ref, iref);
        }
       }
      }
    }
  }

  /* update variables held in memory */
  icol = -1;
  for (i = 0; i < ncol; ++i) {
    if (lookup[i]) {
      /* update reflection for in-memory mode */
      if (mtz->refs_in_memory) {
        lookup[i]->ref[iref-1] = adata[i];
      } 
      /* update column ranges */
      if (iref == 1) {
        lookup[i]->min = FLT_MAX;
        lookup[i]->max = -FLT_MAX;
      }
      if (!ccp4_ismnf(mtz, adata[i])) {
        if (adata[i] < lookup[i]->min) lookup[i]->min = adata[i];
        if (adata[i] > lookup[i]->max) lookup[i]->max = adata[i];
      }
    }
  }

  /* write reflection for on-disk mode */
  if (!mtz->refs_in_memory) {

    icol = -1;
    /* Loop over all active columns */
    for (i = 0; i < mtz->nxtal; ++i) 
     for (j = 0; j < mtz->xtal[i]->nset; ++j) 
      for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) 
        if (mtz->xtal[i]->set[j]->col[k]->active) {
          ++icol;
          /* for each active column, see if value to write */
          for (l = 0; l < ncol; ++l) 
            if (lookup[l] == mtz->xtal[i]->set[j]->col[k]) {
              refldata[icol] = adata[l];
	      break;
	    }
	}

    if (MtzWrefl(mtz->fileout, icol+1, refldata) != icol+1 )
      return 0;

     /* Update resolution limits. For in-memory mode, this is done in MtzPut. */
     /* Check if HKL are first 3 columns */
     if (lookup[0]->type[0] == 'H' && lookup[1]->type[0] == 'H' && 
         lookup[2]->type[0] == 'H') {
       ind[0] = (int) adata[0];
       ind[1] = (int) adata[1];
       ind[2] = (int) adata[2];
     } else {
       MtzFindInd(mtz,&ind_xtal,&ind_set,ind_col);
       for (l = 0; l < ncol; ++l) {
         if (lookup[l] == mtz->xtal[ind_xtal]->set[ind_set]->col[ind_col[0]])
           ind[0] = (int) adata[l];
         if (lookup[l] == mtz->xtal[ind_xtal]->set[ind_set]->col[ind_col[1]])
           ind[1] = (int) adata[l];
         if (lookup[l] == mtz->xtal[ind_xtal]->set[ind_set]->col[ind_col[2]])
           ind[2] = (int) adata[l];
       }
     }

     for (i = 0; i < mtz->nxtal; ++i) {
      if (mtz->xtal[i]->cell[0] > 0.001) {
       MtzHklcoeffs(mtz->xtal[i]->cell, coefhkl);
       res = MtzInd2reso(ind, coefhkl);
       if (res > 0.0) {
         if (res > mtz->xtal[i]->resmax) mtz->xtal[i]->resmax = res;
         if (res < mtz->xtal[i]->resmin) mtz->xtal[i]->resmin = res;
         if (res > mtz->resmax_out) mtz->resmax_out = res;
         if (res < mtz->resmin_out) mtz->resmin_out = res;
       }
      }
     }
  }

  /* increment nref if we are adding new reflections */
  if (iref > mtz->nref)
    mtz->nref = iref;

  return 1;
}

int MtzPut(MTZ *mtz, const char *logname)

{ char hdrrec[81],symline[81],spgname[MAXSPGNAMELENGTH+3];
 CCP4File *fileout;
 int i, j, k, l, hdrst, icol, numbat, isort[5], debug=0;
 int ind[3],ind_xtal,ind_set,ind_col[3],length,glob_cell_written=0;
 double coefhkl[6];
 float res,refldata[MCOLUMNS];
 int nwords=NBATCHWORDS,nintegers=NBATCHINTEGERS,nreals=NBATCHREALS;
 float buf[NBATCHWORDS];
 int *intbuf = (int *) buf;
 float *fltbuf = buf + NBATCHINTEGERS;
 MTZBAT *batch, *lastoldbatch;
 MTZXTAL *xtl;
 char colsource[37], *taskenv;
 int date3[3], time3[3];

 if (debug) 
   printf(" MtzPut: entering \n");

 /* get data to fill out column source information */
 taskenv = getenv( "CCP4_TASK_ID" );
 if ( taskenv != NULL ) {
   strncpy( colsource, taskenv, 36 );
   colsource[36] = '\0';
 } else {
   ccp4_utils_idate( date3 );
   ccp4_utils_itime( time3 );
   sprintf( colsource, "CREATED_%02d/%02d/%04d_%02d:%02d:%02d",
	    date3[0],date3[1],date3[2],time3[0],time3[1],time3[2] );
 }
 for ( i = 0; i < strlen(colsource); i++ )
   if ( colsource[i] == ' ' ) colsource[i] = '_';
 for (i = 0; i < mtz->nxtal; ++i)
   for (j = 0; j < mtz->xtal[i]->nset; ++j)
     for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k)
       if ( mtz->xtal[i]->set[j]->col[k]->source == 0 )
	 strncpy(mtz->xtal[i]->set[j]->col[k]->colsource,colsource,36);

 if (!mtz->fileout) {

   if ( !(fileout = MtzOpenForWrite(logname)) ) return 0;

   if (debug) 
     printf(" MtzPut: file opened \n");

 } else {
   fileout = mtz->fileout;
 }

 if (mtz->refs_in_memory) {
   /* Write all reflections from memory - make this optional? */
   for (l = 0; l < mtz->nref; ++l) {
     icol = 0;
   /* Loop over crystals */
     for (i = 0; i < mtz->nxtal; ++i) {
   /* Loop over datasets for each crystal */
      for (j = 0; j < mtz->xtal[i]->nset; ++j) {
   /* Loop over columns for each dataset */
       for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) {
         if (mtz->xtal[i]->set[j]->col[k]->active) {
           refldata[icol++] = mtz->xtal[i]->set[j]->col[k]->ref[l];
         }
       }
      }
     }
     if (MtzWrefl(fileout, icol, refldata) != icol ) return 0;
   }

   if (debug) 
     printf(" MtzPut: reflections written \n");

 }

 ccp4_file_setmode(fileout,0);
 /* Write header */
 sprintf(hdrrec,"VERS MTZ:V%d.%d",MTZ_MAJOR_VERSN,MTZ_MINOR_VERSN);
 /* if MTZ_MAJOR_VERSN,MTZ_MINOR_VERSN get into double figures,
    adjust following call to MtzWhdrLine */
 MtzWhdrLine(fileout,13,hdrrec);
 strcpy(hdrrec,"TITLE ");
 strncpy(hdrrec+6,mtz->title,70);
 MtzWhdrLine(fileout,76,hdrrec);
 /* if new batch headers have been written, lose the old ones */
 /* mtz->n_orig_bat is original number of batches, MtzNbat(mtz) the current */
 if (MtzNbat(mtz) == mtz->n_orig_bat) {
   numbat = mtz->n_orig_bat;
 } else {
   numbat = MtzNbat(mtz) - mtz->n_orig_bat;
 }
 sprintf(hdrrec,"NCOL %8d %12d %8d",MtzNumActiveCol(mtz),mtz->nref,numbat);
 MtzWhdrLine(fileout,35,hdrrec);
 if (debug) printf(" MtzPut: NCOL just written \n");

 /* Purely for backwards compatibility: output first non-zero cell as
    global cell. Also update base dataset cell. */
 for (i = 0; i < mtz->nxtal; ++i) {
   if ( !strcmp(mtz->xtal[i]->xname,"HKL_base") ) continue;
   if ( (MtzNumActiveSetsInXtal(mtz,mtz->xtal[i]) == 0) ) continue;
   if (mtz->xtal[i]->cell[0] > 0.001) {
     sprintf(hdrrec,"CELL  %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f",mtz->xtal[i]->cell[0],
           mtz->xtal[i]->cell[1],mtz->xtal[i]->cell[2],mtz->xtal[i]->cell[3],
           mtz->xtal[i]->cell[4],mtz->xtal[i]->cell[5]);
     MtzWhdrLine(fileout,65,hdrrec);
     if ((xtl = MtzXtalLookup(mtz,"HKL_base")) != NULL)
       for (j = 0; j < 6; ++j)
         xtl->cell[j] = mtz->xtal[i]->cell[j];
     glob_cell_written=1;
     break;
   }
 }
 /* if no suitable cell found, then try HKL_base cell */
 if (!glob_cell_written) {
   if ((xtl = MtzXtalLookup(mtz,"HKL_base")) != NULL) {
     sprintf(hdrrec,"CELL  %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f",xtl->cell[0],
           xtl->cell[1],xtl->cell[2],xtl->cell[3],xtl->cell[4],xtl->cell[5]);
     MtzWhdrLine(fileout,65,hdrrec);
     glob_cell_written=1;
   }
 }
 if (debug) printf(" MtzPut: CELL just written \n");

 ccp4_lrsort(mtz, isort);
 sprintf(hdrrec,"SORT  %3d %3d %3d %3d %3d",isort[0],isort[1],isort[2],
       isort[3],isort[4]);
 MtzWhdrLine(fileout,25,hdrrec);
 if (debug) printf(" MtzPut: SORT just written \n");

 spgname[0] = '\'';
 length = strlen(mtz->mtzsymm.spcgrpname);
 while ((--length >= 0) && mtz->mtzsymm.spcgrpname[length] == ' ');
 strncpy(spgname+1,mtz->mtzsymm.spcgrpname,length+1);
 spgname[length+2] = '\'';
 spgname[length+3] = '\0';
 sprintf(hdrrec,"SYMINF %3d %2d %c %5d %22s %5s %c",mtz->mtzsymm.nsym,
         mtz->mtzsymm.nsymp,mtz->mtzsymm.symtyp,mtz->mtzsymm.spcgrp,spgname,
         mtz->mtzsymm.pgname,mtz->mtzsymm.spg_confidence);
 MtzWhdrLine(fileout,52,hdrrec);
 if (debug) printf(" MtzPut: SYMINF just written \n");

 for (i = 0; i < mtz->mtzsymm.nsym; ++i) {
     mat4_to_symop(symline,symline+74,(const float (*)[4])mtz->mtzsymm.sym[i]);
     symline[74] = '\0';
     sprintf(hdrrec,"SYMM %74s",symline);
     MtzWhdrLine(fileout,79,hdrrec); 
 }
 if (debug) printf(" MtzPut: symmetry just written \n");

 if (mtz->refs_in_memory) {
  /* Find dataset of indices */
   MtzFindInd(mtz,&ind_xtal,&ind_set,ind_col);

  /* Recalculate crystal resolution limits */
  for (i = 0; i < mtz->nxtal; ++i) {
   mtz->xtal[i]->resmax = 0.0;
   mtz->xtal[i]->resmin = 100.0;
   MtzHklcoeffs(mtz->xtal[i]->cell, coefhkl);
   for (j = 0; j < mtz->nref; ++j) {
      ind[0] = (int) mtz->xtal[ind_xtal]->set[ind_set]->col[ind_col[0]]->ref[j];
      ind[1] = (int) mtz->xtal[ind_xtal]->set[ind_set]->col[ind_col[1]]->ref[j];
      ind[2] = (int) mtz->xtal[ind_xtal]->set[ind_set]->col[ind_col[2]]->ref[j];
      res = MtzInd2reso(ind, coefhkl);
      /* crystal limits */
      if (res > 0.0) {
        if (res > mtz->xtal[i]->resmax) mtz->xtal[i]->resmax = res;
        if (res < mtz->xtal[i]->resmin) mtz->xtal[i]->resmin = res;
        if (res > mtz->resmax_out) mtz->resmax_out = res;
        if (res < mtz->resmin_out) mtz->resmin_out = res;
      }
   }
  }
 }
 /* print enough digits to retain precision. C. Flensburg 20080227 */
 sprintf(hdrrec,"RESO %-20.16f %-20.16f",mtz->resmin_out,mtz->resmax_out);
 MtzWhdrLine(fileout,46,hdrrec);

 if (debug) 
   printf(" MtzPut: resolution limts just written \n");

 if (strncmp (mtz->mnf.amnf,"NAN",3) == 0) {
   sprintf(hdrrec,"VALM NAN");
   MtzWhdrLine(fileout,8,hdrrec);
 } else {
   sprintf(hdrrec,"VALM %-20f",mtz->mnf.fmnf);
   MtzWhdrLine(fileout,25,hdrrec);
 }

 if (debug) 
   printf(" MtzPut: VALM just written \n");

 /* Loop over crystals */
 for (i = 0; i < mtz->nxtal; ++i) {
 /* Loop over datasets for each crystal */
  for (j = 0; j < mtz->xtal[i]->nset; ++j) {
 /* Loop over columns for each dataset */
   for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) {
     if (mtz->xtal[i]->set[j]->col[k]->active) {
       if (strcmp(mtz->xtal[i]->set[j]->col[k]->type,"Y") == 0 && 
           strcmp(mtz->xtal[i]->set[j]->col[k]->label,"M_ISYM") == 0) {
         sprintf(hdrrec,"COLUMN %-30s ","M/ISYM");
       } else {
         sprintf(hdrrec,"COLUMN %-30s ",mtz->xtal[i]->set[j]->col[k]->label);
       }
       /* Check that the column type is set
	  If it is blank then the COLUMN record will be incomplete */
       if (mtz->xtal[i]->set[j]->col[k]->type[0] == '\0') {
         if (ccp4_liberr_verbosity(-1))
  	   printf("From MtzPut: column type for %s is not set, assume type R\n",
		mtz->xtal[i]->set[j]->col[k]->label);
	 strncpy(mtz->xtal[i]->set[j]->col[k]->type,"R",2);
       }
       /* catch case when min and max have not been set*/
       if ( mtz->xtal[i]->set[j]->col[k]->min == FLT_MAX ) mtz->xtal[i]->set[j]->col[k]->min = 0.0f;
       if ( mtz->xtal[i]->set[j]->col[k]->max == -FLT_MAX ) mtz->xtal[i]->set[j]->col[k]->max = 0.0f;

       sprintf(hdrrec+38,"%c %17.9g %17.9g %4d",
                   mtz->xtal[i]->set[j]->col[k]->type[0],
                   mtz->xtal[i]->set[j]->col[k]->min,
                   mtz->xtal[i]->set[j]->col[k]->max,
                   mtz->xtal[i]->set[j]->setid);
       MtzWhdrLine(fileout,MTZRECORDLENGTH,hdrrec);

       if ( mtz->xtal[i]->set[j]->col[k]->colsource[0] != '\0' ) {
         if (strcmp(mtz->xtal[i]->set[j]->col[k]->type,"Y") == 0 && 
           strcmp(mtz->xtal[i]->set[j]->col[k]->label,"M_ISYM") == 0) {
	   sprintf(hdrrec,"COLSRC %-30s %-36s  %4d","M/ISYM",mtz->xtal[i]->set[j]->col[k]->colsource,mtz->xtal[i]->set[j]->setid);
         } else {
	   sprintf(hdrrec,"COLSRC %-30s %-36s  %4d",mtz->xtal[i]->set[j]->col[k]->label,mtz->xtal[i]->set[j]->col[k]->colsource,mtz->xtal[i]->set[j]->setid);
         }
	 MtzWhdrLine(fileout,MTZRECORDLENGTH,hdrrec);
       }

       if ( mtz->xtal[i]->set[j]->col[k]->grpname[0] != '\0' &&
	    mtz->xtal[i]->set[j]->col[k]->grptype[0] != '\0' &&
	    mtz->xtal[i]->set[j]->col[k]->grpposn    >=   0  ) {
         if (strcmp(mtz->xtal[i]->set[j]->col[k]->type,"Y") == 0 && 
           strcmp(mtz->xtal[i]->set[j]->col[k]->label,"M_ISYM") == 0) {
	   sprintf(hdrrec,"COLGRP %-30s %-30s %-4s %1X %4d","M/ISYM",mtz->xtal[i]->set[j]->col[k]->grpname,mtz->xtal[i]->set[j]->col[k]->grptype,mtz->xtal[i]->set[j]->col[k]->grpposn,mtz->xtal[i]->set[j]->setid);
         } else {
	   sprintf(hdrrec,"COLGRP %-30s %-30s %-4s %1X %4d",mtz->xtal[i]->set[j]->col[k]->label,mtz->xtal[i]->set[j]->col[k]->grpname,mtz->xtal[i]->set[j]->col[k]->grptype,mtz->xtal[i]->set[j]->col[k]->grpposn,mtz->xtal[i]->set[j]->setid);
         }
	 MtzWhdrLine(fileout,MTZRECORDLENGTH,hdrrec);
       }
     }
   }
  }
 }

 if (debug) 
   printf(" MtzPut: column info just written \n");

 sprintf(hdrrec,"NDIF %8d",MtzNumActiveSet(mtz));
 MtzWhdrLine(fileout,13,hdrrec);

 if (debug) 
   printf(" MtzPut: about to write dataset info \n");

 /* Loop over crystals */
 for (i = 0; i < mtz->nxtal; ++i) {
 /* Loop over datasets for each crystal */
  for (j = 0; j < mtz->xtal[i]->nset; ++j) {
   /* check if dataset contains any active columns or batches */
   if ( (MtzNumActiveColsInSet(mtz->xtal[i]->set[j]) == 0) &&
        (MtzNbatchesInSet(mtz,mtz->xtal[i]->set[j]) == 0) ) continue;
   sprintf(hdrrec,"PROJECT %7d %-64s",mtz->xtal[i]->set[j]->setid,
                                      mtz->xtal[i]->pname);
   MtzWhdrLine(fileout,MTZRECORDLENGTH,hdrrec);
   sprintf(hdrrec,"CRYSTAL %7d %-64s",mtz->xtal[i]->set[j]->setid,
                                      mtz->xtal[i]->xname);
   MtzWhdrLine(fileout,MTZRECORDLENGTH,hdrrec);
   sprintf(hdrrec,"DATASET %7d %-64s",mtz->xtal[i]->set[j]->setid,
                                      mtz->xtal[i]->set[j]->dname);
   MtzWhdrLine(fileout,MTZRECORDLENGTH,hdrrec);
   sprintf(hdrrec,"DCELL   %7d %10.4f%10.4f%10.4f%10.4f%10.4f%10.4f",
        mtz->xtal[i]->set[j]->setid,mtz->xtal[i]->cell[0],
        mtz->xtal[i]->cell[1],mtz->xtal[i]->cell[2],
        mtz->xtal[i]->cell[3],mtz->xtal[i]->cell[4], 
        mtz->xtal[i]->cell[5]);
   MtzWhdrLine(fileout,76,hdrrec);
   sprintf(hdrrec,"DWAVEL  %7d %10.5f",mtz->xtal[i]->set[j]->setid,
                                       mtz->xtal[i]->set[j]->wavelength);
   MtzWhdrLine(fileout,26,hdrrec);
  }
 }

 if (MtzNbat(mtz) > 0) {
   batch = mtz->batch;
   /* if new batch headers have been written, lose the old ones */
   if (MtzNbat(mtz) > mtz->n_orig_bat) {
     for (i=0; i < mtz->n_orig_bat; ++i) {
       lastoldbatch = batch;
       batch = batch->next;
     }
     numbat = MtzNbat(mtz) - mtz->n_orig_bat;
     batch = sort_batches(batch,numbat);
     if (mtz->n_orig_bat > 0) {
       lastoldbatch->next = batch;
     } else {
       mtz->batch = batch;
     }
   } else {
     numbat = mtz->n_orig_bat;
   }
   if (debug) {
     printf(" MtzPut: original number of batches %d \n",mtz->n_orig_bat);
     printf(" MtzPut: total number of batches %d \n",MtzNbat(mtz));
     printf(" MtzPut: number of batches to be written %d \n",numbat);
   }

   for (i = 0; i < numbat; i += 12) {
     sprintf(hdrrec,"BATCH ");
     l = 6;
     for (j = 0; j < 12 && i+j < numbat; ++j) {
       sprintf(hdrrec+6+6*j,"%6d",batch->num);
       l += 6;
       batch = batch->next;
     }
     MtzWhdrLine(fileout,l,hdrrec);
   }
 }

 /* write out unrecognized headers */
 if ( mtz->unknown_headers )
   for (i = 0; i < mtz->n_unknown_headers; ++i)
     MtzWhdrLine(fileout,MTZRECORDLENGTH,mtz->unknown_headers+i*MTZRECORDLENGTH);

 sprintf(hdrrec,"END ");
 MtzWhdrLine(fileout,4,hdrrec);

 if (debug) 
   printf(" MtzPut: main header written \n");

 if (mtz->histlines > 0) {
   sprintf(hdrrec,"MTZHIST %3d",mtz->histlines);
   MtzWhdrLine(fileout,11,hdrrec);
   for (i = 0; i < mtz->histlines; ++i) {
     strncpy(hdrrec,mtz->hist + MTZRECORDLENGTH*i,MTZRECORDLENGTH);
     MtzWhdrLine(fileout,MTZRECORDLENGTH,hdrrec);
   }
 }

 if (MtzNbat(mtz) > 0) {
   batch = mtz->batch;
   /* if new batch headers have been written, lose the old ones */
   if (MtzNbat(mtz) > mtz->n_orig_bat)
     for (i=0; i < mtz->n_orig_bat; ++i)
       batch = batch->next;
   sprintf(hdrrec,"MTZBATS");
   MtzWhdrLine(fileout,7,hdrrec);
   while (batch != NULL) {
     sprintf(hdrrec,"BH %8d%8d%8d%8d",batch->num,nwords,nintegers,nreals);
     MtzWhdrLine(fileout,35,hdrrec);
     strcpy(hdrrec,"TITLE ");
     strncpy(hdrrec+6,batch->title,70);
     MtzWhdrLine(fileout,76,hdrrec);
     MtzBatchToArray(batch,intbuf,fltbuf);
     intbuf[0] = nwords;
     intbuf[1] = nintegers;
     intbuf[2] = nreals;
     ccp4_file_setmode(fileout,2);
     ccp4_file_write(fileout, (uint8 *) buf, nwords);
     ccp4_file_setmode(fileout,0);
     if (batch->gonlab[0][0] != '\0') {
       sprintf(hdrrec,"BHCH %8s%8s%8s",batch->gonlab[0],batch->gonlab[1],batch->gonlab[2]);
     } else {
       sprintf(hdrrec,"BHCH                         ");
     }
     MtzWhdrLine(fileout,29,hdrrec);
     batch = batch->next;
   }
 }

 if (debug) 
   printf(" MtzPut: batch headers written \n");

 sprintf(hdrrec,"MTZENDOFHEADERS ");
 MtzWhdrLine(fileout,16,hdrrec);

 /* write XML data block */
 if ( mtz->xml != NULL ) {
   ccp4_file_setmode(fileout,0);
   ccp4_file_writechar(fileout,(const uint8 *)mtz->xml,strlen(mtz->xml));
 }

 /* go back and correct hdrst */
 ccp4_file_setmode(fileout,0);
 ccp4_file_seek(fileout, 4, SEEK_SET); 
 hdrst = mtz->nref * MtzNumActiveCol(mtz) + SIZE1 + 1;
 ccp4_file_setmode(fileout,2);
 ccp4_file_write(fileout,(uint8 *) &hdrst,1);

 /* And close the mtz file: */
 if (!mtz->fileout) 
   ccp4_file_close(fileout);

 if (debug) 
   printf(" MtzPut: bye bye \n");

 return 1;
}

MTZBAT *sort_batches(MTZBAT *batch, int numbat)

{ int debug=0; 
 int i,max_num_bat,isort=0,nmerges,sublistsize=1,sublist1,sublist2;
 MTZBAT *cur_batch1, *cur_batch2, *sorted_batch, *tail;
 MTZBAT *tmp;

 /* first check if already in order */
 cur_batch1 = batch;
 max_num_bat = cur_batch1->num;
 for (i=0; i < numbat; ++i) {
   cur_batch1 = cur_batch1->next;
   /* reached end of list */
   if (!cur_batch1) return batch;
   if (cur_batch1->num < max_num_bat) {
     isort=1;
     break;
   } else {
     max_num_bat = cur_batch1->num;
   }
 }
 if (!isort) return batch;

 if (ccp4_liberr_verbosity(-1))
   printf("\n Note: Sorting batch headers prior to writing to file... \n\n");

 /* Sort */
 /* This is Simon Tatham's algorithm, implemented for batches. */

 if (debug) {
   tmp = batch;
   for (i=0; i < numbat; ++i) {
     printf(" %d",tmp->num);
     tmp = tmp->next;
   }
   printf(" \n");
 }

 while (1) {

   if (debug) printf(" sort_batches: pass with sublist size %d \n",sublistsize);
   cur_batch1 = batch;
   tail = NULL;
   batch = NULL;
   nmerges = 0; 
   while (cur_batch1) {

     ++nmerges;
     cur_batch2 = cur_batch1;
     sublist1 = 0;
     while (cur_batch2 && sublist1 < sublistsize) {
       ++sublist1;
       cur_batch2 = cur_batch2->next;
     } 
     sublist2 = sublistsize;

     while (sublist1 > 0 || (sublist2 > 0 && cur_batch2)) {

       /* decide whether next batch comes from cur_batch1 or cur_batch2 */
       if (sublist1 == 0) {
	 /* cur_batch1 is empty; batch must come from cur_batch2. */
	 sorted_batch = cur_batch2; cur_batch2 = cur_batch2->next; sublist2--;
       } else if (sublist2 == 0 || !cur_batch2) {
	 /* cur_batch2 is empty; batch must come from cur_batch1. */
	 sorted_batch = cur_batch1; cur_batch1 = cur_batch1->next; sublist1--;
       } else if (cur_batch1->num <= cur_batch2->num ) {
	 /* cur_batch1 number is lower (or same); batch must come from cur_batch1. */
	 sorted_batch = cur_batch1; cur_batch1 = cur_batch1->next; sublist1--;
       } else {
	 /* cur_batch2 number is lower; batch must come from cur_batch2. */
	 sorted_batch = cur_batch2; cur_batch2 = cur_batch2->next; sublist2--;
       }

       /* add the next element to the merged list */
       if (tail) { 
         tail->next = sorted_batch;
       } else {
         batch = sorted_batch;
       }
       tail = sorted_batch;

     }

     /* sorted this sub-list - move to next */
     cur_batch1 = cur_batch2;

   }

   tail->next = NULL;

   if (debug) {
     tmp = batch;
     for (i=0; i < numbat; ++i) {
       printf(" %d",tmp->num);
       tmp = tmp->next;
     }
     printf(" \n");
   }

   /* If we have done only one merge, we're finished. */
   if (nmerges <= 1)   /* allow for nmerges==0, the empty list case */
     return batch;

    /* Otherwise repeat, merging lists twice the size */
    sublistsize *= 2;
 }

}

CCP4File *MtzOpenForWrite(const char *logname)

{ CCP4File *fileout;
 int debug=0;
 int hdrst;
 char *filename;

 if (debug) printf(" MtzOpenForWrite: entering \n");
 /* Open the mtz file: */
 if (getenv(logname) != NULL) {
   filename = strdup(getenv(logname));
 } else {
   filename = strdup(logname);
 }
 fileout = ccp4_file_open(filename,O_RDWR | O_TRUNC);
 if (! fileout ) {
   ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_CantOpenFile),"MtzOpenForWrite",NULL);
   return NULL;
 }
 if (debug) printf(" MtzOpenForWrite: file opened \n");

 /* Write initial info */
 ccp4_file_setmode(fileout,0);
 ccp4_file_writechar(fileout, (uint8 *) "MTZ ",4);
 ccp4_file_setmode(fileout,2);
 hdrst = SIZE1 + 1;
 ccp4_file_write(fileout,(uint8 *) &hdrst,1);

 ccp4_file_setstamp(fileout,2);
/* Write architecture */
 ccp4_file_warch(fileout);
 if (debug) printf(" MtzOpenForWrite: stamp written \n");

 /* Position at start of reflections - intervening gap should be filled
    with zeros */
 ccp4_file_seek(fileout, SIZE1, SEEK_SET); 

 free(filename); 

 if (debug) printf(" MtzOpenForWrite: bye bye \n");
 return fileout;
}

int MtzBatchToArray(MTZBAT *batch, int *intbuf, float *fltbuf)

{  int i;

  if (!batch) return 0;

  for (i = 0; i < NBATCHINTEGERS; ++i)
    intbuf[i] = 0;
  for (i = 0; i < NBATCHREALS; ++i)
    fltbuf[i] = 0.0;

  intbuf[3] = batch->iortyp;
  for (i = 0; i < 6; ++i)
    intbuf[4+i] = batch->lbcell[i];
  intbuf[10] = batch->misflg;
  intbuf[11] = batch->jumpax;
  intbuf[12] = batch->ncryst;
  intbuf[13] = batch->lcrflg;
  intbuf[14] = batch->ldtype;
  intbuf[15] = batch->jsaxs;
  intbuf[16] = batch->nbscal;
  intbuf[17] = batch->ngonax;
  intbuf[18] = batch->lbmflg;
  intbuf[19] = batch->ndet;
  intbuf[20] = batch->nbsetid;

  for (i = 0; i < 6; ++i)
    fltbuf[i] = batch->cell[i];
  for (i = 0; i < 9; ++i)
    fltbuf[6 + i] = batch->umat[i];
  for (i = 0; i < 3; ++i) 
    fltbuf[15 + i] = batch->phixyz[0][i];
  for (i = 0; i < 3; ++i) 
    fltbuf[18 + i] = batch->phixyz[1][i];
  for (i = 0; i < 12; ++i) 
    fltbuf[21 + i] = batch->crydat[i];
  for (i = 0; i < 3; ++i) 
    fltbuf[33 + i] = batch->datum[i];
  fltbuf[36] = batch->phistt;
  fltbuf[37] = batch->phiend;
  for (i = 0; i < 3; ++i) 
    fltbuf[38 + i] = batch->scanax[i];
  fltbuf[41] = batch->time1;
  fltbuf[42] = batch->time2;
  fltbuf[43] = batch->bscale;
  fltbuf[44] = batch->bbfac;
  fltbuf[45] = batch->sdbscale;
  fltbuf[46] = batch->sdbfac;
  fltbuf[47] = batch->phirange;
  for (i = 0; i < 3; ++i)
    fltbuf[59 + i] = batch->e1[i];
  for (i = 0; i < 3; ++i)
    fltbuf[62 + i] = batch->e2[i];
  for (i = 0; i < 3; ++i)
    fltbuf[65 + i] = batch->e3[i];
  for (i = 0; i < 3; ++i)
    fltbuf[80 + i] = batch->source[i];
  for (i = 0; i < 3; ++i)
    fltbuf[83 + i] = batch->so[i];
  fltbuf[86] = batch->alambd;
  fltbuf[87] = batch->delamb;
  fltbuf[88] = batch->delcor;
  fltbuf[89] = batch->divhd;
  fltbuf[90] = batch->divvd;
  for (i = 0; i < batch->ndet; ++i)
  { fltbuf[111 + (i * 6)] = batch->dx[i];
    fltbuf[112 + (i * 6)] = batch->theta[i];
    fltbuf[113 + (i * 6)] = batch->detlm[i][0][0];
    fltbuf[114 + (i * 6)] = batch->detlm[i][0][1];
    fltbuf[115 + (i * 6)] = batch->detlm[i][1][0];
    fltbuf[116 + (i * 6)] = batch->detlm[i][1][1];}

  return 1;
}

int MtzWhdrLine(CCP4File *fileout, int nitems, char buffer[]) {

  /* write header record to fileout. Record is filled from
     nitems to MTZRECORDLENGTH by blanks.

     If a C-style null terminator is encountered before nitems
     are copied then the null terminator character is not
     copied and the string is padded with blanks from that
     point onwards */

 char hdrrec[MTZRECORDLENGTH];
 int i,j;

 for (i = 0; i < nitems; ++i) {
   /* Trap for C-style null character */
   if (buffer[i] == '\0') {
     break;
   }
   hdrrec[i] = buffer[i];
 }
 for (j = i; j < MTZRECORDLENGTH; ++j)
   hdrrec[j] = ' ';

 return (ccp4_file_writechar(fileout, (uint8 *) hdrrec,MTZRECORDLENGTH));

}

int MtzWrefl(CCP4File *fileout, int ncol, float *refldata) {

  if (!fileout)  {
    ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_NoFile),"MtzWrefl",NULL);
    return 0;
  }
  return (ccp4_file_write(fileout, (uint8 *) refldata, ncol));

}

MTZ *MtzMalloc(int nxtal, int nset[])

{ MTZ *mtz;
  int i,j,itime[3];
  float zerocell[6]={0.0};
  char dummy_xname[17];

  ccp4_utils_itime(itime);
  sprintf(dummy_xname,"NULL_xname%2.2d%2.2d%2.2d",itime[0],itime[1],itime[2]);
  dummy_xname[16]='\0';

  /* Allocate main header and symmetry */
  mtz = (MTZ *) ccp4_utils_malloc(sizeof(MTZ));
  if (mtz == NULL) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_AllocFail),"MtzMalloc",NULL);
    return NULL;
  }
  memset(mtz,'\0',sizeof(MTZ));

  mtz->nxtal=0;
  ccp4array_new_size(mtz->xtal,5);
  /* Allocate crystals and datasets. */
  if (nxtal == 0) {
    mtz->xtal[0] = NULL;
  } else {
    for (i = 0; i < nxtal; ++i) {
      /* This adds mtz->xtal[i] */
      if ( ! MtzAddXtal(mtz,dummy_xname,"NULL_pname",zerocell) ) return NULL;
      mtz->xtal[i]->nset = 0;
      for (j = 0; j < nset[i]; ++j) {
        /* This adds mtz->xtal[i]->set[j] */
        if ( ! MtzAddDataset(mtz,mtz->xtal[i],"NULL_dname",0.0) ) return NULL;
      }
    }
  }

  /* initialise main header */
  mtz->filein = NULL;
  mtz->fileout = NULL;
  mtz->title[0] = '\0';
  mtz->hist = NULL;
  mtz->histlines = 0;
  mtz->nxtal = nxtal;
  mtz->ncol_read = 0;
  mtz->nref = 0;
  mtz->nref_filein = 0;
  mtz->refs_in_memory = 1;
  mtz->n_orig_bat = 0;
  mtz->resmax_out = 0.0f;
  mtz->resmin_out = 999.0f;
  sprintf(mtz->mnf.amnf,"NAN");
  mtz->mtzsymm.spcgrp = 0;
  mtz->mtzsymm.spcgrpname[0] = '\0';
  mtz->mtzsymm.nsym = 0;
  mtz->mtzsymm.nsymp = 0;
  mtz->mtzsymm.symtyp = '\0';
  mtz->mtzsymm.pgname[0] = '\0';
  mtz->mtzsymm.spg_confidence = '\0';
  mtz->batch = NULL;
  for (i = 0; i < 5; ++i) {
    mtz->order[i] = NULL;
  }
  mtz->xml = NULL;
  mtz->unknown_headers = NULL;
  mtz->n_unknown_headers = 0;

  return(mtz);
}

int MtzFree(MTZ *mtz)

/* Frees the memory reserved for 'mtz' */

{ int i,j,k;

 /* Close attached mtz files */
 if (mtz->filein) {
    ccp4_file_close(mtz->filein);
    mtz->filein = NULL;
 }
 if (mtz->fileout) {
    ccp4_file_close(mtz->fileout);
    mtz->fileout = NULL;
 }

  /* Loop over crystals */
  for (i = 0; i < mtz->nxtal; ++i) {
  /* Loop over datasets for each crystal */
   for (j = 0; j < mtz->xtal[i]->nset; ++j) {
  /* Loop over columns for each dataset */
    for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) {
      MtzFreeCol(mtz->xtal[i]->set[j]->col[k]);
    }
    ccp4array_free(mtz->xtal[i]->set[j]->col);
    free((void *) mtz->xtal[i]->set[j]);
   }
   ccp4array_free(mtz->xtal[i]->set);
   free((void *) mtz->xtal[i]);
  }
  ccp4array_free(mtz->xtal);

  if (mtz->batch) {
    MtzFreeBatch(mtz->batch);
    mtz->batch = NULL;
  }

  if (mtz->hist != NULL) 
    MtzFreeHist(mtz->hist);

  if (mtz->xml != NULL)
    free(mtz->xml);

  if (mtz->unknown_headers != NULL)
    free(mtz->unknown_headers);

  free((void *) mtz);
  return 1;
}

MTZBAT *MtzMallocBatch()

/* Allocates memory for a single batch header */

{ MTZBAT *batch;

  batch = (MTZBAT *) ccp4_utils_malloc(sizeof(MTZBAT));
  if (batch == NULL) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_AllocFail),"MtzMallocBatch",NULL);
    return NULL;
  }
  memset(batch,'\0',sizeof(MTZBAT));
  batch->next = NULL;

  return(batch);
}

int MtzFreeBatch(MTZBAT *batch) 

/* Frees the memory reserved for 'batch' */

{
  if (batch != NULL) {
    MtzFreeBatch(batch->next);
    batch->next = NULL;
    free(batch);
  }
  return 1;
}

MTZCOL *MtzMallocCol(MTZ *mtz, int nref)

{ MTZCOL *col;

  col = (MTZCOL *) ccp4_utils_malloc(sizeof(MTZCOL));
  if (col == NULL) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_AllocFail),"MtzMallocCol",NULL);
    return NULL;
  }
  memset(col,'\0',sizeof(MTZCOL));

  col->ref = NULL;
  if (mtz->refs_in_memory) {
    ccp4array_new_size(col->ref,nref);
    if (col->ref == NULL) {
      ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_AllocFail),"MtzMallocCol",NULL);
      return NULL;
    }
  }

  return(col);
}

int MtzFreeCol(MTZCOL *col)

{ if (col->ref) ccp4array_free(col->ref);
  free((void *) col);
  return 1;
}

char *MtzCallocHist(int nhist)

/* Allocates memory for the mtz history with 'nhist' lines */

{ char *hist;

 hist = (char *) ccp4_utils_calloc(nhist, sizeof(char)*MTZRECORDLENGTH);
 return(hist);
}

int MtzFreeHist(char *hist)

/* Frees the memory reserved for 'hist' */

{ 
  free((void *) hist);
  return 1;
}

MTZXTAL *MtzAddXtal(MTZ *mtz, const char *xname, const char *pname,
              const float cell[6])
{
  /* add a new crystal to the mtz */
  int i,x;
  MTZXTAL *xtal;

  xtal = (MTZXTAL *) ccp4_utils_malloc( sizeof(MTZXTAL) );
  if (! xtal ) { 
    ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_AllocFail),"MtzAddXtal",NULL);
    return NULL;
  }
  memset(xtal,'\0',sizeof(MTZXTAL));

  /* fill out the data */
  strncpy( xtal->xname, xname, 64 );
  xtal->xname[64] = '\0';
  strncpy( xtal->pname, pname, 64 );
  xtal->pname[64] = '\0';
  xtal->resmin = 100.0;
  xtal->resmax = 0.0;
  for (i = 0; i < 6; i++) xtal->cell[i] = cell[i];
  /* make new xtalid */
  for (i = x = 0; x < mtz->nxtal; x++)
    if (mtz->xtal[x]->xtalid > i) i = mtz->xtal[x]->xtalid;
  xtal->xtalid = ++i;
  xtal->nset = 0;
  /* create initial array of 10 pointers to datasets */
  ccp4array_new_size(xtal->set,10);

  /* add pointer to mtz */
  if ( ++mtz->nxtal > ccp4array_size(mtz->xtal))
    ccp4array_resize(mtz->xtal, mtz->nxtal + 2);
  mtz->xtal[ mtz->nxtal - 1 ] = xtal;

  return xtal;
}

MTZSET *MtzAddDataset(MTZ *mtz, MTZXTAL *xtl, const char *dname,
                 const float wavelength)
{
  /* add a new dataset to the xtal */
  int i,x,s;
  MTZSET *set;

  set = (MTZSET *) ccp4_utils_malloc( sizeof(MTZSET) );
  if ( ! set ) { 
    ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_AllocFail),"MtzAddDataset",NULL);
    return NULL;
  }
  memset(set,'\0',sizeof(MTZSET));

  /* fill out the data */
  strncpy( set->dname, dname, 64 );
  set->dname[64] = '\0';
  set->wavelength = wavelength;

  /* New setid is one more than greatest current setid.
     It must be at least 1, unless it is the base dataset setid=0. */
  if (!strcmp(set->dname,"HKL_base")) {
    set->setid = 0;
  } else {
    i = 0;
    for (x = 0; x < mtz->nxtal; x++)
      for (s = 0; s < mtz->xtal[x]->nset; s++)
        if (mtz->xtal[x]->set[s]->setid > i) i = mtz->xtal[x]->set[s]->setid;
    set->setid = ++i;
  }

  set->ncol = 0;
  /* create initial array of 20 pointers to columns */
  ccp4array_new_size(set->col,20);

  /* add pointer to xtal */
  if ( ++xtl->nset > ccp4array_size(xtl->set))
    ccp4array_resize(xtl->set, xtl->nset + 4);
  xtl->set[ xtl->nset - 1 ] = set;

  return set;
}

MTZCOL *MtzAddColumn(MTZ *mtz, MTZSET *set, const char *label,
                const char *type)
{
  /* add a new column to the dataset */
  int i,nref;
  union float_uint_uchar uf;
  MTZCOL *col;

  if (set->ncol == MCOLUMNS) { 
    printf("MtzAddColumn: No more columns! \n");
    return NULL;
  }

  /* allocate some memory for first column */
  if (!mtz->refs_in_memory) {
    nref = 0;
  } else if (mtz->nref == 0) {
    nref = 2000;
  } else {
    nref = mtz->nref;
  }
  col = MtzMallocCol(mtz, nref);
  if (col == NULL) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CMTZ_ERRNO(CMTZERR_AllocFail),"MtzAddColumn",NULL);
    return NULL;
  }

  /* fill out the data */
  strncpy( col->label, label, 30 );
  col->label[30] = '\0';
  strncpy( col->type, type, 2);
  col->type[2] = '\0';
  col->active = 1;
  col->source = 0;
  col->min = 1.e06;
  col->max = -1.e06;
  col->colsource[0] = '\0';
  col->grpname[0] = '\0';
  col->grptype[0] = '\0';
  col->grpposn = -1;
  /* add pointer to set */
  if ( ++set->ncol > ccp4array_size(set->col))
    ccp4array_resize(set->col, set->ncol + 9);
  set->col[ set->ncol - 1 ] = col;
  /* initialise column to MNF */
  if (strncmp (mtz->mnf.amnf,"NAN",3) == 0) {
   uf = ccp4_nan();
  } else {
   uf.f = mtz->mnf.fmnf;
  }
  for (i=0; i < nref; i++) col->ref[i] = uf.f;

  return col;
}

int MtzToggleColumn(MTZCOL *col)
{
  /* Toggle active flag of column */
  if (col->active) {
    col->active = 0;
  } else {
    col->active = 1;
  }
  return col->active;
}

MTZSET *MtzColSet(const MTZ *mtz, const MTZCOL *col)
{
  int x,s,c;
  for (x=0; x < mtz->nxtal; x++)
    for (s=0; s < mtz->xtal[x]->nset; s++)
      for (c=0; c < mtz->xtal[x]->set[s]->ncol; c++)
      if (mtz->xtal[x]->set[s]->col[c] == col)
        return mtz->xtal[x]->set[s];
  printf ("MtzColSet: no such column. \n"); 
  return NULL;
}

MTZXTAL *MtzSetXtal(const MTZ *mtz, const MTZSET *set)
{
  int x,s;
  for (x=0; x < mtz->nxtal; x++)
    for (s=0; s < mtz->xtal[x]->nset; s++)
      if (mtz->xtal[x]->set[s] == set)
        return mtz->xtal[x];
  printf ("MtzSetXtal: no such dataset. \n"); 
  return NULL;
}

int MtzNxtal(const MTZ *mtz)
{
  return mtz->nxtal;
}

int MtzNumActiveXtal(const MTZ *mtz)
{
  int k,ixtal=0;

  for (k=0; k < mtz->nxtal; k++)
     if (MtzNumActiveSetsInXtal(mtz,mtz->xtal[k])) ++ixtal;
  return ixtal;
}

MTZXTAL **MtzXtals(MTZ *mtz)
{
  return mtz->xtal;
}

MTZXTAL *MtzIxtal(const MTZ *mtz, const int ixtal)
{
  return mtz->xtal[ixtal];
}

int MtzNsetsInXtal(const MTZXTAL *xtal)
{
  return xtal->nset;
}

int MtzNumActiveSetsInXtal(const MTZ *mtz, const MTZXTAL *xtal)
{
  int k,iset=0;

  for (k=0; k < xtal->nset; k++)
     if (MtzNumActiveColsInSet(xtal->set[k]) ||
         MtzNbatchesInSet(mtz, xtal->set[k])) ++iset;
  return iset;
}

MTZSET **MtzSetsInXtal(MTZXTAL *xtal)
{
  return xtal->set;
}

MTZSET *MtzIsetInXtal(const MTZXTAL *xtal, const int iset)
{
  return xtal->set[iset];
}

int MtzNcolsInSet(const MTZSET *set)
{
  return set->ncol;
}

int MtzNumActiveColsInSet(const MTZSET *set)
{
  int k,icol=0;

  for (k=0; k < set->ncol; k++)
     icol += set->col[k]->active;
  return icol;
}

int MtzNumSourceColsInSet(const MTZSET *set)
{
  int k,icol=0;

  for (k=0; k < set->ncol; k++)
    if (set->col[k]->source) ++icol;
  return icol;
}

int MtzNbatchesInSet(const MTZ *mtz, const MTZSET *set)
{
  int i,ibatch=0;
  MTZBAT *batch;

  batch = mtz->batch;

  /* if new batch headers have been written, lose the old ones */
  if (MtzNbat(mtz) > mtz->n_orig_bat) 
    for (i=0; i < mtz->n_orig_bat; ++i)
       batch = batch->next;

  while (batch) {
    if (batch->nbsetid == set->setid) ++ibatch;
    batch = batch->next;
  }

  return ibatch;
}

MTZCOL **MtzColsInSet(MTZSET *set)
{
  return set->col;
}

MTZCOL *MtzIcolInSet(const MTZSET *set, const int icol)
{
  return set->col[icol];
}

char *MtzColType(MTZCOL *col)
{
  return col->type;
}

int MtzNset(const MTZ *mtz)
{
  int x,iset=0;
  for (x=0; x < mtz->nxtal; x++)
    iset += MtzNsetsInXtal(mtz->xtal[x]);
  return iset;
}

int MtzNumActiveSet(const MTZ *mtz)
{
  int x,iset=0;
  for (x=0; x < mtz->nxtal; x++)
    iset += MtzNumActiveSetsInXtal(mtz,mtz->xtal[x]);
  return iset;
}
  
int MtzNcol(const MTZ *mtz)
{
  int x,s,icol=0;
  for (x=0; x < mtz->nxtal; x++)
    for (s=0; s < mtz->xtal[x]->nset; s++)
      icol += MtzNcolsInSet(mtz->xtal[x]->set[s]);
  return icol;
}
  
int MtzNumActiveCol(const MTZ *mtz)
{
  int x,s,icol=0;
  for (x=0; x < mtz->nxtal; x++)
    for (s=0; s < mtz->xtal[x]->nset; s++)
      icol += MtzNumActiveColsInSet(mtz->xtal[x]->set[s]);
  return icol;
}
  
int MtzNumSourceCol(const MTZ *mtz)
{
  int x,s,icol=0;
  for (x=0; x < mtz->nxtal; x++)
    for (s=0; s < mtz->xtal[x]->nset; s++)
      icol += MtzNumSourceColsInSet(mtz->xtal[x]->set[s]);
  return icol;
}
   
int MtzNref(const MTZ *mtz)
{
  /* get the number of reflections in the mtz */

  return mtz->nref;
}
 
int MtzNbat(const MTZ *mtz)
{
  /* get the number of batches in the mtz */
  int cnt=0;
  MTZBAT *batch;

  for (batch = mtz->batch ; batch != NULL ; batch = batch->next ) 
    ++cnt;

  return cnt;
}
  
char *MtzXtalPath(const MTZXTAL *xtal)
{
  /* Return the full path name of a crystal */
  char *path;
  size_t length;

  length = strlen(xtal->xname)+2;
  path = (char *) ccp4_utils_malloc(length*sizeof(char));
  strcpy( path, "/" );
  strcat( path, xtal->xname );
  path[length-1] = '\0';
  return ( path );
}

char *MtzSetPath(const MTZ *mtz, const MTZSET *set)
{
  /* Return the full path name of a dataset */
  char *path, *path1;
  size_t length;

  path1 = MtzXtalPath( MtzSetXtal( mtz, set ) );
  length = strlen(path1)+strlen(set->dname)+2;
  path = ccp4_utils_malloc(length*sizeof(char));
  strcpy( path, path1 );
  free (path1);  
  strcat( path, "/" );
  strcat( path, set->dname );
  path[length-1] = '\0';
  return ( path );
}

char *MtzColPath(const MTZ *mtz, const MTZCOL *col)
{
  /* Return the full path name of a column */
  char *path, *path1;
  size_t length;

  path1 = MtzSetPath( mtz, MtzColSet( mtz, col ) );
  length = strlen(path1)+strlen(col->label)+2;
  path = (char *) ccp4_utils_malloc(length*sizeof(char));
  strcpy( path, path1 );
  free (path1);  
  strcat( path, "/" );
  strcat( path, col->label );
  path[length-1] = '\0';
  return ( path );
}

int MtzRJustPath(char *path, const char *partial, const int njust)
{
  /* Complete a right-justified path by prefixing with wildcards */
  int i, j;
  /* count the slashes */
  for ( i = j = 0; i < strlen(partial); i++ ) if ( partial[i] == '/' ) j++;

  strcpy( path, "");
  if ( j++ < njust ) strcat( path, "/" );
  while ( j++ < njust ) strcat( path, "*/" );
  strcat( path, partial );

  return 1;
}

int MtzPathMatch(const char *path1, const char *path2)
{
  /* test for match between two paths, including wildcards */
  /* this version only handles wildcards at the end of name components */
  int p1 = 0, p2 = 0;
  while ( path1[p1] != '\0' && path2[p2] != '\0' ) {    /* search both paths */
    if ( path1[p1] != path2[p2] ) {
      if ( path1[p1] != '*' && path2[p2] != '*' )
      return FALSE;                       /* non-wild mismatch is terminal */
      while ( path1[p1] != '/' && path1[p1] != '\0' ) p1++;   /* skip compnt */
      while ( path2[p2] != '/' && path2[p2] != '\0' ) p2++;   /* skip compnt */
    } else {
      p1++; p2++;
    }
  }
  return (path1[p1] == path2[p2]);          /* true only if both paths ended */
}


MTZCOL *MtzColLookup(const MTZ *mtz, const char *label)
{
  /* Returns a pointer to the column of mtz with the given `label`, or NULL */
  int x,s,c;
  char *path1, path2[200];

  /* complete the right-justified path */
  MtzRJustPath( path2, label, 3 );
  /* now find the matching column */
  for (x=0; x < mtz->nxtal; x++)  /* not much point in optimising this */
    for (s=0; s < mtz->xtal[x]->nset; s++)
      for (c=0; c < mtz->xtal[x]->set[s]->ncol; c++) {
        path1 = MtzColPath(mtz, mtz->xtal[x]->set[s]->col[c]);
        if ( MtzPathMatch( path1, path2 ) ) {
          free (path1);
          return mtz->xtal[x]->set[s]->col[c];
	}
        free (path1);
      }
  return NULL;
}

MTZSET *MtzSetLookup(const MTZ *mtz, const char *label)
{
  int x,s;
  char *path1, path2[200];

  /* complete the right-justified path */
  MtzRJustPath( path2, label, 2 );
  /* now find the matching column */
  for (x=0; x < mtz->nxtal; x++)  
    for (s=0; s < mtz->xtal[x]->nset; s++) {
      path1 = MtzSetPath(mtz, mtz->xtal[x]->set[s]);
      if ( MtzPathMatch( path1, path2 ) ) {
        free (path1);
        return mtz->xtal[x]->set[s];
      }
    free (path1);
    }
  return NULL;
}

MTZXTAL *MtzXtalLookup(const MTZ *mtz, const char *label)
{
  /* Returns a pointer to the crystal of mtz with the given `label`, or NULL */
  int x;
  char *path1, path2[200];

  /* complete the right-justified path */
  MtzRJustPath( path2, label, 1 );
  /* now find the matching column */
  for (x=0; x < mtz->nxtal; x++) { 
    path1 = MtzXtalPath(mtz->xtal[x]);
    if ( MtzPathMatch( path1, path2 ) ) {
      free (path1);
      return mtz->xtal[x];
    }
    free (path1);
  }
  return NULL;
}

