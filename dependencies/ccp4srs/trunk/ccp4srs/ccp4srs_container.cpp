//  $Id: ccp4srs_container.h, v 1.00 2010/03/11 15:28:14 Eugene Exp $
//  =================================================================
//
//   CCP4 SRS Library: Storage, Retrieval and Search support for
//   CCP4 ligand data.
//
//   Copyright (C) Eugene Krissinel 2010.
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
//    29.04.10   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  ccp4srs_container  <implementation>
//       ~~~~~~~~~
//  **** Classes :  CCP4SRSContainer  - template container class
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010
//
//  =================================================================
//

/*
#include "ccp4srs_container.h"

template <class T>
CCP4SRSContainer<T>::CCP4SRSContainer()  {
  object    = NULL;
  n_objects = 0;
  n_alloc   = 0;
}

template <class T>
CCP4SRSContainer<T>::~CCP4SRSContainer()  {
  empty();
}

template <class T>
void CCP4SRSContainer<T>::empty()  {
int i;
  for (i=0;i<n_alloc;i++)
    if (object[i])  delete object[i];
  delete[] object;
  object    = NULL;
  n_objects = 0;
  n_alloc   = 0;
}

template <class T>
void CCP4SRSContainer<T>::add ( T* obj )  {
T**  tobject;
int  i;
  if (n_objects>=n_alloc)  {
    n_alloc = n_objects+10;
    tobject = new T*[n_alloc];
    for (i=0;i<n_objects;i++)
      tobject[i] = object[i];
    for (i=n_objects;i<n_alloc;i++)
      tobject[i] = NULL;
    delete[] object;
    object = tobject;
  }
  obj[n_objects++] = obj;
}

template <class T>
void  CCP4SRSContainer<T>::write_mem ( PCMemIO memIO, int version )  {
int i;
  memIO->put_integer ( n_objects );
  for (i=0;i<n_objects;i++)
    object[i]->write_mem ( memIO,version );
}

template <class T>
Boolean CCP4SRSContainer<T>::read_mem ( PCMemIO memIO, int version,
                                        Boolean * Ok )  {
int     i;
Boolean success;
  empty();
  if (Ok)  success = *Ok;
     else  success = True;
  memIO->get_integer ( n_objects,&success );
  if (n_objects>0)  {
    object = new T*[n_alloc];
    for (i=0;i<n_objects;i++)  {
      object[i] = new T();
      object[i]->read_mem ( memIO,version,success );
    }
  }
  if (Ok)  *Ok = success;
  return success;
}

*/
