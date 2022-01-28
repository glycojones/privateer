//  $Id: ccp4srs_container.h $
//  =================================================================
//
//   CCP4 SRS Library: Storage, Retrieval and Search support for
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
//  **** Module  :  ccp4srs_container  <interface>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Container  - template container class
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#ifndef CCP4SRS_CONTAINER_H
#define CCP4SRS_CONTAINER_H

#include "memio_.h"

namespace ccp4srs  {

  template <class T>

  class Container  {

    protected:
      int  n_objects;
      T  **object;

    private:
      int  n_alloc;

    public:

      Container()  {
        object    = NULL;
        n_objects = 0;
        n_alloc   = 0;
      }

      virtual ~Container()  {
        empty();
      }

      void empty()  {
      int i;
        if (object)  {
          for (i=0;i<n_alloc;i++)
            if (object[i])  delete object[i];
          delete[] object;
          object = NULL;
        }
        n_objects = 0;
        n_alloc   = 0;
      }

      inline int numberOf()      { return n_objects;   }
      inline T*  at ( int pos )  { return object[pos]; }

      int  index ( mmdb::cpstr id )  {
      int  i,k;
        k = -1;
        if (id)  {
          for (i=0;(i<n_objects) && (k<0);i++)
            if (!strcmp(id,object[i]->id()))  k = i;
        }
        return k;
      }

      void add ( T* obj )  {
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
        object[n_objects++] = obj;
      }

      void  write_mem ( PMemIO memIO, int version )  {
      int i;
        memIO->put_integer ( n_objects );
        for (i=0;i<n_objects;i++)
          object[i]->write_mem ( memIO,version );
      }

      bool read_mem ( PMemIO memIO, int version, bool * Ok )  {
      int  i;
      bool success;
        empty();
        if (Ok)  success = *Ok;
           else  success = true;
        memIO->get_integer ( n_objects,&success );
        n_alloc = n_objects;
        if (n_objects>0)  {
          object = new T*[n_alloc];
          for (i=0;i<n_objects;i++)  {
            object[i] = new T();
            object[i]->read_mem ( memIO,version,&success );
          }
        }
        if (Ok)  *Ok = success;
        return success;
      }

  };

}  // namespace ccp4srs

#endif // CCP4SRS_CONTAINER_H
