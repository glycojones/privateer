//  $Id: mpfile_.cpp $
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
//  **** Module  :  mpfile_  <implementation>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::MPFile  - multi-part file support
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#include "mpfile_.h"

namespace ccp4srs  {

  //  ==========================  MPFile  ============================

  void MPFilePos::Copy ( RMPFilePos dbf_pos )  {
    chunkNo = dbf_pos.chunkNo;
    offset  = dbf_pos.offset;
  }

  void MPFilePos::write ( mmdb::io::RFile f )  {
    f.WriteByte ( &chunkNo );
    f.WriteLong ( &offset  );
  }

  void MPFilePos::read ( mmdb::io::RFile f )  {
    f.ReadByte ( &chunkNo );
    f.ReadLong ( &offset  );
  }


  MPFile::MPFile()  {
    f         = NULL;
    fpath     = NULL;
    chunkPath = NULL;
    uniBinary = true;
    readOnly  = false;
    chunkSize = 1000000000;
    chunkNo   = 0;
    nchunks   = 0;
  }

  MPFile::~MPFile()  {
    shut();
    if (fpath)  {
      delete[] fpath;
      fpath = NULL;
    }
  }

  void MPFile::shut()  {
  int i;
    if (f)  {
      for (i=0;i<nchunks;i++)
        if (f[i])  f[i]->shut();
      delete[] f;
      f = NULL;
    }
    if (chunkPath)  {
      delete[] chunkPath;
      chunkPath = NULL;
    }
    uniBinary = true;
    readOnly    = false;
    chunkSize = 1000000000;
    chunkNo   = 0;
    nchunks   = 0;
  }

  void MPFile::assign ( mmdb::cpstr FName, bool UniB )  {
    shut();
    mmdb::CreateCopy ( fpath,FName );
    uniBinary = UniB;
  }

  mmdb::pstr MPFile::getChunkPath()  {
  char N[50];
    sprintf ( N,".fp%i",chunkNo+1 );
    mmdb::CreateCopCat ( chunkPath,fpath,N );
    return chunkPath;
  }

  bool MPFile::reset ( bool rdOnly )  {
    shut();
    readOnly = rdOnly;
    chunkNo  = 0;
    nchunks  = 1;
    f    = new mmdb::io::PFile[nchunks];
    f[0] = new mmdb::io::File();
    f[0]->assign ( getChunkPath(),false,uniBinary,mmdb::io::GZM_NONE );
    if (f[0]->reset(readOnly))
         return f[0]->ReadLong ( &chunkSize );
    else return false;
  }

  bool MPFile::rewrite ( long chunk_size )  {

    shut();

    readOnly  = false;
    chunkSize = chunk_size;
    chunkNo   = 1;
    nchunks   = 1;
    f    = new mmdb::io::PFile[nchunks];
    f[0] = new mmdb::io::File();
    while (chunkNo>0)  {
      f[0]->assign ( getChunkPath(),false,uniBinary,mmdb::io::GZM_NONE );
      if (f[0]->exists())  {
        f[0]->erase();
        chunkNo++;
      } else
        chunkNo = 0;
    }
    f[0]->assign ( getChunkPath(),false,uniBinary,mmdb::io::GZM_NONE );
    if (f[0]->rewrite())  return f[0]->WriteLong ( &chunkSize );
                    else  return false;

  }

  bool MPFile::append ( long chunk_size )  {
  mmdb::io::PFile fa;
  int    i;

    shut();

    readOnly = false;
    chunkNo  = 0;
    nchunks  = 1;

    fa = new mmdb::io::File();
    fa->assign ( getChunkPath(),false,uniBinary,mmdb::io::GZM_NONE );

    if (!fa->exists())  {

      f = new mmdb::io::PFile[nchunks];
      chunkSize = chunk_size;
      f[0] = fa;
      if (f[0]->rewrite())
           return f[0]->WriteLong ( &chunkSize );
      else return false;

    } else if (fa->reset(readOnly))  {

      if (!fa->ReadLong(&chunkSize))  {
        delete fa;
        return false;
      }

    } else
      return false;

    while (fa->exists())  {
      chunkNo++;
      fa->assign ( getChunkPath(),false,uniBinary,mmdb::io::GZM_NONE );
    }

    nchunks = chunkNo;
    chunkNo--;

    f = new mmdb::io::PFile[nchunks];
    for (i=0;i<nchunks;i++)
      f[i] = NULL;

    f[chunkNo] = fa;
    f[chunkNo]->assign ( getChunkPath(),false,uniBinary,mmdb::io::GZM_NONE );
    return f[chunkNo]->append();

  }

  void MPFile::getPosition ( RMPFilePos DBFPos )  {
    if (f)  {
      if (f[chunkNo])  {
        DBFPos.chunkNo = chunkNo;
        DBFPos.offset  = f[chunkNo]->Position();
        return;
      }
    }
    DBFPos.chunkNo = 0;
    DBFPos.offset  = 0;
  }

  bool MPFile::getChunk ( int toChunkNo, bool onRead )  {
  mmdb::io::PPFile f1;
  int              i;
  bool             B;
    if (toChunkNo>=nchunks)  {
      f1 = new mmdb::io::PFile[toChunkNo+1];
      for (i=0;i<nchunks;i++)
        f1[i] = f[i];
      for (i=nchunks;i<=toChunkNo;i++)
        f1[i] = NULL;
      if (f)  delete[] f;
      f = f1;
      nchunks = toChunkNo+1;
    }
    chunkNo = toChunkNo;
    if (!f[chunkNo])  {
      f[chunkNo] = new mmdb::io::File();
      f[chunkNo]->assign ( getChunkPath(),false,uniBinary,
                           mmdb::io::GZM_NONE );
    }
    if (!f[chunkNo]->isOpen())  {
      if (onRead)  {
        B = f[chunkNo]->reset ( readOnly );
        if (B && (chunkNo==0))
          f[chunkNo]->ReadLong ( &chunkSize );
      } else  {
        B = f[chunkNo]->rewrite();
        if (B && (chunkNo==0))
          f[chunkNo]->WriteLong ( &chunkSize );
      }
      if (!B)  {
        delete f[chunkNo];
        f[chunkNo] = NULL;
      }
      return B;
    }
    return true;
  }

  bool MPFile::seek ( RMPFilePos DBFPos )  {
    if (getChunk(DBFPos.chunkNo,true))
         return f[chunkNo]->seek ( DBFPos.offset );
    else return false;
  }

  bool MPFile::MPFileEnd()  {
    if (f)  {
      if (f[chunkNo])
        return f[chunkNo]->FileEnd();
    }
    return true;
  }

  bool MPFile::Success()  {
    if (f)  {
      if (f[chunkNo])
        return f[chunkNo]->Success();
    }
    return false;
  }

  bool MPFile::checkChunkEnd ( bool onRead )  {
    if (getChunk(chunkNo,onRead))  {
      if (f[chunkNo]->Position()>=chunkSize)
        return getChunk ( chunkNo+1,onRead );
      return true;
    }
    return false;
  }


  mmdb::io::PFile MPFile::getChunkRead()  {
    if (checkChunkEnd(true))  return f[chunkNo];
                        else  return NULL;
  }

  mmdb::io::PFile MPFile::getChunkWrite()  {
    if (checkChunkEnd(false))  return f[chunkNo];
                         else  return NULL;
  }

}  // namespace ccp4srs
