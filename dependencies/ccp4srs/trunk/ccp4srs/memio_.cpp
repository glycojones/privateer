//  $Id: memio_.cpp $
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
//  **** Module  :  memio_  <implementation>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::MemIO  - buffered I/O
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#include <string.h>
#include <zlib.h>

#include "memio_.h"

namespace ccp4srs  {

  MemIO::MemIO()  {

    compression_level  = Z_NO_COMPRESSION;
    compression_buffer = NULL;
    compression_alloc  = 0;

    buffer        = NULL;
    buffer_length = 0;
    buffer_pos    = 0;
    alloc_length  = 0;
    chunk_size    = 10000;

  }

  MemIO::~MemIO()  {
    free();
  }

  void MemIO::setCompressionLevel ( int compressionLevel )  {
    switch (compressionLevel)  {
      case 0:  compression_level = Z_NO_COMPRESSION;      break;
      case 1:  compression_level = Z_BEST_SPEED;          break;
      case 2:  compression_level = Z_DEFAULT_COMPRESSION; break;
      default: compression_level = Z_BEST_COMPRESSION;
    }
  }


  int MemIO::read ( mmdb::cpstr fileName, mmdb::io::GZ_MODE gzMode )  {
  mmdb::io::File f;
  int            rc;

    buffer_pos = 0;
    rc = MemIO::Ok;

    f.assign ( fileName,false,true,gzMode );

    if (f.reset(true))  rc = read ( f );
                  else  rc = MemIO::CantOpenFile;

    if (rc!=MemIO::Ok)
      free();

    f.shut();

    return rc;

  }

  int MemIO::read ( mmdb::io::RFile f )  {
  int    rc,len;
  uLongf dsize;

    buffer_pos = 0;

    if (f.ReadInt(&buffer_length))  {

      if (buffer_length==0)
        rc = MemIO::EmptyFile;

      else if (buffer_length>0)  { // no compression

        if (buffer_length>alloc_length)  {
          mmdb::FreeVectorMemory ( buffer,0 );
          alloc_length = buffer_length;
          mmdb::GetVectorMemory ( buffer,alloc_length,0 );
        }
        if (f.ReadFile(buffer,buffer_length))
              rc = MemIO::Ok;
        else  rc = MemIO::ReadError;

      } else  { // decompress

        buffer_length = -buffer_length - 1;
        if (buffer_length>alloc_length)  {
          mmdb::FreeVectorMemory ( buffer,0 );
          alloc_length = buffer_length;
          mmdb::GetVectorMemory ( buffer,alloc_length,0 );
        }

        if (f.ReadInt(&len))  {

          if (compression_alloc<=len)  {
            mmdb::FreeVectorMemory ( compression_buffer,0 );
            compression_alloc = len + 100;
            mmdb::GetVectorMemory  ( compression_buffer,compression_alloc,0 );
          }

          if (f.ReadFile(compression_buffer,len))  {
            dsize = alloc_length;
            rc = uncompress ( buffer,&dsize,compression_buffer,len );
            if (rc==Z_OK)  {
              if (int(dsize)!=buffer_length)  rc = MemIO::SizeError;
                                        else  rc = MemIO::Ok;
            } else
              rc = MemIO::UncompressionError;
          } else
            rc = MemIO::ReadError;

        } else
          rc = MemIO::ReadError;

      }

    } else
      rc = MemIO::ReadError;

    return rc;

  }

  int MemIO::write ( mmdb::cpstr fileName, mmdb::io::GZ_MODE gzMode )  {
  mmdb::io::File f;
  int            rc;

    rc = MemIO::Ok;

    f.assign ( fileName,false,true,gzMode );

    if (f.rewrite()) rc = write(f);
                else rc = MemIO::CantOpenFile;

    f.shut();

    return rc;

  }

  int MemIO::write ( mmdb::io::RFile f )  {
  int    rc,len;
  uLongf csize;

    rc = MemIO::Ok;

    if ((compression_level==Z_NO_COMPRESSION) || (buffer_length<=0))  {

      if (f.WriteInt(&buffer_length))  {
        if (!f.WriteFile(buffer,buffer_length))
          rc = MemIO::WriteError;
      }

    } else  {

      len = -buffer_length - 1;
      if (f.WriteInt(&len))  {

        csize = (uLongf)(1.1*buffer_length);
        if (compression_alloc<=int(csize))  {
          mmdb::FreeVectorMemory ( compression_buffer,0 );
          mmdb::GetVectorMemory  ( compression_buffer,csize,0 );
        }

        rc = compress2 ( compression_buffer,&csize,
                         buffer,buffer_length,compression_level );

        if (rc==Z_OK)  {
          len = csize;
          if (f.WriteInt(&len))  {
            if (!f.WriteFile(compression_buffer,csize))
                  rc = MemIO::WriteError;
            else  rc = MemIO::Ok;
          } else
            rc = MemIO::WriteError;
        } else
          rc = MemIO::CompressionError;

      } else
        rc = MemIO::WriteError;

    }

    return rc;

  }

  void MemIO::reset()  {
    buffer_length = 0;
  }

  void MemIO::free()  {
    mmdb::FreeVectorMemory ( compression_buffer,0 );
    mmdb::FreeVectorMemory ( buffer            ,0 );
    compression_alloc = 0;
    buffer_length     = 0;
    buffer_pos        = 0;
    alloc_length      = 0;
  }

  void *MemIO::get_buffer ( const int length )  {
    if (alloc_length<length)  {
      mmdb::FreeVectorMemory ( buffer,0 );
      alloc_length = length;
      mmdb::GetVectorMemory ( buffer,alloc_length,0 );
    }
    return buffer;
  }

  void MemIO::get_buffer ( mmdb::pstr * buf, int * length )  {
    *buf    = mmdb::pstr(buffer);
    *length = buffer_length;
  }


  void MemIO::write_buffer ( const void * src, const int length  )  {
  mmdb::bvector b1;
    if (buffer_length+length>=alloc_length)  {
      alloc_length += length+chunk_size;
      mmdb::GetVectorMemory ( b1,alloc_length,0 );
      memcpy ( b1,buffer,buffer_length );
      mmdb::FreeVectorMemory ( buffer,0 );
      buffer = b1;
    }
    memcpy ( &(buffer[buffer_length]),src,length );
    buffer_length += length;
  }


  void MemIO::put_integer ( int I )  {
  mmdb::intUniBin iUB;
    mmdb::int2UniBin ( I,iUB );
    write_buffer ( iUB,sizeof(mmdb::intUniBin) );
  }

  void MemIO::put_ishort ( int I )  {
  mmdb::shortUniBin sUB;
    mmdb::short2UniBin ( short(I),sUB );
    write_buffer ( sUB,sizeof(mmdb::shortUniBin) );
  }

  void MemIO::put_ibyte ( int I )  {
  // here we do not rely only on little-endian architectures
  mmdb::byte B;
    B = mmdb::byte(I);
    write_buffer ( &B,sizeof(mmdb::byte) );
  }

  void MemIO::put_byte ( mmdb::byte B )  {
    write_buffer ( &B,sizeof(mmdb::byte) );
  }

  void MemIO::put_short ( short I )  {
  mmdb::shortUniBin sUB;
    mmdb::short2UniBin ( I,sUB );
    write_buffer ( sUB,sizeof(mmdb::shortUniBin) );
  }

  void MemIO::put_long ( long I )  {
  mmdb::longUniBin lUB;
    mmdb::long2UniBin ( I,lUB );
    write_buffer ( lUB,sizeof(mmdb::longUniBin) );
  }

  void MemIO::put_word ( mmdb::word W )  {
  mmdb::wordUniBin wUB;
    mmdb::word2UniBin ( W,wUB );
    write_buffer ( wUB,sizeof(mmdb::wordUniBin) );
  }

  void MemIO::put_real ( mmdb::realtype R )  {
  mmdb::realUniBin rUB;
    mmdb::real2UniBin ( R,rUB );
    write_buffer ( rUB,sizeof(mmdb::realUniBin) );
  }

  void MemIO::put_float ( mmdb::realtype R )  {
  mmdb::floatUniBin fUB;
    mmdb::float2UniBin ( R,fUB );
    write_buffer ( fUB,sizeof(mmdb::floatUniBin) );
  }

  void MemIO::put_shortreal ( mmdb::shortreal R )  {
  mmdb::shortrealUniBin srUB;
    mmdb::shortreal2UniBin ( R,srUB );
    write_buffer ( srUB,sizeof(mmdb::shortrealUniBin) );
  }

  void MemIO::put_line ( mmdb::cpstr L )  {
    write_buffer ( L,strlen(L)+1 );  // null-terminated
  }

  void MemIO::put_string ( mmdb::cpstr S )  {
  int len;
    if (S)  len = strlen ( S );
      else  len = 0;
    put_integer ( len );
    if (len>0)
      write_buffer ( S,len );  // not null-terminated
  }

  void MemIO::put_bool ( bool B )  {
  char c;
    if (B)  c = 'Y';
      else  c = 'N';
    write_buffer ( &c,1 );
  }


  bool MemIO::read_buffer ( void * dest, const int length,
                                      bool * Ok )  {
    if (buffer_pos+length<=buffer_length)  {
      memcpy ( dest,&(buffer[buffer_pos]),length );
      buffer_pos += length;
      return true;
    }
    if (Ok) *Ok = false;
    return false;
  }

  bool MemIO::get_integer ( int & I, bool * Ok )  {
  mmdb::intUniBin iUB;
    if (read_buffer(iUB,sizeof(mmdb::intUniBin),Ok))  {
      mmdb::UniBin2int ( iUB,I );
      return true;
    } else  {
      I = 0;
      return false;
    }
  }

  bool MemIO::get_ishort ( int & I, bool * Ok )  {
  mmdb::shortUniBin sUB;
  short             Ishort;
    if (read_buffer(sUB,sizeof(mmdb::shortUniBin),Ok))  {
      mmdb::UniBin2short ( sUB,Ishort );
      I = Ishort;
      return true;
    } else  {
      I = 0;
      return false;
    }
  }

  bool MemIO::get_ibyte ( int & I, bool * Ok )  {
  mmdb::byte B;
    if (read_buffer(&B,sizeof(mmdb::byte),Ok))  {
      I = B;
      return true;
    } else  {
      I = 0;
      return false;
    }
  }

  bool MemIO::get_byte ( mmdb::byte & B, bool * Ok )  {
    return read_buffer ( &B,sizeof(mmdb::byte),Ok );
  }

  bool MemIO::get_short ( short & I, bool * Ok )  {
  mmdb::shortUniBin sUB;
    if (read_buffer(sUB,sizeof(mmdb::shortUniBin),Ok))  {
      mmdb::UniBin2short ( sUB,I );
      return true;
    } else  {
      I = 0;
      return false;
    }
  }

  bool MemIO::get_long ( long & I, bool * Ok )  {
  mmdb::longUniBin lUB;
    if (read_buffer(lUB,sizeof(mmdb::longUniBin),Ok))  {
      mmdb::UniBin2long ( lUB,I );
      return true;
    } else  {
      I = 0;
      return false;
    }
  }

  bool MemIO::get_word ( mmdb::word & W, bool * Ok )  {
  mmdb::wordUniBin wUB;
    if (read_buffer(wUB,sizeof(mmdb::wordUniBin),Ok))  {
      mmdb::UniBin2word ( wUB,W );
      return true;
    } else  {
      W = 0;
      return false;
    }
  }

  bool MemIO::get_real ( mmdb::realtype & R, bool * Ok )  {
  mmdb::realUniBin rUB;
    if (read_buffer(rUB,sizeof(mmdb::realUniBin),Ok))  {
      mmdb::UniBin2real ( rUB,R );
      return true;
    } else  {
      R = 0.0;
      return false;
    }
  }

  bool MemIO::get_float ( mmdb::realtype & R, bool * Ok )  {
  mmdb::floatUniBin fUB;
    if (read_buffer(fUB,sizeof(mmdb::floatUniBin),Ok))  {
      mmdb::UniBin2float ( fUB,R );
      return true;
    } else  {
      R = 0.0;
      return false;
    }
  }

  bool MemIO::get_shortreal ( mmdb::shortreal & R, bool * Ok )  {
  mmdb::shortrealUniBin srUB;
    if (read_buffer(srUB,sizeof(mmdb::shortrealUniBin),Ok))  {
      mmdb::UniBin2shortreal ( srUB,R );
      return true;
    } else  {
      R = 0.0;
      return false;
    }
  }

  bool MemIO::get_line ( mmdb::pstr L, bool * Ok )  {
  int  i;
  char c;
    i = 0;
    c = ' ';
    while (c && (buffer_pos<buffer_length))  {
      c      = buffer[buffer_pos++];
      L[i++] = c;
    }
    if (!c)
      return true;
    else  {
      if (Ok)  *Ok = false;
      return false;
    }
  }

  bool MemIO::get_string ( mmdb::pstr & L, bool * Ok )  {
  int len;
    if (L)  delete[] L;
    L = NULL;
    if (get_integer(len,Ok))  {
      if (len>0)  {
        L = new char[len+1];
        L[len] = char(0);
        return read_buffer ( L,len,Ok );
      }
      return true;
    }
    return false;
  }

  bool MemIO::get_bool ( bool & B, bool * Ok )  {
  char c;
    if (read_buffer(&c,1,Ok))  {
      B = (c=='Y');
      return true;
    }
    return false;
  }

}  // namespace ccp4srs
