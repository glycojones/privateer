/*
     ccp4_array.c: implementation file for resizable array implementation.
     Copyright (C) 2002  Kevin Cowtan

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

/** @file ccp4_array.c
 *  implementation file for resizable array implementation.
 *  Kevin Cowtan
 */

#include "ccp4_array.h"
/* rcsid[] = "$Id$" */

ccp4_ptr ccp4array_new_(ccp4_ptr *p)
{
  ccp4array_base *v;
  v = (ccp4array_base *)malloc(sizeof(ccp4array_base));
  v->size = v->capacity = 0;
  *p = (ccp4_ptr *)((ccp4_byteptr)(v)+sizeof(ccp4array_base));
  return *p;
}

ccp4_ptr ccp4array_new_size_(ccp4_ptr *p, const int size, const size_t reclen)
{
  ccp4array_base *v;
  int capacity = (size * 12) / 10 + 2;
  v = (ccp4array_base *)malloc(sizeof(ccp4array_base) + capacity * reclen);
  v->size = size;
  v->capacity = capacity;
  *p = (ccp4_ptr *)((ccp4_byteptr)(v)+sizeof(ccp4array_base));
  return *p;
}

void ccp4array_resize_(ccp4_ptr *p, const int size, const size_t reclen)
{
  ccp4array_base *v;
  v = (ccp4array_base *)((ccp4_byteptr)(*p)-sizeof(ccp4array_base));
  if (size > v->capacity) {
    v->capacity = (size * 12) / 10 + 2;
    v = (ccp4array_base *)realloc(v, sizeof(ccp4array_base) + v->capacity * reclen);
    *p = (ccp4_ptr *)((ccp4_byteptr)(v)+sizeof(ccp4array_base));
  }
  v->size = size;
}

void ccp4array_reserve_(ccp4_ptr *p, const int size, const size_t reclen)
{
  ccp4array_base *v;
  v = (ccp4array_base *)((ccp4_byteptr)(*p)-sizeof(ccp4array_base));
  v->capacity = size;
  if ( v->size > size ) v->size = size;
  v = (ccp4array_base *)realloc(v, sizeof(ccp4array_base) + v->capacity * reclen);
  *p = (ccp4_ptr *)((ccp4_byteptr)(v)+sizeof(ccp4array_base));
}

void ccp4array_append_(ccp4_ptr *p, ccp4_constptr data, const size_t reclen)
{
  ccp4array_base *v;
  int osize;
  v = (ccp4array_base *)((ccp4_byteptr)(*p)-sizeof(ccp4array_base));
  osize = v->size;
  ccp4array_resize_(p, osize+1, reclen);
  memcpy((ccp4_byteptr)(*p)+osize*reclen, data, reclen);
}

void ccp4array_append_n_(ccp4_ptr *p, ccp4_constptr data, const int n, const size_t reclen)
{
  ccp4array_base *v;
  ccp4_byteptr newdata;
  int osize, i;
  v = (ccp4array_base *)((ccp4_byteptr)(*p)-sizeof(ccp4array_base));
  osize = v->size;
  ccp4array_resize_(p, osize+n, reclen);
  newdata = (ccp4_byteptr)(*p)+osize*reclen;
  for ( i = 0; i < n; i++ ) {
    memcpy(newdata, data, reclen);
    newdata += reclen;
  }
}

void ccp4array_append_list_(ccp4_ptr *p, ccp4_constptr data, const int n, const size_t reclen)
{
  ccp4array_base *v;
  int osize;
  v = (ccp4array_base *)((ccp4_byteptr)(*p)-sizeof(ccp4array_base));
  osize = v->size;
  ccp4array_resize_(p, osize+n, reclen);
  memcpy((ccp4_byteptr)(*p)+osize*reclen, data, n * reclen);
}

void ccp4array_insert_(ccp4_ptr *p, const int i, ccp4_constptr data, const size_t reclen)
{
  ccp4array_base *v;
  int osize;
  v = (ccp4array_base *)((ccp4_byteptr)(*p)-sizeof(ccp4array_base));
  osize = v->size;
  ccp4array_resize_(p, osize+1, reclen);
  memmove((ccp4_byteptr)(*p)+(i+1)*reclen, (ccp4_byteptr)(*p)+i*reclen, (osize-i)*reclen);
  memcpy((ccp4_byteptr)(*p)+i*reclen, data, reclen);  
}

void ccp4array_delete_ordered_(ccp4_ptr *p, const int i, const size_t reclen)
{
  ccp4array_base *v;
  int nsize;
  v = (ccp4array_base *)((ccp4_byteptr)(*p)-sizeof(ccp4array_base));
  nsize = v->size - 1;
  memmove((ccp4_byteptr)(*p)+i*reclen, (ccp4_byteptr)(*p)+(i+1)*reclen, (nsize-i)*reclen);
  v->size--; /* ccp4array_resize_(p, nsize, reclen); */
}

void ccp4array_delete_(ccp4_ptr *p, const int i, const size_t reclen)
{
  ccp4array_base *v;
  int nsize;
  v = (ccp4array_base *)((ccp4_byteptr)(*p)-sizeof(ccp4array_base));
  nsize = v->size - 1;
  memcpy((ccp4_byteptr)(*p)+i*reclen, (ccp4_byteptr)(*p)+nsize*reclen, reclen);
  v->size--; /* ccp4array_resize_(p, size, reclen); */
}

void ccp4array_delete_last_(ccp4_ptr *p, const size_t reclen)
{
  ccp4array_base *v;
  v = (ccp4array_base *)((ccp4_byteptr)(*p)-sizeof(ccp4array_base));
  v->size--; /* ccp4array_resize_(p, v->size-1, reclen); */
}

int ccp4array_size_(ccp4_constptr *p)
{
  ccp4array_base *v;
  v = (ccp4array_base *)((ccp4_byteptr)(*p)-sizeof(ccp4array_base));
  return v->size;
}

void ccp4array_free_(ccp4_ptr *p)
{
  ccp4array_base *v;
  v = (ccp4array_base *)((ccp4_byteptr)(*p)-sizeof(ccp4array_base));
  free(v);
}
