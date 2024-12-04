/*
     ccp4_array.h: header file for resizable array implementation. 
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

/** @file ccp4_array.h
 *  header file for resizable array implementation.
 *  Kevin Cowtan
 */

/*
CCP4 resizable array implementation.

This defines an object and methods which looks just like a simple C
array, but can be resized at will without incurring excessive
overheads.

A pointer to the desired type is created. Array elements are accessed
from this pointer as normal. Other operations depend on macros which
extract the stored type from the type of the pointer.

The array is managed with a header, which is positioned before the
beginning of the array. The malloc'ed memory area starts at the
beginning of this header. However the pointer to the header is not
stored, but rather derived from the array pointer whenever it is
required.

Arrays have a size and a capacity. The size is the number of elements
in use, and the capacity is the number of elements available before a
new memory allocation is required. When new memory is required, an
excess is reqested to allow the array to grow further before
performing another malloc.

If the precise amount of memory is known, the capacity can be
controlled directly using the 'reserve' macro.

Example: to handle an array of type mytype:

\code
  int i;
  mytype x,y;
  mytype *array;

  ccp4array_new(array);

  ccp4array_append_n(array, x, 3);

  for ( i = 0; i < 3; i++ ) y = array[i];

  ccp4array_free(array);
\endcode
*/

#ifndef __CCP4_ARRAY_INC
#define __CCP4_ARRAY_INC

#ifdef  __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <string.h>
/* rcsidha[] = "$Id$" */

/*! constant pointer type */
typedef const void *ccp4_constptr;
/*! byte pointer type */
typedef char *ccp4_byteptr;
/*! pointer type */
typedef void *ccp4_ptr;

/*! struct definition for the array pre-header */
typedef struct ccp4array_base_ {
  int size, capacity;
} ccp4array_base;

/*! Macro to allocate a new array.
  The array is allocated with a size and capacity of 0
  \param v   The array pointer
  \return    The new array pointer (redundent)
*/
#define ccp4array_new(v) ccp4array_new_((ccp4_ptr*)(&v))

/*! Macro to allocate a new array with non-zero size.
  The array is allocated with a size of s and capacity of at least s
  \param v   The array pointer
  \param s   The new size
  \return    The new array pointer (redundent)
*/
#define ccp4array_new_size(v,s) ccp4array_new_size_((ccp4_ptr*)(&v),s,sizeof(*v))

/*! Macro to resize an array.
  This changes the size. Memory allocation only takes place if the new
  size is greater than the capacity. If that occurs, the new capacity
  will be slightly greater than the requested size, to allow room for
  expansion.
  \param v   The array pointer
  \param s   The new size
*/
#define ccp4array_resize(v,s) ccp4array_resize_((ccp4_ptr*)(&v),s,sizeof(*v))

/*! Macro to reserve space for an array.
  This forces a memory reallocation. The size of the array is
  unchanged, unless the new capacity is less than the current size, in
  which case the size is set to the new capacity. Unlike resize, the
  new allocation will be exactly the size of the array.
  \param v   The array pointer
  \param s   The new capacity
*/
#define ccp4array_reserve(v,s) ccp4array_reserve_((ccp4_ptr*)(&v),s,sizeof(*v))

/*! Macro to append an element to an array.
  This increments the size. Memory allocation only takes place if the new
  size is greater than the capacity.
  \param v   The array pointer
  \param d   The new element (may not be a literal)
*/
#define ccp4array_append(v,d) ccp4array_append_((ccp4_ptr*)(&v),(ccp4_constptr)(&d),sizeof(*v))

/*! Macro to append n copies of an element to an array.
  This increments the size by n. Memory allocation only takes place if the new
  size is greater than the capacity.
  \param v   The array pointer
  \param d   The new element (may not be a literal)
  \param n   The number of copies to append
*/
#define ccp4array_append_n(v,d,n) ccp4array_append_n_((ccp4_ptr*)(&v),(ccp4_constptr)(&d),n,sizeof(*v))

/*! Macro to append n elements from another list to an array.
  This increment the size by n. Memory allocation only takes place if the new
  size is greater than the capacity.
  \param v   The array pointer
  \param l   Pointer to the list
  \param n   The number of copies to append
*/
#define ccp4array_append_list(v,l,n) ccp4array_append_list_((ccp4_ptr*)(&v),(ccp4_constptr)l,n,sizeof(*v))

/*! Macro to insert an element before the element[i] of an array.
  This increments the size. All subsequent elements are moved up. As a
  result this method is slow.
  \param v   The array pointer
  \param d   The new element (may not be a literal)
  \param i   The element before which the insertion is to be made.
*/
#define ccp4array_insert(v,i,d) ccp4array_insert_((ccp4_ptr*)(&v),i,(ccp4_constptr)(&d),sizeof(*v))

/*! Macro to delete element[i] of an array, preserving order.
  This decrements the size. All subsequent elements are moved down. As a
  result this method is slow.
  \param v   The array pointer
  \param i   The element to be deleted
*/
#define ccp4array_delete_ordered(v,i) ccp4array_delete_ordered_((ccp4_ptr*)(&v),i,sizeof(*v))

/*! Macro to delete element[i] of an array without preserving order. The
  last element is moved into the gap, and the size is decremented.
  \param v   The array pointer
  \param i   The element to be deleted
*/
#define ccp4array_delete(v,i) ccp4array_delete_((ccp4_ptr*)(&v),i,sizeof(*v))

/*! Macro to delete the last element of an array.
  This decrements the size.
  \param v   The array pointer
*/
#define ccp4array_delete_last(v) ccp4array_delete_last_((ccp4_ptr*)(&v),sizeof(*v))

/*! Macro to return the size of the array.
  \param v   The array pointer
  \return    The size (int)
*/
#define ccp4array_size(v) ccp4array_size_((ccp4_constptr*)(&v))

/*! Macro free the array.
  All memory, including the header, is freed.
  \param v   The array pointer
*/
#define ccp4array_free(v) ccp4array_free_((ccp4_ptr*)(&v))

/** 
 * See macro ccp4array_new 
*/
ccp4_ptr ccp4array_new_(ccp4_ptr *p);
/** 
 * See macro ccp4array_new_size 
*/
ccp4_ptr ccp4array_new_size_(ccp4_ptr *p, const int size, const size_t reclen);
/** 
 * See macro ccp4array_resize 
*/
void ccp4array_resize_(ccp4_ptr *p, const int size, const size_t reclen);
/** 
 * See macro ccp4array_reserve 
*/
void ccp4array_reserve_(ccp4_ptr *p, const int size, const size_t reclen);
/** 
 * See macro ccp4array_append 
*/
void ccp4array_append_(ccp4_ptr *p, ccp4_constptr data, const size_t reclen);
/** 
 * See macro ccp4array_append_n 
*/
void ccp4array_append_n_(ccp4_ptr *p, ccp4_constptr data, const int n, const size_t reclen);
/** 
 * See macro ccp4array_append_list 
*/
void ccp4array_append_list_(ccp4_ptr *p, ccp4_constptr data, const int n, const size_t reclen);
/** 
 * See macro ccp4array_insert  
*/
void ccp4array_insert_(ccp4_ptr *p, const int i, ccp4_constptr data, const size_t reclen);
/** 
 * See macro ccp4array_delete_ordered 
*/
void ccp4array_delete_ordered_(ccp4_ptr *p, const int i, const size_t reclen);
/** 
 * See macro ccp4array_delete 
*/
void ccp4array_delete_(ccp4_ptr *p, const int i, const size_t reclen);
/** 
 * See macro ccp4array_delete_last 
*/
void ccp4array_delete_last_(ccp4_ptr *p, const size_t reclen);
/** 
 * See macro ccp4array_size 
*/
int ccp4array_size_(ccp4_constptr *p);
/** 
 * See macro ccp4array_free
*/
void ccp4array_free_(ccp4_ptr *p);

#ifdef __cplusplus
} 
#endif

#endif /* __CCP4_ARRAY_INC */
