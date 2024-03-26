/*
  Copyright (C) 2003 PWSCF group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include <c_ftn_defs.h>
#include <stdio.h>
#ifndef __MACH__
#include <malloc.h> 
#else
#include <stdlib.h> 
#endif
FORTRAN(void,MEMSTAT,memstat,(int *kilobytes),(kilobytes)) 
{
/*#if defined(HAVE_MALLINFO) && !defined(__QK_USER__) */
#ifndef __MACH__
  struct mallinfo info;
  info = mallinfo();
  *kilobytes = info.arena / 1024 ;
#else 
  *kilobytes = -1;
#endif 
}

FORTRAN(void,SHOWMEM,showmem,(),()) 
{
/*#if defined(HAVE_MALLINFO) && !defined(__QK_USER__) */
#ifndef __MACH__
  int kb0,kb1,kb2,kb3,kb4;
  struct mallinfo info;
  info = mallinfo();
  kb0 = info.arena / 1024 ;
  kb1 = info.hblkhd/1024;
  kb2 = info.uordblks/1024;
  kb3 = info.fordblks/1024;
  kb4 = info.keepcost/1024;
  printf("%d , %d, %d ,%d, %d\n",kb0,kb1,kb2,kb3,kb4);
  
#else 
  printf("mallinfo not supported!\n");
#endif 

}
