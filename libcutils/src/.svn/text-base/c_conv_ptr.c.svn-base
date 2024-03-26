/*
  Copyright (C) 2003 PWSCF group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <time.h>
#include <c_ftn_defs.h>


FORTRAN(int,C_IPTR,c_iptr,(int* source, int * target ),(source,target)) 
{

   int retval=1;
   int **q;
   q=&source;
   *q=target;
   return retval;
}
