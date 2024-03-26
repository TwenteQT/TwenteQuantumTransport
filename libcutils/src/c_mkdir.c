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


FORTRAN(int,C_MKDIR,c_mkdir,( const char * dirname , const int * length ),(dirname,length)) 
{
   int retval = -1 ;
   mode_t mode = 0777 ;
   char * ldir = ( char * ) malloc( (*length) + 1 ) ;
   memcpy( ldir , dirname , *length ) ;
   ldir[*length] = '\0' ;
   retval = mkdir( ldir , mode ) ;
   if ( retval == -1  && errno != EEXIST )
     fprintf( stderr , "mkdir fail: [%d] %s\n" , errno , strerror( errno ) ) ;
   free( ldir ) ;
   return retval ;
}

FORTRAN(int,C_RENAME,c_rename,( const char * oldname , const int * len1, const char * newname, const int * len2 ),(oldname,len1,newname,len2)) 
{
   int retval = -1 ;
   mode_t mode = 0777 ;
   char * oname = ( char * ) malloc( (*len1) + 1 ) ;
   char * nname = ( char * ) malloc( (*len2) + 1 ) ;
   
   
   memcpy( oname , oldname , *len1 ) ;
   memcpy( nname , newname , *len2 ) ;
   
   oname[*len1] = '\0' ;
   nname[*len2] = '\0' ;
   
   retval = rename( oname , nname ) ;
   if ( retval == -1  && errno != EEXIST )
     fprintf( stderr , "rename fail: [%d] %s\n" , errno , strerror( errno ) ) ;
   free( oname ) ;
   free( nname ) ;
   return retval ;
}
