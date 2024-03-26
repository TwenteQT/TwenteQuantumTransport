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


FORTRAN(FILE *,C_FOPEN,c_fopen,( const char * filename , const int *flen,const char * mode),(filename,flen,mode)) 
{
    FILE * retval = NULL ;
    char lfile[*flen+1],lmode[2];
    memcpy( lfile , filename , *flen ) ;
    lfile[*flen] = '\0' ;

    lmode[0]=*mode;
    lmode[1] = '\0' ;

    retval = fopen( lfile , lmode ) ;
    if ( retval == NULL ){
        fprintf( stderr , "fopen fail: [%d] %s (file: %s)\n" , errno , strerror( errno ), lfile ) ;
            exit(1);
        }
    return retval ;
}
FORTRAN(int, C_FCLOSE,c_fclose,( FILE ** const fl),(fl)){
    int retval=-1;
    retval=fclose(*fl);
    if ( retval !=0 ){
        fprintf( stderr , "fclose fail: [%d] %s\n" , errno , strerror( errno ) ) ;
        exit(1);
    }
    return retval ;
}

FORTRAN(int, C_FWRITE,c_fwrite,(const void *ptr, const int * size,const int *nitems, FILE ** const fl),(ptr,size,nitems,fl))
{
    int retval=0;
    retval=fwrite(ptr, *size, *nitems, *fl);
    if ( retval < (*nitems) ){
        fprintf( stderr , "fwrite fail: [%d] %s\n" , errno , strerror( errno ) ) ;    
            exit(1);
        }
    return retval;
}

FORTRAN(int, C_FREAD,c_fread,( void  * const ptr, const int * size,const int *nitems, FILE ** const fl),(ptr,size,nitems,fl))
{
    int retval=0;
    retval=fread(ptr, *size, *nitems, *fl);
    if ( retval < (*nitems) ){
        fprintf( stderr , "fread fail: [%d] %s\n" , errno , strerror( errno ) ) ;    
            exit(1);
        }
    return retval;
}


