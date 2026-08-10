    /* This program reads input files (in the format of I.Turek interface code) from ./base subdirectory       */
    /* and sets up the inputs for the system with a number of pricipal layers (PLs) repeated. The range of PLs */
    /* to be repeated and the number of repetitions are given as arguments at the command line.                */
#include <sys/stat.h>   
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/* Various macros/functions defined below for reading/writing lines of the files etc. */

#define scanv3l(f,a,lb) fgets(buf,400,f);sscanf(buf,"%le %le %le %s",(a)+0,(a)+1,(a)+2,lb);

#define scanline(f,...) fgets(buf,400,f);sscanf(buf,""__VA_ARGS__);

#define scanv3(f,a) fgets(buf,400,f);sscanf(buf,"%le %le %le",(a)+0,(a)+1,(a)+2);
#define scanv2(f,a) fgets(buf,400,f);sscanf(buf,"%le %le",(a)+0,(a)+1);
#define fpvec3(f,a,cm) \
fprintf(f," %#15.10f %#15.10f %#15.10f   %s\n",(a)[0],(a)[1],(a)[2],cm)
#define fpvec2(f,a,cm) \
    fprintf(f," %#15.10g %#15.10g                  %s\n",(a)[0],(a)[1],cm)

#define openinp(a,nm) a=fopen(nm,"r");\
    if (a==NULL) \
{ \
    printf("Do you remeber to put \"%s\" file in current dir? :)\n",nm); \
    exit(1); \
} 


#define startwrite(fp,fn,na,nb,nmtr,c1,sf,t1,t2,t3,t4) \
fp=fopen(fn,"w");\
    fprintf(fp,"-------===========  GEOMETRY  FILE ===========-------\n");\
    fprintf(fp," %5d %3d %3d     %45s\n",na,nb,nmtr,"# NP, NB, NMTR");\
    fprintf(fp," %#15.10g %40s\n",c1,"# CUTRAT");\
    fpvec3(fp,sf,"# SCALING FACTORS");\
    fpvec2(fp,t1,"# 1ST TRANSL. VECTOR");\
    fpvec2(fp,t2,"# 2ND TRANSL. VECTOR");\
    fpvec3(fp,t3,"# PERP.TR. VECTOR");\
    fpvec3(fp,t4,"# PERP.TR. VECTOR");	
#define endwrite(fp)  \
fprintf(fp,"-------=========== END OF GEOMETRY ===========-------\n");\
    fclose(fp);
#define atfmt " %-15s   %s\n"
void getcoordz(double *pos,double *coord,double *tr, int l);
inline void getcoordz(double *pos,double *coord,double *tr, int l) {
    pos[0]=coord[0]+l*tr[0];
    pos[1]=coord[1]+l*tr[1];
    pos[2]=coord[2]+l*tr[2];
}


#define wrconc(ind) \
fprintf(FPco,"%6d\n",natk[ind]);\
for (__j=0;__j<natk[ind];++__j){\
    fprintf(FPco,"    %-10s  %-10s\n",label[ind][__j],label[ind][__j]);\
    fprintf(FPco,"%8.2f\n",conc[ind][__j]);\
}


#include <sys/time.h>
#include <math.h>


int main(int argc, char *argv[])
{
    FILE *FP,*FPg,*FPc,*FPco;
    int N,**kind,nn,i,p,nat,at,nx,ny,pind,j,l,mna,ord;
    double *coord,pos[3],pr,disp,trans[3],shift;
    char buf[400],*dirn,buf1[40],***label, **glabel;
    int __i,__j;

    int np,nb,nmtr,*natk,extbf;
    int nstart,nend,nrep;
    double scx[3],tr1[2],tr2[2],ltr[3],rtr[3],cr,mtr1[2],mtr2[2],**conc,scl;

    int rbl,lbl;

    if (argc<4) {             /* Check number of arguments. argv[0] is the executable name. */
        printf("Usage: %s Nstart Nnum Nrep\n",argv[0]);
        exit(1);
    }
    nstart=atoi(argv[1])-1;                 /*  nstart is the number of PLs before the start of the repeated slice */
    nend=atoi(argv[2])+nstart;              /*  The last PL to be repeated */
    nrep=atoi(argv[3])+1;	            /* How many repetitions */
    if (nrep<1) {
        printf("Error! Nrep must be >= 0\n");
        exit(1);
    }
    dirn=argv[3];    


    openinp(FP,"base/inpge");                          /* reading the initial geometry number of PL, sites per PL, prymitive vectors etc. */
    scanline(FP,"");
    scanline(FP,"%d %d %d",&np,&nb,&nmtr);  
    scanline(FP,"%le",&cr);
    nat=np*nb;                                         /* nat = total number of at. sites;    */  
    scanv3(FP,scx);
    scanv2(FP,tr1);
    scanv2(FP,tr2);
    scanv3(FP,ltr);
    scanv3(FP,rtr);
    coord=(double *)malloc(sizeof(double)*nat*3);            /* Coordinates */


    natk=(int *)calloc(nat,sizeof(int));                    /* natk[i]=Number of atoms per site (in case of CPA) */  
    conc=(double **)calloc(nat,sizeof(double*));            /* Concentrations. Here a vector (length nat) of pointers is allocated.*/
    label=(char ***)calloc(nat,sizeof(char **));            /* Ditto for atomic labels */
    glabel=(char **)calloc(nat,sizeof(char *));

    openinp(FPc,"base/inpch");                             /* Reading chemical composition */       
    scanline(FPc,"");
    mna=0;
    for (i=0;i<nat;i++){                                   /* Loop over atomic sites */
      scanline(FPc,"%d",natk+i);                           /* Read the number of elements per site */
      conc[i]=(double *)calloc(natk[i],sizeof(double));    /* Allocate appropriate vector holding concentrations */
      label[i]=(char **)calloc(natk[i],sizeof(char *));
        glabel[i]=(char *)calloc(30,sizeof(char));
        glabel[i][0]=0;
        for (j=0;j<natk[i];j++){                              /* Loop over different elements populating site i  */
            label[i][j]=(char *)calloc(30,sizeof(char));    
            scanline(FPc,"%s",label[i][j]);
            // printf("%s\n",label[i][j]);
            scanline(FPc,"%le",conc[i]+j);                   /* Read the concentration, conc[i]+j is pointer arithmetics, basicaly equals conc[i][j] */
            strcat(glabel[i],label[i][j]);
            mna++;
        }
    }  
    fclose(FPc);

     /* Begins to write a new inpge. So far the only change is the number of principal layers */
    startwrite(FPg,"inpge",(np+(nrep-1)*(nend-nstart)),nb,nmtr,cr,scx,tr1,tr2,ltr,rtr);

    FPco=fopen("inpch","w"); /* Opens new inpch */

    fprintf(FPco,"------ CHEMICAL OCCUPATIONS:\n");
    
    for (i=0;i<nb;i++) {     /* Copies the coordinates of the PL of the left lead.  */       
        scanv3l(FP,pos,buf1);
        fpvec3(FPg,pos,buf1);
    }

    for (i=0;i<nat;i++) {   /* Read coordinates of the central part */
        scanv3(FP,(coord+3*i));
    }


    for (i=0;i<3;++i) {     /* Compute the translation from the beginning to the end of the repeated section. */
        trans[i]=coord[3*nb*nend+i]-coord[3*nb*nstart+i];
    }

    
    for (i=0;i<nstart;i++) {  /* The coordinates of the sites in PLs before the repeated region.  */
      for (j=0;j<nb;++j){     /* Unchanged from original */
	fpvec3(FPg,(coord+3*(nb*i+j)),glabel[nb*i+j]);
	wrconc(nb*i+j);

        }
    }

    for (l=0;l<nrep;l++) {               /* The PLs from nstart to nend are repeated (with translation) nrep times. */
        for (i=nstart;i<nend;i++) {
  	  /* fprintf(FPg,"l=%d , i= %d : ",l,i); */
            for (j=0;j<nb;++j){
                getcoordz(pos,coord+3*(nb*i+j),trans,l);
                fpvec3(FPg,pos,glabel[nb*i+j]);
                wrconc(nb*i+j);

            }
        }
    }

    /* fprintf(FPg, " End of the extended region. \n "); */
    	
    for (i=nend;i<np;i++) {            /* The remaining PLs of the scattering region translated appropriately */         
        for (j=0;j<nb;++j){
            getcoordz(pos, (coord+3*(nb*i+j)),trans,nrep-1);
            fpvec3(FPg,pos,glabel[nb*i+j]);
            wrconc(nb*i+j);
        }
    }


    for (i=0;i<nb;i++) {              /* And the right lead  */
        scanv3l(FP,pos,buf1);
        getcoordz(pos, pos,trans,nrep-1);		
        fpvec3(FPg,pos,buf1);
    }
    endwrite(FPg);
    fprintf(FPco,"-------=========== END OF CHEM ===========-------\n");	
    fclose(FPco);

    openinp(FP,"base/inpbu");

    FPco=fopen("inpbu","w");                       /* inpbu describe the chemical composition of the leads */
    extbf=(fgetc(FP)=='+')?1:0;                    /* seemingly copied line by line from the original */
    fputc((extbf)?'+':' ',FPco);
    scanline(FP,"");
    fprintf(FPco,"%s",buf);	
    for (l=0;l<((extbf)?nb:1);l++){
        scanline(FP,"%d",&lbl);
        fprintf(FPco,"%s",buf);		
        for (i=0;i<lbl;++i){
            scanline(FP,""); /* read concentration here if need */
            fprintf(FPco,"%s",buf);			
            scanline(FP,"");
            fprintf(FPco,"%s",buf);			
        }
    }
    scanline(FP,"");
    fprintf(FPco,"%s",buf);		
    for (l=0;l<((extbf)?nb:1);l++){
        scanline(FP,"%d",&rbl);
        fprintf(FPco,"%s",buf);		
        for (i=0;i<rbl;++i){
            scanline(FP,""); /* read concentration here if need  */
            fprintf(FPco,"%s",buf);			
            scanline(FP,"");
            fprintf(FPco,"%s",buf);			
        }
    }
    fclose(FP);

    fprintf(FPco,"-------=========== END OF BULK ===========-------\n");
    fclose(FPco);
    /////////////////////injection//////////// /* The rest of the code copies the inpsk and inpsk2 files line by line */
    openinp(FP,"base/inpsk");                  /* from the files in ./base. I see no function here. It's likely a placeholder to be modified when needed   */
    
    FPco=fopen("inpsk","w");                 
    extbf=(fgetc(FP)=='+')?1:0;
    fputc((extbf)?'+':' ',FPco);
    scanline(FP,"");
    fprintf(FPco,"%s",buf);
    for (l=0;l<((extbf)?nb:1);l++){
        scanline(FP,"%d",&lbl);
        fprintf(FPco,"%s",buf);
        for (i=0;i<lbl;++i){
            scanline(FP,""); /* read concentration here if need */
            fprintf(FPco,"%s",buf);
            scanline(FP,"");
            fprintf(FPco,"%s",buf);
        }
    }
    fclose(FP);
    
    fprintf(FPco,"-------=========== END OF SINK ===========-------\n");
    fclose(FPco);
 /////////////////////injection////////////
    openinp(FP,"base/inpsk2");
    
    FPco=fopen("inpsk2","w");
    extbf=(fgetc(FP)=='+')?1:0;
    fputc((extbf)?'+':' ',FPco);
    scanline(FP,"");
    fprintf(FPco,"%s",buf);
    for (l=0;l<((extbf)?nb:1);l++){
        scanline(FP,"%d",&lbl);
        fprintf(FPco,"%s",buf);
        for (i=0;i<lbl;++i){
            scanline(FP,""); /* read concentration here if need */
            fprintf(FPco,"%s",buf);
            scanline(FP,"");
            fprintf(FPco,"%s",buf);
        }
    }
    fclose(FP);
    
    fprintf(FPco,"-------=========== END OF SINK number two ===========-------\n");
    fclose(FPco);

}
