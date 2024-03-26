#include <sys/stat.h>   
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <string.h>

#define scanv3l(f,a,lb) fgets(buf,400,f);sscanf(buf,"%le %le %le %s",(a)+0,(a)+1,(a)+2,lb);

#define scanline(f,...) fgets(buf,400,f);sscanf(buf,""__VA_ARGS__);
#define scanv3(f,a) fgets(buf,400,f);sscanf(buf,"%le %le %le",(a)+0,(a)+1,(a)+2);
#define scanv2(f,a) fgets(buf,400,f);sscanf(buf,"%le %le",(a)+0,(a)+1);
#define fpvec3(f,a,cm) \
fprintf(f," %#15.10g%#15.10g%#15.10g   %s\n",(a)[0],(a)[1],(a)[2],cm)
#define fpvec2(f,a,cm) \
    fprintf(f," %#15.10g%#15.10g                  %s\n",(a)[0],(a)[1],cm)

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
void getshift(double *shift,double *pos,double *pos1);


inline void getcoordz(double *pos,double *coord,double *tr, int l) {
    pos[0]=coord[0]+l*tr[0];
    pos[1]=coord[1]+l*tr[1];
    pos[2]=coord[2]+l*tr[2];
}

inline void getshift(double *shift,double *pos,double *pos1) {
    shift[0]=pos[0]-pos1[0];
    shift[1]=pos[1]-pos1[1];
    shift[2]=pos[2]-pos1[2];
}


#define wrconc(ind) \
fprintf(FPco,"%6d\n",natk[ind]);\
for (__j=0;__j<natk[ind];++__j){\
    fprintf(FPco,"    %-10s  %-10s\n",label[ind][__j],ilabel[ind][__j]);\
    fprintf(FPco,"%8.2f        %5.2f \n",conc[ind][__j],zval[ind][__j]);\
}


#include <sys/time.h>
#include <math.h>
#if !defined(stpcpy)                                                            
char *stpcpy(char *dest, const char *src){                                      
    strcpy(dest,src);                                                           
    return(dest+strlen(src));                                                   
}                                                                               
#endif

int main(int argc, char *argv[])
{
    FILE *FP,*FPg,*FPc,*FPco;
    int i,nat,nat1,j,mna;
    double *coord,pos[3],trans[3];
    char buf[400],buf1[40],***label, **glabel,***ilabel;
    int __i,__j;
    int side,atofs;
    char cside,dirn[20],bufr[1024],*bufp,bulklbl[20],abulklbl[20];
    double *posl,*posr,*posptr,*shift,**zval;
    int np,nb,nmtr,*natk,extbf;
    double scx[3],tr1[2],tr2[2],ltr[3],rtr[3],cr,**conc;
    double mbound[3],mbcorr[3];
    int rbl,lbl;

    if (argc<2 || !(argv[1][0]=='r' || argv[1][0]=='l')) {
        printf("Usage: %s [r|l]\n",argv[0]);
        exit(1);
    }

    if (argc==4){
        mbcorr[0]=atoi(argv[2]);
        mbcorr[1]=atoi(argv[3]);
        mbcorr[2]=atoi(argv[4]);
        
    } else{
        mbcorr[0]=0.0;
        mbcorr[1]=0.0;
        mbcorr[2]=0.0;
    }
    
    side=-1;
    if (argv[1][0]=='r') side=1;
    cside=argv[1][0];

    sprintf(dirn,"%c-mirror",cside);
    mkdir(dirn,0777);

    openinp(FP,"inpbu");

    sprintf(buf,"%s/inpbu",dirn);
    FPco=fopen(buf,"w");
    bufp=bufr;

    extbf=(fgetc(FP)=='+')?1:0;
    if (extbf==1) {
        printf("This program doesn't support\nextended bulk format!");
        exit(0);

    }
    fputc((extbf)?'+':' ',FPco);
    if (side==-1) {
        scanline(FP,"");
        scanline(FP,"%d",&lbl);
        for (i=0;i<lbl;++i){
            scanline(FP,""); /* read concentration here if need */
            scanline(FP,"%s",abulklbl);
        }
        fgetc(FP);
    }
    bufp[0]=' ';
    bufp++;        

    scanline(FP,"");
    fprintf(FPco,buf);		
    bufp=stpcpy(bufp,buf)+1;
    scanline(FP,"%d",&rbl);		
    fprintf(FPco,buf);	
    bufp=stpcpy(bufp,buf)+1;		
    for (i=0;i<rbl;++i){
        scanline(FP,""); /* read concentration here if need  */
        fprintf(FPco,buf);	
        bufp=stpcpy(bufp,buf)+1;						
        // scanline(FP);
        scanline(FP,"%s",bulklbl);
        fprintf(FPco,buf);			
        bufp=stpcpy(bufp,buf)+1;
        // strncpy(bulklbl,buf,20);             
        // bulklbl[19]=0;
    }
    bufp=bufr;
    bufp+=(fprintf(FPco,bufp)+1);		
    bufp+=(fprintf(FPco,bufp)+1);		
    for (i=0;i<rbl;++i){
        bufp+=(fprintf(FPco,bufp)+1);		
        bufp+=(fprintf(FPco,bufp)+1);		
    }

    if (side==1) {
        scanline(FP,"");
        scanline(FP,"%d",&lbl);
        for (i=0;i<lbl;++i){
            scanline(FP,""); /* read concentration here if need */
            scanline(FP,"%s",abulklbl);
        }
    }

    fclose(FP);

    fprintf(FPco,"-------=========== END OF BULK ===========-------\n");
    fclose(FPco);

    openinp(FP,"inpge");
    scanline(FP,"");
    scanline(FP,"%d %d %d",&np,&nb,&nmtr);  
    scanline(FP,"%le",&cr);
/*	if (extbf &&(nb!=nblbl[0] || nb!=nblbl[1])) {
    printf("Error! For extended bulk file format\n");
    printf("number of atoms in leads should be\n");
    printf("equal nuber of atoms in planes\n");
    printf("Stop.\n");		
    exit(1);
}	
*/
nat=np*nb;  
scanv3(FP,scx);
scanv2(FP,tr1);
scanv2(FP,tr2);
scanv3(FP,ltr);
scanv3(FP,rtr);

nat1=(np*2+1)*nb;
coord=(double *)malloc(sizeof(double)*nat1*3);
zval=(double **)malloc(sizeof(double*)*nat1);
shift=(double *)malloc(sizeof(double)*nat*3);


natk=(int *)calloc(nat1,sizeof(int));
conc=(double **)calloc(nat1,sizeof(double*));
label=(char ***)calloc(nat1,sizeof(char **));
ilabel=(char ***)calloc(nat1,sizeof(char **));
glabel=(char **)calloc(nat1,sizeof(char *));

openinp(FPc,"inpch");
scanline(FPc,"");
atofs=0;
if (side==-1) atofs=(np+1)*nb;

mna=0;
for (i=0;i<nat;i++){
    scanline(FPc,"%d",natk+i+atofs);
    conc[i+atofs]=(double *)calloc(natk[i+atofs],sizeof(double));
    zval[i+atofs]=(double *)calloc(natk[i+atofs],sizeof(double));
    label[i+atofs]=(char **)calloc(natk[i+atofs],sizeof(char *));
    ilabel[i+atofs]=(char **)calloc(natk[i+atofs],sizeof(char *));
    glabel[i+atofs]=(char *)calloc(30,sizeof(char));
    glabel[i+atofs][0]=0;
    for (j=0;j<natk[i+atofs];j++){
        label[i+atofs][j]=(char *)calloc(30,sizeof(char));      
        ilabel[i+atofs][j]=(char *)calloc(30,sizeof(char));      
        scanline(FPc,"%s %s",label[i+atofs][j],ilabel[i+atofs][j]);
        scanline(FPc,"%le %le",conc[i+atofs]+j,zval[i+atofs]+j);
        strcat(glabel[i+atofs],label[i+atofs][j]);
        mna++;
    }
}  
fclose(FPc);


if (side==1){
    // mbound[0]=rtr[0];
    // mbound[1]=rtr[1];
    // mbound[2]=rtr[2];
    rtr[0]=-ltr[0];
    rtr[1]=-ltr[1];
    rtr[2]=-ltr[2];    
}
else
{
    // mbound[0]=rtr[0];
    // mbound[1]=rtr[1];
    // mbound[2]=rtr[2];
    ltr[0]=-rtr[0];
    ltr[1]=-rtr[1];
    ltr[2]=-rtr[2];
}
sprintf(buf,"%s/inpge",dirn);    
startwrite(FPg,buf,np*2+1,nb,nmtr,cr,scx,tr1,tr2,ltr,rtr);
sprintf(buf,"%s/inpch",dirn);
FPco=fopen(buf,"w");

fprintf(FPco,"------ CHEMICAL OCCUPATIONS:\n");

posl=(double*)malloc(sizeof(double)*nb*3);
posr=(double*)malloc(sizeof(double)*nb*3);

for (i=0;i<nb;i++) {
    scanv3(FP,(posl+3*i));
}

for (i=0;i<nat;i++) {
    scanv3(FP,(coord+3*(i+atofs)));
}

for (i=0;i<nb;i++) {
    scanv3(FP,(posr+3*i));
}

    // for (i=0;i<nat;i++) {
    //  scanv3(FP,(coord+3*(i+atofs));
    // }

    // printf("%d , %le, %le\n",nb,posr[0],posr[1]);
if (side==1){
    mbound[0]=mbound[1]=mbound[2]=0.0;
    for (i=0;i<nb;i++) {
        fpvec3(FPg,(posl+3*i),bulklbl);
    }
    for (i=0;i<nb;i++) {
        getcoordz((coord+3*(nat+i)),(posr+3*i),shift,0);
        glabel[i+nat]=abulklbl;
        natk[i+nat]=1;
        conc[i+nat]=(double *)calloc(1,sizeof(double));
        conc[i+nat][0]=1.0;
        zval[i+nat]=(double *)calloc(1,sizeof(double));
        zval[i+nat][0]=0.0;
        label[i+nat]=(char **)calloc(1,sizeof(char *));
        ilabel[i+nat]=(char **)calloc(1,sizeof(char *));
        label[i+nat][0]=abulklbl;
        ilabel[i+nat][0]=abulklbl;
        mbound[0]+=(posr+3*i)[0];
        mbound[1]+=(posr+3*i)[1];
        mbound[2]+=(posr+3*i)[2];

    }
    mbound[0]/=(double)nb+mbcorr[0];
    mbound[1]/=(double)nb+mbcorr[1];
    mbound[2]/=(double)nb+mbcorr[2];
        // printf("%d , %le, %le\n",nb,posr[0],mbound[2]);

    for (i=nat+nb;i<nat1;i++){
        atofs=nat-(i-nat-nb)-1;
        getshift(shift,mbound,coord+3*(atofs));
            // printf("%le, %le, %le\n",shift[0],shift[1],shift[2]);            
        getcoordz((coord+3*(i)),mbound,shift,1);
        label[i]=label[atofs];
        ilabel[i]=ilabel[atofs];
        natk[i]=natk[atofs];
        conc[i]=conc[atofs];
        zval[i]=zval[atofs];
        glabel[i]=glabel[atofs];
    }

        // for (i=nat+nb-1;i<nat1;i++){
        //     atofs=nat-(i-nat-nb)-1;
        //     getshift(shift,coord+3*(atofs),coord+3*(atofs+1));
        //     getcoordz((coord+3*(i)),(coord+3*(i-1)),shift,-1);
        //     label[i]=label[atofs];
        //     natk[i]=natk[atofs];
        //     conc[i]=conc[atofs];
        //     glabel[i]=glabel[atofs];
        // }

            //         for (i=nat+nb;i<;i++) {
            //          for (j=0;j<nb;++j){
            //              
            //              
            //              fpvec3(FPg,(coord+3*(nb*i+j)),glabel[nb*i+j]);
            //                     wrconc(nb*i+j);
            //                 }
            //             }
            // }
    for (i=0;i<np*2+1;i++) {
        for (j=0;j<nb;++j){
            fpvec3(FPg,(coord+3*(nb*i+j)),glabel[nb*i+j]);
            wrconc(nb*i+j);

        }
    }
    for (i=nb-1;i>=0;i--) {
        getshift(shift,mbound,(posl+3*i));
        getcoordz((coord),mbound,shift,1);
        fpvec3(FPg,(coord),bulklbl);
    }

}
else{
    mbound[0]=mbound[1]=mbound[2]=0.0;

    for (i=0;i<nb;i++) {
        getcoordz((coord+3*(nat+i)),(posl+3*i),shift,0);
        glabel[i+nat]=abulklbl;
        natk[i+nat]=1;
        conc[i+nat]=(double *)calloc(1,sizeof(double));
        conc[i+nat][0]=1.0;
        zval[i+nat]=(double *)calloc(1,sizeof(double));
        zval[i+nat][0]=1.0;        
        label[i+nat]=(char **)calloc(1,sizeof(char *));
        ilabel[i+nat]=(char **)calloc(1,sizeof(char *));
        label[i+nat][0]=abulklbl;
        ilabel[i+nat][0]=abulklbl;
        mbound[0]+=(posl+3*i)[0];
        mbound[1]+=(posl+3*i)[1];
        mbound[2]+=(posl+3*i)[2];

    }
    mbound[0]/=(double)nb;
    mbound[1]/=(double)nb;
    mbound[2]/=(double)nb;
    for (i=nb-1;i>=0;i--) {
        getshift(shift,mbound,(posr+3*i));
        getcoordz((coord),mbound,shift,1);
        fpvec3(FPg,(coord),bulklbl);
    }

            // printf("%d , %le, %le\n",nb,posr[0],mbound[2]);

    for (i=nat+nb;i<nat1;i++){
        atofs=nat-(i-nat-nb)-1;
        getshift(shift,mbound,coord+3*(i));
                // printf("%le, %le, %le\n",shift[0],shift[1],shift[2]);            
        getcoordz((coord+3*(atofs)),mbound,shift,1);
        ilabel[atofs]=ilabel[i];
        label[atofs]=label[i];
        natk[atofs]=natk[i];
        conc[atofs]=conc[i];
        zval[atofs]=zval[i];
        glabel[atofs]=glabel[i];
    }
    // exit(0);

            // for (i=nat+nb-1;i<nat1;i++){
            //     atofs=nat-(i-nat-nb)-1;
            //     getshift(shift,coord+3*(atofs),coord+3*(atofs+1));
            //     getcoordz((coord+3*(i)),(coord+3*(i-1)),shift,-1);
            //     label[i]=label[atofs];
            //     natk[i]=natk[atofs];
            //     conc[i]=conc[atofs];
            //     glabel[i]=glabel[atofs];
            // }

                //         for (i=nat+nb;i<;i++) {
                //          for (j=0;j<nb;++j){
                //              
                //              
                //              fpvec3(FPg,(coord+3*(nb*i+j)),glabel[nb*i+j]);
                //                     wrconc(nb*i+j);
                //                 }
                //             }
                // }
    for (i=0;i<np*2+1;i++) {
        for (j=0;j<nb;++j){
            fpvec3(FPg,(coord+3*(nb*i+j)),glabel[nb*i+j]);
            wrconc(nb*i+j);

        }
    }
    for (i=0;i<nb;i++) {
        fpvec3(FPg,(posr+3*i),bulklbl);
    }
    

}




    //     for (i=0;i<np*2+1;i++) {
    //      for (j=0;j<nb;++j){
    //          fpvec3(FPg,(coord+3*(nb*i+j)),glabel[nb*i+j]);
    //                 wrconc(nb*i+j);
    // 
    //  }
    // }

    // 
    // for (i=0;i<nb;i++) {
    //  scanv3l(FP,pos,buf1);
    //  getcoordz(pos, pos,trans,nrep-1);       
    //  fpvec3(FPg,pos,buf1);
    // }
endwrite(FPg);
fprintf(FPco,"-------=========== END OF CHEM ===========-------\n");	
fclose(FPco);
sprintf(buf,"%s/atoms",dirn);    
symlink("../atoms",buf);
}
