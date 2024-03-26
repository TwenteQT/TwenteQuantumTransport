#include <sys/stat.h>   
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define scanv3l(f,a,lb) fgets(buf,400,f);sscanf(buf,"%le %le %le %s",(a)+0,(a)+1,(a)+2,lb);

#define scanline(f,...) fgets(buf,400,f);sscanf(buf,""__VA_ARGS__);

#define fpvec3(f,a,cm) \
fprintf(f," %#15.11g %#15.11g %#15.11g   %s\n",(a)[0],(a)[1],(a)[2],cm)

#define getcoordz(pos,coord,tr,l) \
    pos[0]=(coord)[0]+(l)*((tr)[0]);\
    pos[1]=(coord)[1]+(l)*((tr)[1]);\
    pos[2]=(coord)[2]+(l)*((tr)[2]);
#define endwrite(fp)  \
fprintf(fp,"-------=========== END OF GEOMETRY ===========-------\n");\


#include <sys/time.h>
#include <math.h>


    int main(int argc, char *argv[])
{
    FILE *FPi,*FPo;
    char buf[400];
    int na,nmtr,tp,na1,nadd,nsz,nra;
    int nstart,nend,nrep;
    char **lbl,lblnext[20];
    double **pos,posnext[3], trv[3],cpos[3];
    size_t i,j;
    // double scx[3],tr1[2],tr2[2],ltr[3],rtr[3],cr,mtr1[2],mtr2[2],**conc,scl;

    // int rbl,lbl;

    if (argc<4) {
        printf("Usage: %s Nstart Nnum Nrep\n",argv[0]);
        exit(1);
    }
    nstart=atoi(argv[1])-1;
    nend=atoi(argv[2])+nstart;
    nrep=atoi(argv[3]);	
    if (nrep<0) {
        printf("Error! Nrep must be >=0\n");
        exit(1);
    }

    nsz=atoi(argv[2]);	
    nadd=nsz*atoi(argv[3]);

    pos=(double**)malloc(sizeof(double*)*nsz);
    lbl=(char **)malloc(sizeof(char *)*nsz);
    for( i = 0; i < nsz; ++i)
    {
        pos[i]=(double*)malloc(sizeof(double)*3);
        lbl[i]=(char*)malloc(sizeof(char)*20);
    }

    FPi=fopen("base/atomlist","r");
    FPo=fopen("atomlist","w");

    while (!feof(FPi)){
        fgets(buf,400,FPi);
        fputs(buf,FPo);
    }

    fclose(FPi);
    fclose(FPo);


    FPi=fopen("base/geom_l","r");
    FPo=fopen("geom_l","w");

    while (feof(FPi)==0){
        fgets(buf,400,FPi);
        fputs(buf,FPo);
    }

    fclose(FPi);
    fclose(FPo);

    FPi=fopen("base/geom_m","r");
    FPo=fopen("geom_m","w");

    fgets(buf,400,FPi);
    fputs(buf,FPo);

    scanline(FPi,"%d %d %d", &na,&nmtr,&tp);
    na1=na+nadd;
    fprintf(FPo," %5d   %5d   %5d %45s\n",na1,nmtr,tp,"# NP,  NMTR,TR_PERP");
    for( i = 0; i < 6; ++i)
    {
        fgets(buf,400,FPi);
        fputs(buf,FPo);
    }


    for( i = 0; i < nstart; ++i)
    {
        scanv3l(FPi,posnext,lblnext);
        fpvec3(FPo,posnext,lblnext);            
    }
    // fprintf(FPo,"\n");

    for( i = 0; i < nsz; ++i)
    {
        scanv3l(FPi,pos[i],lbl[i]);
        fpvec3(FPo,pos[i],lbl[i]);
    }
    // fprintf(FPo,"\n");

    scanv3l(FPi,posnext,lblnext);
    trv[0]=posnext[0]-pos[0][0];
    trv[1]=posnext[1]-pos[0][1];
    trv[2]=posnext[2]-pos[0][2];
    printf("trans.vec=[%lg,%lg,%lg]\n",trv[0],trv[1],trv[2]);
    for( i = 0; i < nrep; ++i)
    {
        for( j = 0; j < nsz; ++j)
        {
            getcoordz(cpos,pos[j],trv,i+1);
            fpvec3(FPo,cpos,lbl[j]);            
        }
    }

    // fprintf(FPo,"\n");
    getcoordz(cpos,posnext,trv,nrep);
    fpvec3(FPo,cpos,lblnext);    

    for( i = nend+1; i < na; ++i)
    {
        scanv3l(FPi,posnext,lblnext);
        getcoordz(cpos,posnext,trv,nrep);
        fpvec3(FPo,cpos,lblnext);            
    }

    endwrite(FPo);

    fclose(FPi);
    fclose(FPo);


    FPi=fopen("base/geom_r","r");
    FPo=fopen("geom_r","w");

    fgets(buf,400,FPi);
    fputs(buf,FPo);

    scanline(FPi,"%d %d %d", &nra,&nmtr,&tp);
    fprintf(FPo," %5d   %5d   %5d %45s\n",nra,nmtr,tp,"# NP,  NMTR,TR_PERP");

    for( i = 0; i < 6; ++i)
    {
        fgets(buf,400,FPi);
        fputs(buf,FPo);
    }



    for( i = 0; i < nra; ++i)
    {
        scanv3l(FPi,posnext,lblnext);
        getcoordz(cpos,posnext,trv,nrep);
        fpvec3(FPo,cpos,lblnext);            

    }	
    endwrite(FPo);

    fclose(FPi);
    fclose(FPo);


}