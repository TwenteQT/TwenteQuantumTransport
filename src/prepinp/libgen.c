#include "libgen.h"
#include "ldl.h"
#include <stdarg.h>


int LOG_MAX=0;

#define COMPLAIN_EXIT {LOGMSG(0,"Wrong input!");exit(1);};

inline void scanv3(char *buf,FILE *f,double *a ){
    fgets(buf,400,f);if (sscanf(buf,"%le %le %le",a+0,a+1,a+2)!=3) COMPLAIN_EXIT;
}

inline void scanv2(char *buf,FILE *f,double *a ){
    fgets(buf,400,f);if(sscanf(buf,"%le %le",a+0,a+1)!=2) COMPLAIN_EXIT;    
}

void loadmasks(t_scgeom *geom,int irm,int ivm){
    FILE *FP;
    int i,n1,n2;
    char buf[400];

    // n1=geom->nb;
    // n2=geom->tnat=(geom->np+1)*geom->nb;  
    /* printf("%d, %d\n",ivm,irm);*/
    
    if (ivm){
    FP=openread("ivmask");
    for(i = 0; i < geom->tnat; ++i)
    {
        fgets(buf,400,FP);
        if (sscanf(buf,"%le",&(geom->site[i].vmask))!=1) {
            COMPLAIN_EXIT;
        }
/*        printf("ms[%d]=%le\n",i,geom->site[i].vmask);*/
    }
    fclose(FP);
    }

    if (irm){
    FP=openread("irmask");
    for(i = 0; i < geom->tnat; ++i)
    {
        fgets(buf,400,FP);
        if (sscanf(buf,"%le",&(geom->site[i].rmask))!=1) {
            COMPLAIN_EXIT;
        }
    }
    fclose(FP);
    }
    
}
int readinput(t_scgeom *geom, t_atomset *atoms){    
    FILE *FP;
    char buf[400], label[30];
    int i,n, extbf;
    t_scsite *site, *site0, *site1;
    double rawsr,dawsr,vol;

    atoms->ptr=NULL;
    atoms->num=0;
// ADDD EOF DETECTION EVERYWHERE!!!!!!!!!!!!!!!!
//     Read "inpge"
    LOGMSG(2,"Reading INPGE...");

    FP=openread("inpge");
    fgets(buf,400,FP);
    fgets(buf,400,FP);
    if(sscanf(buf,"%d %d %d",&(geom->np),&(geom->nb),&(geom->nmtr))!=3) COMPLAIN_EXIT;
    fgets(buf,400,FP);
    if (sscanf(buf,"%le",&(geom->cr))!=1) COMPLAIN_EXIT;
    geom->nat=geom->np*geom->nb;  
    geom->tnat=(geom->np+2)*geom->nb;  

    scanv3(buf,FP,geom->scx);
    scanv2(buf,FP,geom->trv[0]);
    scanv2(buf,FP,geom->trv[1]);
    scanv3(buf,FP,geom->ltrans);
    scanv3(buf,FP,geom->rtrans);
    geom->ltrans[0]=-geom->ltrans[0];
    geom->ltrans[1]=-geom->ltrans[1];
    geom->ltrans[2]=-geom->ltrans[2];

    geom->site=(t_scsite *)malloc(sizeof(t_scsite)*geom->tnat);
    geom->plane=(t_scsite **)malloc(sizeof(t_scsite *)*(geom->np+2));

    for (i=0;i<(geom->tnat);i++){
        scanv3(buf,FP,geom->site[i].coord);
        geom->site[i].vmask=1.0;
        geom->site[i].rmask=1.0;
    };   
    fclose(FP);
    geom->plane[0]=geom->site;
    for (i=0;i<=(geom->np);i++){
        geom->plane[i+1]=(geom->plane[i])+geom->nb;
    };               

    LOGMSG(2,"Done with INPGE!");


//     Read "inpch"
    LOGMSG(2,"Reading INPCH...");
    FP=openread("inpch");
    fgets(buf,400,FP);

    for (site=(geom->site+geom->nb);site<(geom->site+geom->nb+geom->nat);site++){
        LOGMSG(3,"Processing atom %d",site);
        fgets(buf,400,FP);
        if(sscanf(buf,"%d",&(site->nc))!=1) COMPLAIN_EXIT;
        LOGMSG(3,"Number of chemical elements on site: %d",site->nc);        
        site->conc=(double *)calloc(site->nc,sizeof(double));
        site->at=(t_atom **)calloc(site->nc,sizeof(t_atom *));
        for (i=0;i<site->nc;i++){
            LOGMSG(3,"Loading %d-th concentration...",i);        
            fgets(buf,400,FP);
            if(sscanf(buf,"%s",label)!=1) COMPLAIN_EXIT;
            LOGMSG(3,"Label = %s",label);        
            fgets(buf,400,FP);
            if(sscanf(buf,"%le",(site->conc)+i)!=1) COMPLAIN_EXIT;
            LOGMSG(3,"Concentration = %le",site->conc[i]);
            site->at[i]=getatom(atoms, label);
        }

    };   
    fclose(FP);
    LOGMSG(2,"Done with INPCH!");

//     Read "inpbu"
    LOGMSG(2,"Loading INPBU...");
    FP=openread("inpbu");
    extbf=(fgetc(FP)=='+')?1:0;

    for (n=0;n<2;n++){
        fgets(buf,400,FP); //lead name line
        site0=NULL;
        for (site=geom->plane[(geom->np+1)*n];site<(geom->plane[(geom->np+1)*n]+geom->nb);site++){
            if ( (extbf==1) || (site0==NULL)){
                fgets(buf,400,FP);
                if(sscanf(buf,"%d",&(site->nc))!=1) COMPLAIN_EXIT;
                site->conc=(double *)calloc(site->nc,sizeof(double));
                site->at=(t_atom **)calloc(site->nc,sizeof(t_atom *));
                for (i=0;i<site->nc;i++){
                    fgets(buf,400,FP);
                    if(sscanf(buf,"%le",(site->conc)+i)!=1) COMPLAIN_EXIT;
                    fgets(buf,400,FP);
                    if(sscanf(buf,"%s",label)!=1) COMPLAIN_EXIT;
                    site->at[i]=getatom(atoms, label);
                };
                site0=site;
            } else {
                site->at=site0->at;  
                site->conc=site0->conc;
                site->nc=site0->nc;                
            };
        };
    };
    fclose(FP);    
    LOGMSG(2,"Done with INPBU!");

    // for (i=0;i<geom->tnat;i++){
    //         for (n=0;n<geom->site[i].nc;n++){
    //             printf("%s ", geom->site[i].at[n]->label);
    //         }
    //         printf("\n");
    // }
    //calculate scaling factor
    rawsr=0.0;
//    LOGMSG(7,"test vars: %d %d ",geom->site,geom->site+geom->nb);
    for(site=(geom->site);site<((geom->site)+(geom->nb));site++){
//	LOGMSG(7,"test vars2: %d %le %s",site,site->at[0]->wsr,site->at[0]->label);
        rawsr += pow(site->at[0]->wsr,3);        
    }
    rawsr = pow((rawsr/geom->nb), (1.0/3.0));
    vol = fabs ((geom->trv[0][0]*geom->trv[1][1]-geom->trv[1][0]*geom->trv[0][1])*geom->ltrans[2]*geom->scx[0]*geom->scx[1]*geom->scx[2]);
    dawsr=pow((0.75*vol/(geom->nb*M_PI)),(1.0/3.0));
    geom->dimcf=rawsr/dawsr;
    // printf("%f\n", geom->dimcf);
    // for(site=geom->site;site<geom->site+geom->tnat;site++){
    //     printf("%s\n",site->at[0]->label);
    // }

    // for (i=0;i<=(geom->np+1);i++){
        // geom->plane[i+1]=(geom->plane[i])+geom->nb;
        // printf("%s\n",(geom->plane[i]+1)->at[0]->label);
    // };               
    
    return(1);
    
}
/////////////// FOR INJECTION
int readinput2(t_scgeom *geom, t_atomset *atoms){
    FILE *FP;
    char buf[400], label[30];
    int i,n, extbf;
    t_scsite *site, *site0;
    t_sksite *sksite, *site1;
    double rawsr,dawsr,vol;
    
    atoms->ptr=NULL;
    atoms->num=0;
    // ADDD EOF DETECTION EVERYWHERE!!!!!!!!!!!!!!!!
    //     Read "inpge"
    LOGMSG(2,"Reading INPGE...");
    
    FP=openread("inpge");
    fgets(buf,400,FP);
    fgets(buf,400,FP);
    if(sscanf(buf,"%d %d %d",&(geom->np),&(geom->nb),&(geom->nmtr))!=3) COMPLAIN_EXIT;
    fgets(buf,400,FP);
    if (sscanf(buf,"%le",&(geom->cr))!=1) COMPLAIN_EXIT;
    geom->nat=geom->np*geom->nb;
    geom->tnat=(geom->np+2)*geom->nb;
    scanv3(buf,FP,geom->scx);
    scanv2(buf,FP,geom->trv[0]);
    scanv2(buf,FP,geom->trv[1]);
    scanv3(buf,FP,geom->ltrans);
    scanv3(buf,FP,geom->rtrans);
    geom->ltrans[0]=-geom->ltrans[0];
    geom->ltrans[1]=-geom->ltrans[1];
    geom->ltrans[2]=-geom->ltrans[2];
    
    geom->site=(t_scsite *)malloc(sizeof(t_scsite)*geom->tnat);
    geom->sitesk=(t_sksite *)malloc(sizeof(t_sksite)*1);
    geom->sitesk2=(t_sksite *)malloc(sizeof(t_sksite)*1);
    geom->plane=(t_scsite **)malloc(sizeof(t_scsite *)*(geom->np+2));
    
    for (i=0;i<(geom->tnat);i++){
        scanv3(buf,FP,geom->site[i].coord);
        geom->site[i].vmask=1.0;
        geom->site[i].rmask=1.0;
    };
    fclose(FP);
    geom->plane[0]=geom->site;
    for (i=0;i<=(geom->np);i++){
        geom->plane[i+1]=(geom->plane[i])+geom->nb;
    };
    
    LOGMSG(2,"Done with INPGE!");
    
    
    //     Read "inpch"
    LOGMSG(2,"Reading INPCH...");
    FP=openread("inpch");
    fgets(buf,400,FP);
    
    for (site=(geom->site+geom->nb);site<(geom->site+geom->nb+geom->nat);site++){
        LOGMSG(3,"Processing atom %d",site);
        fgets(buf,400,FP);
        if(sscanf(buf,"%d",&(site->nc))!=1) COMPLAIN_EXIT;
        LOGMSG(3,"Number of chemical elements on site: %d",site->nc);
        site->conc=(double *)calloc(site->nc,sizeof(double));
        site->at=(t_atom **)calloc(site->nc,sizeof(t_atom *));
        for (i=0;i<site->nc;i++){
            LOGMSG(3,"Loading %d-th concentration...",i);
            fgets(buf,400,FP);
            if(sscanf(buf,"%s",label)!=1) COMPLAIN_EXIT;
            LOGMSG(3,"Label = %s",label);
            fgets(buf,400,FP);
            if(sscanf(buf,"%le",(site->conc)+i)!=1) COMPLAIN_EXIT;
            LOGMSG(3,"Concentration = %le",site->conc[i]);
            site->at[i]=getatom(atoms, label);
        }
        
    };
    fclose(FP);
    LOGMSG(2,"Done with INPCH!");
    
    //     Read "inpbu"
    LOGMSG(2,"Loading INPBU...");
    FP=openread("inpbu");
    extbf=(fgetc(FP)=='+')?1:0;
    
    for (n=0;n<2;n++){
        fgets(buf,400,FP); //lead name line
        site0=NULL;
        for (site=geom->plane[(geom->np+1)*n];site<(geom->plane[(geom->np+1)*n]+geom->nb);site++){
            if ( (extbf==1) || (site0==NULL)){
                fgets(buf,400,FP);
                if(sscanf(buf,"%d",&(site->nc))!=1) COMPLAIN_EXIT;
                site->conc=(double *)calloc(site->nc,sizeof(double));
                site->at=(t_atom **)calloc(site->nc,sizeof(t_atom *));
                for (i=0;i<site->nc;i++){
                    fgets(buf,400,FP);
                    if(sscanf(buf,"%le",(site->conc)+i)!=1) COMPLAIN_EXIT;
                    fgets(buf,400,FP);
                    if(sscanf(buf,"%s",label)!=1) COMPLAIN_EXIT;
                    site->at[i]=getatom(atoms, label);
                };
                site0=site;
            } else {
                site->at=site0->at;
                site->conc=site0->conc;
                site->nc=site0->nc;
            };
        };
    };
    fclose(FP);
    LOGMSG(2,"Done with INPBU!");
    //////////////for spin injection//////////
    //     Read "inpsk"
    LOGMSG(2,"Loading INPSK...");
    FP=openread("inpsk");
        extbf=(fgetc(FP)=='+')?1:0;
        fgets(buf,400,FP); //sink name line
        site1=NULL;
   sksite= geom->sitesk;
            if ( (extbf==1) || (site1==NULL)){
                fgets(buf,400,FP);
                if(sscanf(buf,"%d",&(sksite->nc))!=1) COMPLAIN_EXIT;
                sksite->conc=(double *)calloc(sksite->nc,sizeof(double));
                sksite->at=(t_atom **)calloc(sksite->nc,sizeof(t_atom *));
                for (i=0;i<sksite->nc;i++){
                    fgets(buf,400,FP);
                    if(sscanf(buf,"%le",(sksite->conc)+i)!=1) COMPLAIN_EXIT;
                    fgets(buf,400,FP);
                    if(sscanf(buf,"%s",label)!=1) COMPLAIN_EXIT;
                    sksite->at[i]=getatom(atoms, label);
//                    printf("%s\t %s\n","sink name:",sitesk->at->label);
                };
            }
//                printf("%s\n","stop1");
//                site1->at=sksite->at;
//                site1->conc=sksite->conc;
//                site1->nc=sksite->nc;
//            } else {
//                sksite->at=site1->at;
//                sksite->conc=site1->conc;
//                sksite->nc=site1->nc;
//            };
    fclose(FP);
    LOGMSG(2,"Done with INPSK!");
    //////////////////////////////////////////
    
    // for (i=0;i<geom->tnat;i++){
    //         for (n=0;n<geom->site[i].nc;n++){
    //             printf("%s ", geom->site[i].at[n]->label);
    //         }
    //         printf("\n");
    // }
    //calculate scaling factor
    rawsr=0.0;
    //    LOGMSG(7,"test vars: %d %d ",geom->site,geom->site+geom->nb);
    for(site=(geom->site);site<((geom->site)+(geom->nb));site++){
        //    LOGMSG(7,"test vars2: %d %le %s",site,site->at[0]->wsr,site->at[0]->label);
        rawsr += pow(site->at[0]->wsr,3);
    }
    rawsr = pow((rawsr/geom->nb), (1.0/3.0));
    vol = fabs ((geom->trv[0][0]*geom->trv[1][1]-geom->trv[1][0]*geom->trv[0][1])*geom->ltrans[2]*geom->scx[0]*geom->scx[1]*geom->scx[2]);
    dawsr=pow((0.75*vol/(geom->nb*M_PI)),(1.0/3.0));
    geom->dimcf=rawsr/dawsr;
    // printf("%f\n", geom->dimcf);
    // for(site=geom->site;site<geom->site+geom->tnat;site++){
    //     printf("%s\n",site->at[0]->label);
    // }
    
    // for (i=0;i<=(geom->np+1);i++){
    // geom->plane[i+1]=(geom->plane[i])+geom->nb;
    // printf("%s\n",(geom->plane[i]+1)->at[0]->label);
    // };
    
    return(1);
}

///////////////////////////////
/// MSR: Having these duplicate readinput functions is a disaster, I'm sorry. But RSN started it....
/////////////// FOR TRILAYER INJECTION
int readinput3(t_scgeom *geom, t_atomset *atoms){
    FILE *FP;
    char buf[400], label[30];
    int i,n, extbf;
    t_scsite *site, *site0;
    t_sksite *sksite, *site1;
    double rawsr,dawsr,vol;
    
    atoms->ptr=NULL;
    atoms->num=0;
    // ADDD EOF DETECTION EVERYWHERE!!!!!!!!!!!!!!!!
    //     Read "inpge"
    LOGMSG(2,"Reading INPGE...");
    
    FP=openread("inpge");
    fgets(buf,400,FP);
    fgets(buf,400,FP);
    if(sscanf(buf,"%d %d %d",&(geom->np),&(geom->nb),&(geom->nmtr))!=3) COMPLAIN_EXIT;
    fgets(buf,400,FP);
    if (sscanf(buf,"%le",&(geom->cr))!=1) COMPLAIN_EXIT;
    geom->nat=geom->np*geom->nb;
    geom->tnat=(geom->np+2)*geom->nb;
    scanv3(buf,FP,geom->scx);
    scanv2(buf,FP,geom->trv[0]);
    scanv2(buf,FP,geom->trv[1]);
    scanv3(buf,FP,geom->ltrans);
    scanv3(buf,FP,geom->rtrans);
    geom->ltrans[0]=-geom->ltrans[0];
    geom->ltrans[1]=-geom->ltrans[1];
    geom->ltrans[2]=-geom->ltrans[2];
    
    geom->site=(t_scsite *)malloc(sizeof(t_scsite)*geom->tnat);
    geom->sitesk=(t_sksite *)malloc(sizeof(t_sksite)*1);
    geom->sitesk2=(t_sksite *)malloc(sizeof(t_sksite)*1);
    geom->plane=(t_scsite **)malloc(sizeof(t_scsite *)*(geom->np+2));
    
    for (i=0;i<(geom->tnat);i++){
        scanv3(buf,FP,geom->site[i].coord);
        geom->site[i].vmask=1.0;
        geom->site[i].rmask=1.0;
    };
    fclose(FP);
    geom->plane[0]=geom->site;
    for (i=0;i<=(geom->np);i++){
        geom->plane[i+1]=(geom->plane[i])+geom->nb;
    };
    
    LOGMSG(2,"Done with INPGE!");
    
    
    //     Read "inpch"
    LOGMSG(2,"Reading INPCH...");
    FP=openread("inpch");
    fgets(buf,400,FP);
    
    for (site=(geom->site+geom->nb);site<(geom->site+geom->nb+geom->nat);site++){
        LOGMSG(3,"Processing atom %d",site);
        fgets(buf,400,FP);
        if(sscanf(buf,"%d",&(site->nc))!=1) COMPLAIN_EXIT;
        LOGMSG(3,"Number of chemical elements on site: %d",site->nc);
        site->conc=(double *)calloc(site->nc,sizeof(double));
        site->at=(t_atom **)calloc(site->nc,sizeof(t_atom *));
        for (i=0;i<site->nc;i++){
            LOGMSG(3,"Loading %d-th concentration...",i);
            fgets(buf,400,FP);
            if(sscanf(buf,"%s",label)!=1) COMPLAIN_EXIT;
            LOGMSG(3,"Label = %s",label);
            fgets(buf,400,FP);
            if(sscanf(buf,"%le",(site->conc)+i)!=1) COMPLAIN_EXIT;
            LOGMSG(3,"Concentration = %le",site->conc[i]);
            site->at[i]=getatom(atoms, label);
        }
        
    };
    fclose(FP);
    LOGMSG(2,"Done with INPCH!");
    
    //     Read "inpbu"
    LOGMSG(2,"Loading INPBU...");
    FP=openread("inpbu");
    extbf=(fgetc(FP)=='+')?1:0;
    
    for (n=0;n<2;n++){
        fgets(buf,400,FP); //lead name line
        site0=NULL;
        for (site=geom->plane[(geom->np+1)*n];site<(geom->plane[(geom->np+1)*n]+geom->nb);site++){
            if ( (extbf==1) || (site0==NULL)){
                fgets(buf,400,FP);
                if(sscanf(buf,"%d",&(site->nc))!=1) COMPLAIN_EXIT;
                site->conc=(double *)calloc(site->nc,sizeof(double));
                site->at=(t_atom **)calloc(site->nc,sizeof(t_atom *));
                for (i=0;i<site->nc;i++){
                    fgets(buf,400,FP);
                    if(sscanf(buf,"%le",(site->conc)+i)!=1) COMPLAIN_EXIT;
                    fgets(buf,400,FP);
                    if(sscanf(buf,"%s",label)!=1) COMPLAIN_EXIT;
                    site->at[i]=getatom(atoms, label);
                };
                site0=site;
            } else {
                site->at=site0->at;
                site->conc=site0->conc;
                site->nc=site0->nc;
            };
        };
    };
    fclose(FP);
    LOGMSG(2,"Done with INPBU!");
    //////////////for spin injection//////////
    //     Read "inpsk"
    LOGMSG(2,"Loading INPSK...");
    FP=openread("inpsk");
        extbf=(fgetc(FP)=='+')?1:0;
        fgets(buf,400,FP); //sink name line
        site1=NULL;
   sksite= geom->sitesk;
            if ( (extbf==1) || (site1==NULL)){
                fgets(buf,400,FP);
                if(sscanf(buf,"%d",&(sksite->nc))!=1) COMPLAIN_EXIT;
                sksite->conc=(double *)calloc(sksite->nc,sizeof(double));
                sksite->at=(t_atom **)calloc(sksite->nc,sizeof(t_atom *));
                for (i=0;i<sksite->nc;i++){
                    fgets(buf,400,FP);
                    if(sscanf(buf,"%le",(sksite->conc)+i)!=1) COMPLAIN_EXIT;
                    fgets(buf,400,FP);
                    if(sscanf(buf,"%s",label)!=1) COMPLAIN_EXIT;
                    sksite->at[i]=getatom(atoms, label);
//                    printf("%s\t %s\n","sink name:",sitesk->at->label);
                };
            }
//                printf("%s\n","stop1");
//                site1->at=sksite->at;
//                site1->conc=sksite->conc;
//                site1->nc=sksite->nc;
//            } else {
//                sksite->at=site1->at;
//                sksite->conc=site1->conc;
//                sksite->nc=site1->nc;
//            };
    fclose(FP);
    LOGMSG(2,"Done with INPSK!");
    //////////////////////////////////////////
    //////////////for trilayer injection//////////
    //     Read "inpsk2"
    LOGMSG(2,"Loading INPSK2...");
    FP=openread("inpsk2");
        extbf=(fgetc(FP)=='+')?1:0;
        fgets(buf,400,FP); //sink name line
        site1=NULL;
   sksite= geom->sitesk2;
            if ( (extbf==1) || (site1==NULL)){
                fgets(buf,400,FP);
                if(sscanf(buf,"%d",&(sksite->nc))!=1) COMPLAIN_EXIT;
                sksite->conc=(double *)calloc(sksite->nc,sizeof(double));
                sksite->at=(t_atom **)calloc(sksite->nc,sizeof(t_atom *));
                for (i=0;i<sksite->nc;i++){
                    fgets(buf,400,FP);
                    if(sscanf(buf,"%le",(sksite->conc)+i)!=1) COMPLAIN_EXIT;
                    fgets(buf,400,FP);
                    if(sscanf(buf,"%s",label)!=1) COMPLAIN_EXIT;
                    sksite->at[i]=getatom(atoms, label);
//                    printf("%s\t %s\n","sink name:",sitesk->at->label);
                };
            }
//                printf("%s\n","stop1");
//                site1->at=sksite->at;
//                site1->conc=sksite->conc;
//                site1->nc=sksite->nc;
//            } else {
//                sksite->at=site1->at;
//                sksite->conc=site1->conc;
//                sksite->nc=site1->nc;
//            };
    fclose(FP);
    LOGMSG(2,"Done with INPSK2!");
    //////////////////////////////////////////
    
    // for (i=0;i<geom->tnat;i++){
    //         for (n=0;n<geom->site[i].nc;n++){
    //             printf("%s ", geom->site[i].at[n]->label);
    //         }
    //         printf("\n");
    // }
    //calculate scaling factor
    rawsr=0.0;
    //    LOGMSG(7,"test vars: %d %d ",geom->site,geom->site+geom->nb);
    for(site=(geom->site);site<((geom->site)+(geom->nb));site++){
        //    LOGMSG(7,"test vars2: %d %le %s",site,site->at[0]->wsr,site->at[0]->label);
        rawsr += pow(site->at[0]->wsr,3);
    }
    rawsr = pow((rawsr/geom->nb), (1.0/3.0));
    vol = fabs ((geom->trv[0][0]*geom->trv[1][1]-geom->trv[1][0]*geom->trv[0][1])*geom->ltrans[2]*geom->scx[0]*geom->scx[1]*geom->scx[2]);
    dawsr=pow((0.75*vol/(geom->nb*M_PI)),(1.0/3.0));
    geom->dimcf=rawsr/dawsr;
    // printf("%f\n", geom->dimcf);
    // for(site=geom->site;site<geom->site+geom->tnat;site++){
    //     printf("%s\n",site->at[0]->label);
    // }
    
    // for (i=0;i<=(geom->np+1);i++){
    // geom->plane[i+1]=(geom->plane[i])+geom->nb;
    // printf("%s\n",(geom->plane[i]+1)->at[0]->label);
    // };
    
    return(1);
}

///////////////////////////////

t_atom *getatom(t_atomset *atoms,char *label){
    t_atom *ptr;
    FILE *fl;
    int qq;
    for (ptr=atoms->ptr;ptr!=NULL; ptr=(t_atom *)(ptr->next)){
        if (!strcmp(ptr->label,label)) break;
    }
    LOGMSG(2,"Loading atomic file for %s...",label);
    
    // Add new atom
    if (ptr==NULL) {
        char buf[400];
        t_atom *nptr;
        sprintf(buf,"atoms/%s",label);
        fl=fopen(buf,"r");
        nptr=(t_atom *)malloc(sizeof(t_atom));        
        fgets(buf,400,fl);
        fgets(buf,400,fl);
        fgets(buf,400,fl);
        fgets(buf,400,fl);
        if(sscanf(buf,"%le %le",&(nptr->charge),&(nptr->wsr))!=2) COMPLAIN_EXIT;
        LOGMSG(9," checking input:%f %f | %s\n",nptr->charge,nptr->wsr,buf);
        nptr->mass=getmass(nptr->charge);
        // nptr->umsqr=0.0;
        fclose(fl);
        nptr->label=(char *)malloc(sizeof(char)*(strlen(label)+1));
        strcpy(nptr->label,label);
        atoms->num++;
        nptr->next=(void *)(atoms->ptr);
        atoms->ptr=nptr;
        ptr=nptr;
        // atoms->umsqr[i]=(Tdeb<0.)?0.:devmeansqr(TK,Tdeb,atoms->mass[i]);
    }
    LOGMSG(2,"Done loading atomic file for %s!",label);    
    return(ptr);    
}


double getmass(double charge){
    double ncharge;
    double mass;
    int icharge;

    // if (charge<(0.9999)) return -1.0;

    ncharge=floor(charge+1.0e-4);
    icharge=(int)(ncharge+1.0e-4);
    if (fabs(ncharge-charge)<1.0e-6){
        mass=atomic_mass[icharge];
    }else{
        mass=atomic_mass[icharge]+(atomic_mass[icharge+1]-atomic_mass[icharge])*(charge-ncharge);
    }
    return mass;
}


int cmpdist(const void*a,const void*b){
    if (((cstr*)a)->r<((cstr*)b)->r) return(-1);
    if (((cstr*)a)->r>((cstr*)b)->r) return(1);
    return(0);

}

double devmeansqr(double T0, double Tdeb){
    double res,qq,T,T2;
    int i;
    int maxe=550;

    T=T0/Tdeb;
    T2=T*T;
    res=0.;
    qq=1.;
    if (T>2e-3) {
        for (i=1;i<=maxe;i++){
            res-=1./(i*i*exp(i/T));        
        }
        res*=T2;
        qq=log(-1. + exp(1./T))*T;
    }
    res+=(T2*pow(M_PI,2)/6.-.75+qq);
    res=sqrt(res);
    return(res);
}


void prepvisdat(t_trgeom *geom,t_visdat *vis){
    int nuniq=0;
    int nuniqa=0;    
    t_atom *ptr;
    int na=geom->na;
    // int np=geom->np;
    // int nb=geom->nb;
    int j,k,l,i,ip,ia ,vn;
    char buf[400],*lbl;  
    int conn[MAX_NC],nc;  

    LOGMSG(5,"PrepVis1");

    vis->ulabels=(char **)malloc(sizeof(char *)*na);
    vis->ulats=(char **)malloc(sizeof(char *)*na);
    vis->x=(double *)malloc(sizeof(double)*na);   
    vis->y=(double *)malloc(sizeof(double)*na);   
    vis->z=(double *)malloc(sizeof(double)*na);   
    vis->wsr=(double *)malloc(sizeof(double)*na);   
    vis->ii1=(int *)malloc(sizeof(int)*na);   
    vis->ii2=(int *)malloc(sizeof(int)*na);   
    vis->plane=(int *)malloc(sizeof(int)*na);   
    vis->info=(char (*)[5])malloc(sizeof(char[5])*na);  
    vis->num=na;
    for(i=0;i<na;i++){
        vis->wsr[i]=(geom->site+i)->at->wsr*aunit;
        vis->x[i]=geom->rcoord[i][0]*aunit;
        vis->y[i]=geom->rcoord[i][1]*aunit;
        vis->z[i]=geom->rcoord[i][2]*aunit;
        vis->plane[i]=(geom->site+i)->pidx;
        lbl=(geom->site+i)->at->label;
        j=strcspn(lbl, "0123456789_");
        strcpy(buf,lbl);
        buf[j]=0;
        l=1;
        for (k=0;k<nuniq;k++){
            if (!strcmp(vis->ulabels[k],buf)){
                l=0;
                break;                  
            }
        }
        if (l){
            vis->ulabels[nuniq]=(char *)malloc(sizeof(char)*40);
            strcpy(vis->ulabels[nuniq],buf);
            vis->ii1[i]=nuniq;          
            nuniq++;      
        }
        else
        {
            vis->ii1[i]=k;
        };

        buf[2]=0;

        if (buf[0]=='E' | buf[0]=='e') {
            buf[1]=0;
            buf[0]='E';
            vis->ii2[i]=ES_ii2;  
        }
        else{

            l=1;
            for (k=0;k<nuniqa;k++){
                if (!strcmp(vis->ulats[k],buf)){
                    l=0;
                    break;                  
                }
            }
            if (l){
                vis->ulats[nuniqa]=(char *)malloc(sizeof(char)*3);
                strcpy(vis->ulats[nuniqa],buf);
                vis->ii2[i]=nuniqa;
                nuniqa++;                
            }
            else
            {
                vis->ii2[i]=k;  
            };
        };
        strcpy(vis->info[i],buf);
    }

    LOGMSG(5,"PrepVis2");
    vis->nuniqa=nuniqa;
    vis->nuniq=nuniq;

    vis->bonds1=(int *)malloc(sizeof(int)*vis->num*MAX_NUM_BOND);
    vis->bonds2=(int *)malloc(sizeof(int)*vis->num*MAX_NUM_BOND);
    vis->nbond=0;
    LOGMSG(5,"PrepVis5");
    {   
        register double rb,x1,y1,z1,x0,y0,z0,w0;
        vn=vis->num;
        for( i = 0; i < vn ; ++i)
        {
            nc=0;
            x0=vis->x[i];
            y0=vis->y[i];
            z0=vis->z[i];
            w0=vis->wsr[i];
            for( j = i+1; j < vn; ++j)
            {
                rb=vis->wsr[j]+w0;
                if (
                    ((z1=fabs(z0-vis->z[j]))>rb) ||
                    ((x1=fabs(x0-vis->x[j]))>rb) ||
                    ((y1=fabs(y0-vis->y[j]))>rb)
                             ) continue;
                if ((x1*x1+y1*y1+z1*z1)<rb*rb){
                    conn[nc]=j;
                    nc++;
                    if (nc==MAX_NC) break;
                }
            }

            if (nc>0){
                for( j = 0; j < nc; ++j)
                {
                    vis->bonds1[vis->nbond]=i;
                    vis->bonds2[vis->nbond]=conn[j];
                    vis->nbond++;
                }
            }
        }
    }
    LOGMSG(5,"PrepVis6");

}

void dumpvisdat(t_visdat * vis)
{
    FILE *fp2;
    int conn[150],nc,mnc=150,nbond;
    int *bonds1,*bonds2;
    size_t i,j;
    // cstr conn1[150];
    LOGMSG(2,"Writing structure visualization data...");
    fp2=fopen("vis.pdb","w");
    LOGMSG(7,"dumpvisdat: numat to write %d",vis->num);
    for(i = 0; i < vis->num; ++i)
    {
//        LOGMSG(7,"dumpvisdat: ith=%d",i);
        fprintf(fp2,"ATOM  %5i %-4.4s %3.3i  %4i    %8.3f%8.3f%8.3f%6.2f%6.2f %3i  %-4.4i%2.2s\n",(int)(i+1),vis->info[i],vis->ii1[i],vis->ii2[i],vis->x[i],vis->y[i],vis->z[i],1.0,0.0,vis->ii1[i],vis->plane[i],vis->info[i]);
//        LOGMSG(7,"dumpvisdat: ith=%d - done",i);
            // fprintf(fp2,"HETATM%5i %-4.4s %3.3i  %4i    %8.3f%8.3f%8.3f%6.2f%6.2f %3i  %-3.3i %2.2s\n",i,visdat[i].info,visdat[i].ii1,visdat[i].ii2,visdat[i].x,visdat[i].y,visdat[i].z,0.0,0.0,visdat[i].ii1,visdat[i].plane,visdat[i].info);
    };
    LOGMSG(7,"dumpvisdat: loop done");
    fprintf(fp2,"END\n");        
    fclose(fp2);
    
    LOGMSG(5,"pdb is done");


    fp2=fopen("vis.psf","w");
    LOGMSG(7,"dumpvisdat-6");
    fprintf(fp2,"PSF \n\n");        
    LOGMSG(7,"dumpvisdat-7");
    LOGMSG(7,"dumpvisdat-8");
    fprintf(fp2,"%8i !NATOM\n",vis->num);
    for(i = 0; i < vis->num; ++i)
    {
        fprintf(fp2,"%8i %-4.4s %-4i %3.3i  %-4.4s %-4.4s %12.6lf %12.6lf\n",(int)(i+1),vis->info[i],vis->ii2[i],vis->ii1[i],vis->info[i],"AAA",1.0,1.0);
    };


    LOGMSG(7,"dumpvisdat-9");

    fprintf(fp2,"\n");

    fprintf(fp2,"%8i !NBOND: bonds\n",vis->nbond);        

    for( i = 0; i < vis->nbond; ++i)
    {
        fprintf(fp2,"%8i%8i",vis->bonds1[i]+1,vis->bonds2[i]+1);
        if (((i+1)%4)==0) fprintf(fp2,"\n");
    }
    fprintf(fp2,"\n");
    fclose(fp2);
    
    LOGMSG(2,"Done!");


}

void writeatomlist(t_atomset *atoms){
    FILE *fl;
    t_atom *ptr;
    fl=fopen("atomlist","w");
    fprintf(fl,"-------=========== ATOM DEFINITIONS ===========-------\n");
    fprintf(fl,"#FORMAT OF LINES:  LABEL    FILENAME\n");
    fprintf(fl,"%5d  %29s    # mass (a.u)\n",atoms->num,"# NUMBER OF ATOMS");

    for (ptr=atoms->ptr;ptr!=NULL; ptr=(t_atom *)(ptr->next)){
        fprintf(fl,atfmt,ptr->label,ptr->label,ptr->mass);
    }
    fprintf(fl,"-------=========== END OF ATOM FILE ===========-------");
    fclose(fl);
}
inline FILE *openread(char *nm) {
    FILE *a=fopen(nm,"r");
    if (a==NULL) { 
        printf("Did you put \"%s\" file in current dir? :)\n",nm);
        exit(1); 
    }
    return(a);
}

void writegeom(t_trgeom *geom, char *name){
    int i;
    t_trsite *site;
    FILE *fp=fopen(name,"w");
    LOGMSG(5,"writegeom1");
    fprintf(fp,"-------===========  GEOMETRY  FILE ===========-------\n");
    fprintf(fp," %5d   %5d   %5d %45s\n",geom->na,geom->nmtr,geom->tp,"# NP,  NMTR,TR_PERP");
    fprintf(fp," %#15.10g %40s\n",geom->cr,"# CUTRAT");
    fprintf(fp," %5d   %5d    %40s\n",geom->sc[0],geom->sc[1],"# SC_SIZE");
    fpvec3(fp,geom->scx,"# SCALING FACTORS");
    fpvec2(fp,geom->trv[0],"# 1ST TRANSL. VECTOR");
    fpvec2(fp,geom->trv[1],"# 2ND TRANSL. VECTOR");
    fpvec3(fp,geom->trans,"# PERP.TR. VECTOR");
    LOGMSG(5,"writegeom3, NA=%d",geom->na);
    for (site=geom->site;site<(geom->site+geom->na);site++){
        LOGMSG(9,"writegeom, site=%d",site-geom->site);
        LOGMSG(9,"writegeom, coord=(%le,%le,%le)",site->coord[0],site->coord[1],site->coord[2]);        
        LOGMSG(9,"writegeom, label=%s",site->at->label);
        fprintf(fp,fmt3s,site->coord[0],site->coord[1],site->coord[2],site->at->label);
    }
    fprintf(fp,"-------=========== END OF GEOMETRY ===========-------");
    fclose(fp);
}
//--------------------
//void modifygeom(t_trgeom *geom,nlayer){
//    int i;
//    t_trsite *site;
//    for (site=geom->site;site<(geom->site+geom->na);site++){
//        if (site->coord[1]>nlayer) {
//            site->at->label=label2
//        }
//    }
//    fprintf(fp,"-------=========== END OF GEOMETRY ===========-------");
//    fclose(fp);
//}
//------------------

inline void getcoordz(double *pos,double *coord,double *tr, int l) {
    pos[0]=coord[0]+l*tr[0];
    pos[1]=coord[1]+l*tr[1];
    pos[2]=coord[2]+l*tr[2];
}

inline double randn0()
{
    double x1, x2, w, y1;
    do {
        x1 = 2.0 * drand48() - 1.0;
        x2 = 2.0 * drand48() - 1.0;
        w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    // y2 = x2 * w;
    return(y1);
}


inline double randn()
{
    return(ltqnorm(drand48()));
}

void writerotcf(t_trgeom *geom){
    FILE *fp;
    int i;
    fp=fopen("mask","w");
    for (i=0;i<geom->na;i++) fprintf(fp,"%f\n",(double)(geom->site[i].mask));
    fclose(fp);

}

void writerotmask(t_trgeom *geom){
    FILE *fp;
    int i;
    fp=fopen("rot_mask","w");
    fprintf(fp, "#theta phi\n");
    for (i=0;i<geom->na;i++) fprintf(fp,"%#15.10g  %#15.10g\n",geom->site[i].mdir[0],geom->site[i].mdir[1]);
    fclose(fp);

}

void initrnd(){
    struct timeval tmv;
    short unsigned int dval[3]; 
    long int sval;

    gettimeofday(&tmv,(struct timezone *)0);
    sval=(time(NULL)*1e6+tmv.tv_usec);
    srand48(sval);
    dval[0]=lrand48();
    dval[1]=lrand48();
    dval[2]=lrand48();
    seed48(dval);
    gettimeofday(&tmv,(struct timezone *)0);  
    dval[0]=lrand48()+tmv.tv_usec;
    dval[1]=lrand48();
    dval[2]=lrand48();
    seed48(dval);
}


int cholesky(int N,int *Ap,int *Ai,double *Ax,int **Lp,int **Li,double **Lx, double **D){
    int Flag[N],Parent[N], Lnz[N],d,Pattern[N];
    double Y[N];

    (*Lp)=(int *)realloc((*Lp),sizeof(int)*(N+1));

    // Flag=(int *)malloc(sizeof(int)*(N));
    // Parent=(int *)malloc(sizeof(int)*(N));
    // Lnz=(int *)malloc(sizeof(int)*(N));
    // Y=(double *)malloc(sizeof(double)*(N));
    // Pattern=(int *)malloc(sizeof(int)*(N));


    ldl_symbolic (N, Ap, Ai, *Lp, Parent, Lnz, Flag, NULL, NULL) ;
    (*Li)=(int *)realloc((*Li),sizeof(int)*((*Lp)[N]));
    (*Lx)=(double *)realloc((*Lx),sizeof(double)*((*Lp)[N]));
    (*D)=(double *)realloc((*D),sizeof(double)*(N));

    d = ldl_numeric (N, Ap, Ai, Ax, *Lp, Parent, Lnz, *Li, *Lx, *D, Y, Pattern, Flag, NULL, NULL) ;
    // printf ("Nonzeros in L, excluding diagonal: %d , %d\n", (*Lp)[N],d) ;
    if (d!=N) {
        printf("cholesky error!\n");
        exit(0);
    }
    return(d);
}

int get_cf_dists(t_trgeom *geom, int** ApI,int **AiI, double **res, int **idxI){
    int i,j,n,nm,nz,mnz;
    int N=geom->na;
    int *idx;
    int *Ap,*Ai;
    double *dx,*dy,*dz,*dd,swsr[geom->na];
    double rcoord[geom->na][3], c1[3],c2[3], *cp, ct[2];
    double cr,trv[2][2],dist,dz2,wsr1,mdist,dx0,dy0,dz0, dist0,dx1,dy1;
    int ix,iy,nmtr;

    nmtr=geom->nmtr;
    cr=geom->cr;    
    idx=(int *)malloc(sizeof(int)*geom->na);
    for (i=0;i<2;++i) for (j=0;j<2;++j) trv[i][j]=geom->trv[i][j]*geom->scx[j]*geom->dimcf;
    n=0;
    for (i=0;i<N;++i){
        if (geom->site[i].dmask!=0){
            idx[n]=i;
            swsr[n]=sqrt(cr*geom->site[i].at->wsr);
            memcpy(rcoord[n],geom->rcoord[i],sizeof(double)*3);
            n++;
        }
    }
    nm=n;
    mnz=nm*30;
    Ap=(int *)malloc(sizeof(int)*(nm+1));
    Ai=(int *)malloc(sizeof(int)*mnz);
    for (i=0;i<5;++i) res[i]=(double *)malloc(sizeof(double)*mnz);
    // printf("%d,%d\n",nm,N);

    Ap[0]=0;
    nz=0;

// Lower/Upper part
#define CORR_LOWER_PART 1

    for (i=0;i<nm;++i){            
        Ap[i]=nz;
#if CORR_LOWER_PART != 1
        Ai[nz]=i;
        res[0][nz]=res[1][nz]=res[2][nz]=res[3][nz]=0.0;
        res[4][nz]=geom->site[idx[i]].at->mass;
        nz++;
#endif        
        wsr1=swsr[i];
        memcpy(c1,rcoord[i],sizeof(double)*3);
#if CORR_LOWER_PART != 1
        for (j=i+1;j<nm;++j){
#else
            for (j=0;j<i;++j){
#endif
                cp=rcoord[j];
                c2[2]=cp[2];
                dz0=c1[2]-c2[2];
                dz2=pow(dz0,2);
                mdist=wsr1*swsr[j];
                dist0=2*mdist;
                for (ix=-nmtr;ix<=nmtr;++ix){
                    ct[0]=cp[0]+ix*trv[0][0];
                    ct[1]=cp[1]+ix*trv[0][1];
                    for (iy=-nmtr;iy<=nmtr;++iy){
                        c2[0]=ct[0]+iy*trv[1][0];
                        c2[1]=ct[1]+iy*trv[1][1];
                        dx1=c1[0]-c2[0];
                        dy1=c1[1]-c2[1];
                        dist=sqrt(pow(dx1,2)+pow(dy1,2)+dz2);
                        if ((dist<mdist) && (dist<dist0)){
                            dist0=dist;
                            dx0=dx1;
                            dy0=dy1;
                        }
                    }                
                }
                if (dist0<mdist){
                // if (i==1){
                //  // printf("%d, %d : %f : %f, %f, %f : %f, %f, %f : %f, %f, %f\n",i,j,dist0,c1[0],c1[1],c1[2],c2[0],c2[1],c2[2],dx0,dy0,dz0);
                //  printf("%f : %f, %f, %f : %f, %f, %f : %f, %f, %f\n",dist0,c1[0],c1[1],c1[2],c2[0],c2[1],c2[2],dx0,dy0,dz0);
                // }
                    Ai[nz]=j;
                    res[0][nz]=dist0;
                    res[1][nz]=dx0;
                    res[2][nz]=dy0;
                    res[3][nz]=dz0;
                    res[4][nz]=(geom->site[idx[i]].at->mass*geom->site[idx[j]].at->mass)/(geom->site[idx[i]].at->mass+geom->site[idx[j]].at->mass);
                    nz++;                
                }
            }   
#if CORR_LOWER_PART == 1
            Ai[nz]=i;
            res[0][nz]=res[1][nz]=res[2][nz]=res[3][nz]=0.0;
            res[4][nz]=geom->site[idx[i]].at->mass;            
            nz++;
#endif

        };
        Ap[nm]=nz;

        (*AiI)=Ai;
        (*ApI)=Ap;
        (*idxI)=idx;

        return (nm);
    }

    void dump_sparse(int nm,int *Ap, int *Ai, double *Ax, char *name){
        FILE *fp;
        int i,j;
        fp=fopen(name,"w");
        for (i=0;i<nm;i++){
            for (j=Ap[i];j<Ap[i+1];j++){
                fprintf(fp,"%d  %d  %#22.18e %#22.18e\n",i+1, Ai[j]+1,Ax[j],0.);
            }        
        }
        fprintf(fp,"%d  %d  %#22.18e %#22.18e\n",nm, nm,0.,0.);
        fclose(fp);
    }

    void dump_vector(int nm,double *v, char *name){
        FILE *fp;
        int i;
        fp=fopen(name,"w");
        for (i=0;i<nm;i++){
            fprintf(fp,"%d  %d  %#22.18e %#22.18e\n",i+1, 1,v[i],0.);
        };
        fclose(fp);
    }

    void fillcorf(int nnz, double *corf,double *dst, double *dx,double wdeb, double b){
        int i;
        double d0,d1,d2,dov,q;
        q=.5/(b*b);
        for (i=0;i<nnz;++i){
            d0=dst[i];
            if (d0!=0.0){
                dov=d0*wdeb; ///pow(3.,1./3.); // Remove it later!!!
                d1=dx[i];
                corf[i]=q*(1-cos(dov))/pow(dov,2)*(d1*d1/(d0*d0));
            }
            else{
                corf[i]=1.;
            }
        }
    }
    void cholmult(int n,int *Lp, int* Li,double *Lx,double *Dg,double *disp0,double *disp1){
        int i,j;
        double x;
    // memset(disp1,0,sizeof(double)*n);   
        for (i=0;i<n;++i) disp0[i]*=sqrt(abs(Dg[i]));
        memcpy(disp1,disp0,sizeof(double)*n);   

        for (i=0;i<n;++i){
            x=disp0[i];
            for(j=Lp[i];j<Lp[i+1];++j){
                disp1[Li[j]]+=Lx[j]*x;
            }
        }; 

    }
    
    void rotspheric (double *basei, double *anglei,double *nang) {
        double base[2],angle[2];
        int i;
        base[0]=-M_PI*basei[0];
        base[1]=-M_PI*basei[1];
        angle[0]=M_PI*anglei[0];
        angle[1]=M_PI*anglei[1];
        nang[0] = acos (cos(angle[0])*cos(base[0])+cos(angle[1])*sin(angle[0])*sin(base[0])) / M_PI;
        nang[1] = atan2 (cos(base[1])*sin(angle[0])*sin(angle[1])+(
        -(cos(angle[1])*cos(base[0])*sin(angle[0]))+cos(angle[0])*sin(base[0]))*sin(base[1]),
        cos(angle[1])*cos(base[0])*cos(base[1])*sin(angle[0])-cos(angle[0])*cos(base[1])*sin(base[0])+
        sin(angle[0])*sin(angle[1])*sin(base[1])) / M_PI;
    }

    // void geteasymagnons(t_trgeom *geom,double Mtheta, double Mphi,double Mdev){
    //     int i,j;
    //     double basei[2],devang[2],nang[2];
    // 
    //     for (i=0;i<geom->na;++i){
    //         basei[0]=nang[0]=Mtheta;
    //         basei[1]=nang[1]=Mphi;            
    //         if (geom->site[i].mask!=0){
    //             devang[0]=Mdev*fabs(randn());
    //             devang[1]=drand48()*2.0-1.0;
    //             rotspheric(basei,devang,nang);
    //         };  
    //         geom->site[i].mdir[0]=nang[0];      
    //         geom->site[i].mdir[1]=nang[1];
    //     };
    // 
    // }
    void geteasymagnons(t_trgeom *geom,double Mdev){
        int i,j=0;
        double nang[2],q=0,z=0,y=0,x=0;
        geom->mvals[0]=0.0;
        geom->mvals[1]=0.0;
        geom->mvals[2]=0.0;
        for (i=0;i<geom->na;++i){
            nang[0]=0.0;
            nang[1]=0.0;            
            if (geom->site[i].dmask!=0){
                j++;
                nang[0]=Mdev*fabs(randn())*geom->site[i].rmask;
                q=MAX(q,nang[0]);
                nang[1]=(drand48()*2.0-1.0)*geom->site[i].rmask;
                geom->mvals[0]+=(nang[0]*nang[0]);
                geom->mvals[1]+=nang[1];
                z+=cos(nang[0]*M_PI);
                x+=sin(nang[0]*M_PI)*cos(nang[1]*M_PI);
                y+=sin(nang[0]*M_PI)*sin(nang[1]*M_PI);
            };  
            geom->site[i].mdir[0]=nang[0];      
            geom->site[i].mdir[1]=nang[1];
        };
        geom->mvals[0]=sqrt(geom->mvals[0]/((double)j));
        geom->mvals[1]/=j;
        geom->mvals[2]=q;

        q=sqrt(z*z+x*x+y*y);
        z/=q;
        x/=q;
        y/=q;
/*         printf("%f,%f,%f\n",x,y,z);*/
        geom->mvals[3]=acos(z)/M_PI;
        geom->mvals[4]=atan2(y,x)/M_PI;
        geom->mvals[5]=q/j;
    }

    double getMdev(double TTc){
        double fitpar[10];
        double a1,a2,ret;
        int i;
        FILE *fid;
        
        fid=fopen("magfit","r");
        if (fid<=0) {
            printf("input file with fit params of magnetization is missing!\n");
            exit(-1);
        }
        for(i=0;i<10;i++) {
            fscanf(fid,"%le\n",fitpar+i);
            // printf("%le\n",fitpar[i]);
        }
        a1=0.0;
        a2=0.0;
        for (i=0;i<5;i++){
            a1+=(fitpar[i]*pow(TTc,i+1)-fitpar[i]);
            a2+=(fitpar[i+5]*pow(TTc,i+1)-fitpar[i]);
        }
        ret=sqrt(-log(a1/a2)*2.0)/M_PI;
        return(ret);
    }
        
    void geteasyvibrations(t_trgeom *geom,double Tdeb, double Tk,double mdef){
        int i,j,k;
        double cf0,cf1,x,qq,rms;
        cf0=sqrt(3.*debcf/Tdeb)*devmeansqr(Tk, Tdeb)*1e9/bohrrad/geom->dimcf*mdef;
        rms=0;
        k=0;
        for (i=0;i<geom->na;++i){
            if (geom->site[i].dmask!=0){
                cf1=cf0/sqrt(geom->site[i].at->mass)*geom->site[i].vmask;
                for (j=0;j<3;++j){
                    x=cf1*randn();
                    geom->site[i].coord[j]+=x/geom->scx[j];
                    geom->rcoord[i][j]+=x;                
                };
//                printf("%10i\t, %lg\t, %lg\n",i,geom->site[i].at->mass, 3.*cf1*cf1);
                rms+=3.*cf1*cf1;
                k++;
            };
        };
        rms=sqrt(rms/k);
        printf("Dvib=%lg, %lg, %lg\n",rms,rms*geom->dimcf,mdef);
    }
//-----------------------
    void getrms(t_trgeom *geom){
        FILE *fp;
        char ch,buf[400];
        char *str1, *str2;
        int i,k,length;
        double irms;
        length=0;
        fp=fopen("rms", "r");
        if(fp !=NULL){
            fscanf(fp, "%s\n",buf);
            length = atoi(buf);
            if(length<0)COMPLAIN_EXIT;
            printf("Following vibrations are applied from the file rms:\n");
            for (k=0;k<length;++k){
                fscanf(fp, "%s\n",buf);
                str1=buf;
                if(str1 !=NULL){
                    str2=strsep(&str1,",");
                    irms=atof(str1);
                    printf("%s:\t %lg\n",str2,irms);
//-----------------------------------------------
                    for (i=0;i<geom->na;++i){
                        if (geom->site[i].dmask!=0){
                                if(!strcmp(geom->site[i].at->label,str2)){
                                    geom->site[i].at->rms=irms/(geom->dimcf*sqrt(3));
                                };
                        };
                    };
                   //
                };
            };
            fclose(fp);
        }else {
            printf( "Could not open file rms\n" ) ;
            exit(0);
        };
        }
    void getsepvibrations(t_trgeom *geom){
        int i,j;
        double cf,x;
        cf=0;
        getrms(geom);
        for (i=0;i<geom->na;++i){
//            printf("%s\t,%i\n","stop1",i);
                if (geom->site[i].dmask!=0){
                            cf=geom->site[i].at->rms*geom->site[i].vmask;
                            for (j=0;j<3;++j){
                                x=cf*randn();
                                geom->site[i].coord[j]+=x/geom->scx[j];
                                geom->rcoord[i][j]+=x;
                        };
                    };
            };
        }
    
// -----------------

    void getvibrations(t_trgeom *geom,double Tdeb, double Tk,double mdef){
        int *Lp=NULL,*Li=NULL, *idx=NULL;
        double *Lx=NULL,*Dg=NULL;
        int *Ap=NULL,*Ai=NULL;
        double *dd=NULL,*dc[5];
        int nm,i,j,k,nnz;
        double *corf;
        double umbase,wdeb,avol,dst,dx,dov;
        double *disp0,*disp1[3];
        double x;
        umbase=devmeansqr(Tk, Tdeb);

        nm=get_cf_dists(geom, &Ap,&Ai,dc,&idx);
        nnz=Ap[nm];

        avol=0.0;
        for (i=0;i<nm;++i){
            avol+=pow(geom->site[idx[i]].at->wsr,3);
        }
        avol*=4./3.*M_PI/nm;

    // printf("V0=%le, %le\n",avol, umbase);
        wdeb=pow(6*M_PI*M_PI/avol,1.0/3.0);
        printf("Debye wave vector : %lf\n",wdeb);
    // printf("cf=%le\n",devmeansqr(0., 375));

        corf=(double *)malloc(sizeof(double)*nnz);
        disp0=(double *)malloc(sizeof(double)*nm);

        for (i=0;i<3;++i){
            fillcorf(nnz,corf,dc[0],dc[i+1],wdeb,umbase);
            cholesky(nm,Ap,Ai,corf,&Lp,&Li,&Lx,&Dg);
            for (j=0;j<nm;++j) disp0[j]=randn();
            disp1[i]=(double *)malloc(sizeof(double)*nm);        
            cholmult(nm,Lp,Li,Lx,Dg,disp0,disp1[i]);
        }
        dov=sqrt(3*debcf/Tdeb)*umbase*1e9/bohrrad/geom->dimcf;
        printf("Thermo disp. fact.: %lf\n",dov);
        
    // printf("%le , %le , %le, %le\n",dov,debcf,umbase,bohrrad);
        for (j=0;j<nm;++j){
            k=idx[j];
            for (i=0;i<3;++i){
                x=disp1[i][j]*(dov/sqrt(geom->site[k].at->mass));
                geom->site[k].coord[i]+=x/geom->scx[i]*geom->site[k].vmask;
                geom->rcoord[k][i]+=x;
            }
        }

    // dump_sparse(nm,Ap,Ai,corf,"cor.dat");


    // dump_sparse(nm,Lp,Li,Lx,"A.dat");
    // dump_vector(nm,Dg,"D.dat");

    // dump_sparse(nm,Ap,Ai,dc[0],"dd.dat");
    // dump_sparse(nm,Ap,Ai,dc[1],"dx.dat");
    // dump_sparse(nm,Ap,Ai,dc[2],"dy.dat");
    // dump_sparse(nm,Ap,Ai,dc[3],"dz.dat");

    // dump_vector(nm,disp0,"v0.dat");


    // dump_vector(nm,disp1,"v1.dat");



     // for (i=0;i<nm;i++){
     //     printf("%d\t%d\n",i, Ap[i+1]-Ap[i]);
     //     
     // }
     // printf("%d\n",nz);
    }

#define MAXLOGSIZE 2048
    void log_print(const char *file, int line, char *message, ...){
        va_list args; 
        char buf[MAXLOGSIZE];
        
        va_start(args, message);
        vsnprintf(buf,MAXLOGSIZE-1,message,args);
        printf("%s:%d : %s\n",file, line, buf);
        va_end(args);
    }
