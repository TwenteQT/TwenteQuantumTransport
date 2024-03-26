#include "libgen.h"

typedef struct{
    int ord, nx,ny,N,lbnl,rbnl;
    int vis,mask, corr;
    int reps,repn,nreps,skip_dis;
    double scl;
    double dimcf,Tdeb,TK,mdef,Tc,Tm;    
    double Mtheta,Mphi,Mdev;
    int vibmask,rotmask,plord,vib,pin;
    int nlayer, nlayer2;
    } t_opt;


// int genscdist(t_scsite *plane,t_trsite  *isite,double (*tr)[2],double *ofs,int nx,int ny,int nb, int ord, int pidx);
// void maketrgeom(t_scgeom *inge, t_trgeom *geom, t_trgeom *lg, t_trgeom *rg, int nx,int ny,int nlb,int nrb, int ord);
void extrafields(t_trgeom *geom);
static void usage(void);
void readopts(int argc, char *argv[], t_opt *par);

//inline int genscdist(t_scsite *plane,t_trsite  *isite,double (*tr)[2],double *ofs,int nx,int ny,int nb, int ord, int pidx,int **kind, int ko){
int genscdist(t_scsite *plane,t_trsite  *isite,double (*tr)[2],double *ofs,int nx,int ny,int nb, int ord, int pidx,int **kind, int ko){
    int j,at,nn,l,p;
    double pr;
    t_trsite  *site;

    site=isite;
    if (ord!=0)     
    for (at=0;at<nb;at++){
        for (j=0;j<ny;j++){
            for (l=0;l<nx;l++)
            {

                // printf("%d, %s, %s\n",at,(plane+at)->at[kind[at][nx*j+l]]->label,(plane+at)->at[0]->label);
                site->coord[0]=(plane+at)->coord[0]+l*tr[0][0]+j*tr[1][0]+ofs[0];
                site->coord[1]=(plane+at)->coord[1]+l*tr[0][1]+j*tr[1][1]+ofs[1];
                site->coord[2]=(plane+at)->coord[2]+ofs[2];
                site->at=(plane+at)->at[kind[at][nx*j+l+ko]];
                site->pidx=pidx;
                site->mask=1;
                site->dmask=1;
                site->rmask=(plane+at)->rmask;
                site->vmask=(plane+at)->vmask;
                site->ply=j;
                site++;
            };
        };
    }
    else
    for (j=0;j<ny;j++){
        for (l=0;l<nx;l++){
            for (at=0;at<nb;at++)
            {
                site->coord[0]=(plane+at)->coord[0]+l*tr[0][0]+j*tr[1][0]+ofs[0];
                site->coord[1]=(plane+at)->coord[1]+l*tr[0][1]+j*tr[1][1]+ofs[1];
                site->coord[2]=(plane+at)->coord[2]+ofs[2];
                site->at=(plane+at)->at[kind[at][nx*j+l+ko]];
                site->pidx=pidx;
                site->mask=1;
                site->dmask=1;
                site->rmask=(plane+at)->rmask;
                site->vmask=(plane+at)->vmask;
                site->ply=j;
                site++;
            };
        };
    };


    return(nx*ny*nb);    
}
//--------------------------------------------
//inline int genscsepdist(t_scsite *plane,t_trsite  *isite,double (*tr)[2],double *ofs,int nx,int ny,int nb, int ord, int pidx,int **kind, int ko,int nlayer,char *sink){
//    int j,at,nn,l,p;
//    double pr;
//    t_trsite  *site;
//
//    site=isite;
//    if (ord!=0)
//    for (at=0;at<nb;at++){
//        for (j=0;j<ny;j++){
//            for (l=0;l<nx;l++)
//            {
//
//                // printf("%d, %s, %s\n",at,(plane+at)->at[kind[at][nx*j+l]]->label,(plane+at)->at[0]->label);
//                site->coord[0]=(plane+at)->coord[0]+l*tr[0][0]+j*tr[1][0]+ofs[0];
//                site->coord[1]=(plane+at)->coord[1]+l*tr[0][1]+j*tr[1][1]+ofs[1];
//                site->coord[2]=(plane+at)->coord[2]+ofs[2];
//                site->at=(plane+at)->at[kind[at][nx*j+l+ko]];
//                site->pidx=pidx;
//                site->mask=1;
//                site->dmask=1;
//                site->rmask=(plane+at)->rmask;
//                site->vmask=(plane+at)->vmask;
//                if(j > nlayer){
//                    site->at->label=sink;
//                }
//                site++;
//            };
//        };
//    }
//    else
//    for (j=0;j<ny;j++){
//        for (l=0;l<nx;l++){
//            for (at=0;at<nb;at++)
//            {
//                site->coord[0]=(plane+at)->coord[0]+l*tr[0][0]+j*tr[1][0]+ofs[0];
//                site->coord[1]=(plane+at)->coord[1]+l*tr[0][1]+j*tr[1][1]+ofs[1];
//                site->coord[2]=(plane+at)->coord[2]+ofs[2];
//                site->at=(plane+at)->at[kind[at][nx*j+l+ko]];
//                site->pidx=pidx;
//                site->mask=1;
//                site->dmask=1;
//                site->rmask=(plane+at)->rmask;
//                site->vmask=(plane+at)->vmask;
//                if(j > nlayer){
//                    //printf("At:%i\t, before:%s\n",j,site->at->label);
//                    site->at->label=sink;
//                    //printf("after:%s\t %s\n",site->at->label,buff);
//                }
//                site++;
//            };
//        };
//    };
//
//
//    return(nx*ny*nb);
//}
void genpingeom(t_trgeom *geom, t_scgeom *geomsc,int nlayer,int m){
    int i;
    t_trsite *site;
    for (site=geom->site;site<(geom->site+geom->na);site++){
        if(site->ply > nlayer){
          //  printf("%s\t,%s\n","Atom label :",skgeom->site->at->label);
            site->at=*(geomsc->sitesk->at);
        };
//        printf("layercheck :%i\t,%s\n",site->ply,site->label );
    };
};
void genpingeom_trilayer(t_trgeom *geom, t_scgeom *geomsc,int nlayer,int nlayer2,int m){
    int i;
    t_trsite *site;
    for (site=geom->site;site<(geom->site+geom->na);site++){
        if(site->ply > nlayer){
          //  printf("%s\t,%s\n","Atom label :",skgeom->site->at->label);
            site->at=*(geomsc->sitesk->at);
            if (site->ply > nlayer2) {
                site->at=*(geomsc->sitesk2->at);
            }
        };
//        printf("layercheck :%i\t,%s\n",site->ply,site->label );
    };
};

//void genpingeom(t_trgeom *geom,int nlayer,char *sink,int m){
//    int i;
//    t_trsite *site;
//    if(m){
//        for (site=geom->site;site<(geom->site+geom->na);site++){
//            if(site->ply > nlayer){
//                site->label=sink;
//            }
//            else{
//                site->label=site->at->label;
//            }
//            //        printf("layercheck :%i\t,%s\n",site->ply,site->label );
//        };
//    }
//    else
//    {
//        for (site=geom->site;site<(geom->site+geom->na);site++){
//            site->label=site->at->label;
//        }
//
//    };
//};
    //
//    for (i=0;i<geom->na;++i){
//        if(geom ->site[i].ply > nlayer){
//            geom->site[i].at->label="N_B";            //sink;
//  //              printf("layer : %i,%s\n,",geom ->site[i].ply,sink );
//        };
//    };
//}

//---------------------------------------------
//inline void genchemdist(t_scsite *plane, int N, int nb, int **kind, int ko)
void genchemdist(t_scsite *plane, int N, int nb, int **kind, int ko)
{   
    int j,l,at,nn,p;
    double pr;
    LOGMSG(5,"genchemdist start for N=%d, nb=%d, ko=%d",N,nb,ko);
        for (j=0;j<nb;++j){
            memset(kind[j]+ko,0,sizeof(int)*N); 
            LOGMSG(8,"genchemdist2: j=%d, nc=%d",j,((plane+j)->nc));
            for (at=1;at<((plane+j)->nc);at++) {
                nn=round(((plane+j)->conc[at])*N);
//                printf("%d,%d,%d\n",nn,N,at);
                l=0;
                while (l<nn) {
                    pr=drand48()-1e-16;
                    p=floor(N*pr);
                    if (kind[j][p+ko]==0){
                        l++;
                        kind[j][p+ko]=at;
                    };
                };
            };
        };
    LOGMSG(5,"genchemdist done");
};


void maketrgeom(t_scgeom *inge, t_trgeom *geom, t_trgeom *lg, t_trgeom *rg, int nx,int ny,int nlb,int nrb, int ord,int reps,int repn,int nreps,int nskip,int plord){

    int na,N,nn,i,j,k,pl,l;
    double ofs[3]={0,0,0};
    double ofs0[3]={0,0,0};
    double trans[3]={0,0,0};
    t_trsite *site;
    int repe;
    int **kind;
    // int nreps;
    
    // nreps=nr+1;
    // repn=MAX(repn0,1);
    repe=reps+repn;
    N=nx*ny;
    LOGMSG(5,"maketransgen call, N=%d, repe=%d, inge->nb=%d, nreps=%d,repn=%d, reps=%d",N,repe,inge->nb,nreps,repn,reps);
    
    kind=malloc(sizeof(int *)*inge->nb);
    
    LOGMSG(5,"");
    
    for(i = 0; i < inge->nb; ++i)
    {
        kind[i]=malloc(sizeof(int)*N*nreps*MAX(repn,1));
    }
    
    // printf("%d\n",inge->nb*repn*(nreps-1));
    na=(inge->nat+inge->nb*nlb+inge->nb*nrb+inge->nb*repn*(nreps-1))*N;
    // printf("%d, %d, %d, %d, %d, %d\n",reps,repe,repn,nreps,na);

    LOGMSG(5,"maketransgen call, NA=%d",na);
    geom->site=(t_trsite  *)malloc(sizeof(t_trsite )*na);    
    
    nn=0;
    pl=0;
    for (i=0;i<nlb;i++){
        ofs[0]+=inge->ltrans[0];ofs[1]+=inge->ltrans[1];ofs[2]+=inge->ltrans[2];
        genchemdist(inge->plane[0], N, inge->nb, kind,0);
        nn+=genscdist(inge->plane[0],geom->site+nn,inge->trv,ofs,nx,ny,inge->nb, ord,pl,kind,0);
        pl++;
    }        

    LOGMSG(5,"maketransgen plane idx after left replication %d",pl);
    
    for (i=0;i<(reps);i++){
        genchemdist(inge->plane[i+1], N, inge->nb, kind,0);
        nn+=genscdist(inge->plane[i+1],geom->site+nn,inge->trv,ofs,nx,ny,inge->nb, ord,pl,kind,0);
        pl++;
    }

    LOGMSG(5,"maketransgen plane idx after left base %d",pl);
    
    if ((repn>0)) {
        for (i=0;i<3;++i) {
            trans[i]=(inge->plane[repe+1]->coord[i])-(inge->plane[reps+1]->coord[i]);
        }
        LOGMSG(2,"will perform replication with tv=[%le,%le,%le]",trans[0],trans[1],trans[2]);

        for (i=reps;i<repe;i++) {
            if (plord==0){
                genchemdist(inge->plane[i+1], N*nreps, inge->nb, kind, N*nreps*(i-reps));
            }else{
                for(l=0;l<nreps;l++){
                    genchemdist(inge->plane[i+1], N, inge->nb, kind, N*nreps*(i-reps)+N*l);
                }
            }
        }
        for (l=0;l<nreps;l++) {
            if (l>0) {ofs[0]+=trans[0];ofs[1]+=trans[1];ofs[2]+=trans[2];}
            for (i=reps;i<repe;i++) {
                    nn+=genscdist(inge->plane[i+1],geom->site+nn,inge->trv,ofs,nx,ny,inge->nb, ord,pl,kind,N*l+N*nreps*(i-reps));
                    
                    pl++;                    
            }
        }
    }
    
    LOGMSG(5,"maketransgen plane idx after replication %d",pl);
    
    for (i=repe;i<(inge->np);i++){
        genchemdist(inge->plane[i+1], N, inge->nb, kind,0);        
        nn+=genscdist(inge->plane[i+1],geom->site+nn,inge->trv,ofs,nx,ny,inge->nb, ord,pl,kind,0);
        pl++;
    }

    LOGMSG(5,"maketransgen plane idx aftet right base %d",pl);


    for (i=0;i<nrb;i++){
        genchemdist(inge->plane[inge->np+1], N, inge->nb, kind,0);        
        nn+=genscdist(inge->plane[inge->np+1],geom->site+nn,inge->trv,ofs,nx,ny,inge->nb, ord,pl,kind,0);
        ofs[0]+=inge->rtrans[0];ofs[1]+=inge->rtrans[1];ofs[2]+=inge->rtrans[2];
        pl++;
    }        

    LOGMSG(5,"maketransgen plane idx aftet right replication %d",pl);

    LOGMSG(5,"maketransgeom : point 2, pl=%d",pl);
    for (site=geom->site;site<(geom->site+inge->nb*N);site++) site->mask=0;
    for (site=geom->site+nn-1;site>=(geom->site+nn-inge->nb*N);site--) site->mask=0;

    for (site=geom->site;site<(geom->site+inge->nb*nskip*N);site++) site->dmask=0;
    for (site=geom->site+nn-1;site>=(geom->site+nn-inge->nb*nskip*N);site--) site->dmask=0;

    geom->nmtr=ceil(inge->nmtr/((double)(MIN(nx,ny))));
    for (i=0;i<2;i++) {geom->trv[0][i]=inge->trv[0][i]*nx;geom->trv[1][i]=inge->trv[1][i]*ny;}
    geom->cr=inge->cr;
    geom->na=na;
    geom->sc[0]=geom->sc[1]=1;
    memcpy(geom->scx,inge->scx,3*sizeof(double));
    geom->tp=1;
    geom->dimcf=inge->dimcf;
    memset(geom->trans,0,sizeof(double)*3);
    // lead geometries
    lg->site=(t_trsite  *)malloc(sizeof(t_trsite )*inge->nb);    
    rg->site=(t_trsite  *)malloc(sizeof(t_trsite )*inge->nb);    

    genchemdist(inge->plane[0], 1, inge->nb, kind,0); 
    genscdist(inge->plane[0],lg->site,inge->trv,ofs0,1,1,inge->nb, ord,0,kind,0);

    genchemdist(inge->plane[inge->np+1], 1, inge->nb, kind,0); 
    genscdist(inge->plane[inge->np+1],rg->site,inge->trv,ofs,1,1,inge->nb, ord,0,kind,0);

    LOGMSG(5,"maketransgeom : point 4");

    lg->nmtr=inge->nmtr;
    lg->cr=inge->cr;
    lg->na=inge->nb;
    lg->sc[0]=nx;
    lg->sc[1]=ny;
    lg->tp=1;
    lg->dimcf=inge->dimcf;
    memcpy(lg->scx,inge->scx,3*sizeof(double));
    memcpy(lg->trv,inge->trv,4*sizeof(double));
    memcpy(lg->trans,inge->ltrans,3*sizeof(double));

    rg->nmtr=inge->nmtr;
    rg->cr=inge->cr;
    rg->na=inge->nb;
    rg->sc[0]=nx;
    rg->sc[1]=ny;
    rg->tp=1;
    rg->dimcf=inge->dimcf;
    memcpy(rg->trv,inge->trv,4*sizeof(double));
    memcpy(rg->scx,inge->scx,3*sizeof(double));
    memcpy(rg->trans,inge->rtrans,3*sizeof(double));
    LOGMSG(5,"maketransgeom : point 6");

    extrafields(geom);
    extrafields(lg);
    extrafields(rg);
    LOGMSG(5,"maketransgeom : point 7");
}
//----------------------------

//void makeseptrgeom(t_scgeom *inge, t_trgeom *geom, t_trgeom *lg, t_trgeom *rg, int nx,int ny,int nlb,int nrb, int ord,int reps,int repn,
//                   int nreps,int nskip,int plord,int nlayer, char *sink){
//
//    int na,N,nn,i,j,k,pl,l;
//    char *bla;
//    double ofs[3]={0,0,0};
//    double ofs0[3]={0,0,0};
//    double trans[3]={0,0,0};
//    t_trsite *site;
//    int repe;
//    int **kind;
//    // int nreps;
//    bla=sink;
//    printf("atom_input:%s\n",bla);
//    // nreps=nr+1;
//    // repn=MAX(repn0,1);
//    repe=reps+repn;
//    N=nx*ny;
//    LOGMSG(5,"maketransgen call, N=%d, repe=%d, inge->nb=%d, nreps=%d,repn=%d, reps=%d",N,repe,inge->nb,nreps,repn,reps);
//
//    kind=malloc(sizeof(int *)*inge->nb);
//
//    LOGMSG(5,"");
//
//    for(i = 0; i < inge->nb; ++i)
//    {
//        kind[i]=malloc(sizeof(int)*N*nreps*MAX(repn,1));
//    }
//
//    // printf("%d\n",inge->nb*repn*(nreps-1));
//    na=(inge->nat+inge->nb*nlb+inge->nb*nrb+inge->nb*repn*(nreps-1))*N;
//    // printf("%d, %d, %d, %d, %d, %d\n",reps,repe,repn,nreps,na);
//
//    LOGMSG(5,"maketransgen call, NA=%d",na);
//    geom->site=(t_trsite  *)malloc(sizeof(t_trsite )*na);
//
//    nn=0;
//    pl=0;
//    for (i=0;i<nlb;i++){
//        ofs[0]+=inge->ltrans[0];ofs[1]+=inge->ltrans[1];ofs[2]+=inge->ltrans[2];
//        genchemdist(inge->plane[0], N, inge->nb, kind,0);
//        nn+=genscsepdist(inge->plane[0],geom->site+nn,inge->trv,ofs,nx,ny,inge->nb, ord,pl,kind,0,nlayer,sink);
//        pl++;
//    }
//
//    LOGMSG(5,"maketransgen plane idx after left replication %d",pl);
//
//    for (i=0;i<(reps);i++){
//        genchemdist(inge->plane[i+1], N, inge->nb, kind,0);
//        nn+=genscsepdist(inge->plane[i+1],geom->site+nn,inge->trv,ofs,nx,ny,inge->nb, ord,pl,kind,0,nlayer,sink);
//        pl++;
//    }
//
//    LOGMSG(5,"maketransgen plane idx after left base %d",pl);
//
//    if ((repn>0)) {
//        for (i=0;i<3;++i) {
//            trans[i]=(inge->plane[repe+1]->coord[i])-(inge->plane[reps+1]->coord[i]);
//        }
//        LOGMSG(2,"will perform replication with tv=[%le,%le,%le]",trans[0],trans[1],trans[2]);
//
//        for (i=reps;i<repe;i++) {
//            if (plord==0){
//                genchemdist(inge->plane[i+1], N*nreps, inge->nb, kind, N*nreps*(i-reps));
//            }else{
//                for(l=0;l<nreps;l++){
//                    genchemdist(inge->plane[i+1], N, inge->nb, kind, N*nreps*(i-reps)+N*l);
//                }
//            }
//        }
//        for (l=0;l<nreps;l++) {
//            if (l>0) {ofs[0]+=trans[0];ofs[1]+=trans[1];ofs[2]+=trans[2];}
//            for (i=reps;i<repe;i++) {
//                nn+=genscsepdist(inge->plane[i+1],geom->site+nn,inge->trv,ofs,nx,ny,inge->nb, ord,pl,kind,N*l+N*nreps*(i-reps),nlayer,sink);
//
//                pl++;
//            }
//        }
//    }
//
//    LOGMSG(5,"maketransgen plane idx after replication %d",pl);
//
//    for (i=repe;i<(inge->np);i++){
//        genchemdist(inge->plane[i+1], N, inge->nb, kind,0);
//        nn+=genscsepdist(inge->plane[i+1],geom->site+nn,inge->trv,ofs,nx,ny,inge->nb, ord,pl,kind,0,nlayer,sink);
//        pl++;
//    }
//
//    LOGMSG(5,"maketransgen plane idx aftet right base %d",pl);
//
//
//    for (i=0;i<nrb;i++){
//        genchemdist(inge->plane[inge->np+1], N, inge->nb, kind,0);
//        nn+=genscsepdist(inge->plane[inge->np+1],geom->site+nn,inge->trv,ofs,nx,ny,inge->nb, ord,pl,kind,0,nlayer,sink);
//        ofs[0]+=inge->rtrans[0];ofs[1]+=inge->rtrans[1];ofs[2]+=inge->rtrans[2];
//        pl++;
//    }
//
//    LOGMSG(5,"maketransgen plane idx aftet right replication %d",pl);
//
//    LOGMSG(5,"maketransgeom : point 2, pl=%d",pl);
//    for (site=geom->site;site<(geom->site+inge->nb*N);site++) site->mask=0;
//    for (site=geom->site+nn-1;site>=(geom->site+nn-inge->nb*N);site--) site->mask=0;
//
//    for (site=geom->site;site<(geom->site+inge->nb*nskip*N);site++) site->dmask=0;
//    for (site=geom->site+nn-1;site>=(geom->site+nn-inge->nb*nskip*N);site--) site->dmask=0;
//
//    geom->nmtr=ceil(inge->nmtr/((double)(MIN(nx,ny))));
//    for (i=0;i<2;i++) {geom->trv[0][i]=inge->trv[0][i]*nx;geom->trv[1][i]=inge->trv[1][i]*ny;}
//    geom->cr=inge->cr;
//    geom->na=na;
//    geom->sc[0]=geom->sc[1]=1;
//    memcpy(geom->scx,inge->scx,3*sizeof(double));
//    geom->tp=1;
//    geom->dimcf=inge->dimcf;
//    memset(geom->trans,0,sizeof(double)*3);
//    // lead geometries
//    lg->site=(t_trsite  *)malloc(sizeof(t_trsite )*inge->nb);
//    rg->site=(t_trsite  *)malloc(sizeof(t_trsite )*inge->nb);
//
//    genchemdist(inge->plane[0], 1, inge->nb, kind,0);
//    genscsepdist(inge->plane[0],lg->site,inge->trv,ofs0,1,1,inge->nb, ord,0,kind,0,nlayer,sink);
//
//    genchemdist(inge->plane[inge->np+1], 1, inge->nb, kind,0);
//    genscsepdist(inge->plane[inge->np+1],rg->site,inge->trv,ofs,1,1,inge->nb, ord,0,kind,0,nlayer,sink);
//
//    LOGMSG(5,"maketransgeom : point 4");
//
//    lg->nmtr=inge->nmtr;
//    lg->cr=inge->cr;
//    lg->na=inge->nb;
//    lg->sc[0]=nx;
//    lg->sc[1]=ny;
//    lg->tp=1;
//    lg->dimcf=inge->dimcf;
//    memcpy(lg->scx,inge->scx,3*sizeof(double));
//    memcpy(lg->trv,inge->trv,4*sizeof(double));
//    memcpy(lg->trans,inge->ltrans,3*sizeof(double));
//
//    rg->nmtr=inge->nmtr;
//    rg->cr=inge->cr;
//    rg->na=inge->nb;
//    rg->sc[0]=nx;
//    rg->sc[1]=ny;
//    rg->tp=1;
//    rg->dimcf=inge->dimcf;
//    memcpy(rg->trv,inge->trv,4*sizeof(double));
//    memcpy(rg->scx,inge->scx,3*sizeof(double));
//    memcpy(rg->trans,inge->rtrans,3*sizeof(double));
//    LOGMSG(5,"maketransgeom : point 6");
//
//    extrafields(geom);
//    extrafields(lg);
//    extrafields(rg);
//    LOGMSG(5,"maketransgeom : point 7");
//}
//----------------------------

void extrafields(t_trgeom *geom){
    int i,j;
    double d=geom->dimcf;
    geom->rcoord=(double (*)[3])malloc(sizeof(double[3])*geom->na);
    for(i=0;i<geom->na;++i){
        for (j=0;j<3;++j){
            geom->rcoord[i][j]=geom->site[i].coord[j]*geom->scx[j]*d;            
        }
    }    
}


    static void usage(void){
        printf("prandgen.x - input files generator for transport code\n");
        printf("Usage:\n"
            "\tprandgen.x [options]\n"
            "\n");
        printf("Options:\n"
            "\t-s <size>    supercell size, number or NxM\n"
            "\t-b <num>     number of bulk layers to add on both sides\n"
            "\t-l <num>     number of bulk layers to add on LEFT side\n"
            "\t-r <num>     number of bulk layers to add on RIGHT side\n"
            "\t               '-r' and '-l' have priority over '-b'\n"   
            "\t-e <Nstart,Num,Nrep>  define layers to replicate in the scattering region\n"
            "\t-a <a0>      dimensional unit (lattice constant in most cases)\n"
            "\t               in atomic units. By default calculated from atomic Wsr's\n"
            "\t-c <scl>     scale resulting coordinates using <scl> as factor(not used)\n"
            "\t-o           use alternative ordering\n"
            "\t-m           write mask for damping calculations/temperature/complex energy calcs\n"
            "\t-f           do not generate visualisation files\n"        
            "\t-t <Tdeb,T>  used for temperature vibrations\n"
            "\t-d Theta     M2(Theta) for magnetic disorder\n"
            "\t-x <Tc,T>    Magnetic disorder with fitted M(T/Tc)\n"
            "\t-y <nl>      Skip number of layer on the sides for lattice/magnetic disorder\n"
            "\t               <Tdeb> - Debye temperature, <T> - temperature, in K\n"        
            "\t-z           use correlated Debye model for vibrations(not ready)\n"        
            "\t-v <num>     verbose output with <num> level of detalization\n"        
            "\t--vmask      use vibration mask from file \"ivmask\"\n"        
            "\t--rmask      use rotation mask from file \"irmask\"\n"
            "\t--rms        use rms input from file \"rms\"\n"
            "\t-i           configuration for injection (i.e. bilayer) in y direction \n"
            "\t-j           configuration for trilayer in y direction \n"
            "\t-h           Display help.\n"
            "\n");
        exit(0);

    }

    void readopts(int argc, char *argv[], t_opt *par){
        int opt, iproc=0;
        char *st1,*st2;
        char buf[400];
        int fvibmask=0,frotmask=0,fperlayer=0,fvib=0;
        struct option long_options[] =
          {
            /* These options set a flag. */
            {"vmask", no_argument,       &fvibmask, 1},
            {"rmask",   no_argument,     &frotmask, 1},
            {"perlayer",   no_argument,     &fperlayer, 1},
              {"rms",   no_argument,     &fvib, 1},
            /* These options don't set a flag.
               We distinguish them by their indices. */
/*            {"add",     no_argument,       0, 'a'},
            {"append",  no_argument,       0, 'b'},
            {"delete",  required_argument, 0, 'd'},
            {"create",  required_argument, 0, 'c'},
            {"file",    required_argument, 0, 'f'},*/
            {0, 0, 0, 0}
          };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        par->ny=par->nx=1;
        par->lbnl=par->rbnl=0;
        par->ord=0;
        par->vis=1;
        par->mask=0;
        par->scl=1.0;
        par->Tdeb=par->TK=-1.;
        par->Tc=par->Tm=-1.;
        par->Mdev=-1.0;
        par->dimcf=-1.0;
        par->mdef=1.0;
        par->Mtheta=0.0;
        par->Mphi=0.0;
        par->corr=0;
        par->nreps=1;
        par->reps=0;
        par->repn=0;
        par->vibmask=0;
        par->rotmask=0;
        par->skip_dis=1;
        par->plord=0;
        par->vib=0;
        par->pin=0;
        par->nlayer=1000;
        while ((opt=getopt_long(argc,argv,"+s:v:b:e:l:r:x:y:a:t:d:c:i:j:ozmfh",long_options, &option_index))!=-1){
            switch(opt){
                case 0:
                    // if (long_options[option_index].flag != 0) break;
                    // printf ("option %s", long_options[option_index].name);
                    // if (optarg)
                    // printf (" with arg %s", optarg);
                    // printf ("\n");
                break;
                                   
                case 's':                
                strcpy(buf,optarg);
                st1=buf;
                st2=strsep(&st1,"xX");
                if (st1!=NULL){
                    par->nx=atoi(st2);
                    par->ny=atoi(st1);
                    if (par->nx==0 || par->ny==0) usage();
                }else{
                    par->ny=par->nx=atoi(optarg);
                }
                break;

                case 'b':
                if (par->lbnl==0) par->lbnl=atoi(optarg);
                if (par->rbnl==0) par->rbnl=atoi(optarg);
                break;

                case 'l':
                par->lbnl=atoi(optarg);
                break;

                case 'r':
                par->rbnl=atoi(optarg);
                break;

                case 'a':
                par->dimcf=atof(optarg);
                break;

                case 'f':
                par->vis=0;
                break;

                case 'z':
                par->corr=1;
                break;

                case 'm':
                par->mask=1;
                break;

                case 'y':
                par->skip_dis=MAX(atoi(optarg),1);
                break;

                case 't':
                strcpy(buf,optarg);
                st1=buf;
                st2=strsep(&st1,",");
                if (st1!=NULL){
                    par->Tdeb=atof(st2);
                    st2=strsep(&st1,",");
                    if (st1!=NULL){
                        par->mdef=atof(st1);
                    }
                    par->TK=atof(st2);
                }else{
                    usage();
                }                    
                break;

                case 'x':
                strcpy(buf,optarg);
                st1=buf;
                st2=strsep(&st1,",");
                if (st1!=NULL){
                    par->Tc=atof(st2);
                    par->Tm=atof(st1);
                }else{
                    usage();
                }                    
                break;

                case 'e':
                strcpy(buf,optarg);
                st1=buf;
                st2=strsep(&st1,",");
                if (st1!=NULL){
                    par->reps=atoi(st2)-1;
                    st2=strsep(&st1,",");
                    if (st1!=NULL){
                        par->nreps=atoi(st1);
                    }
                    par->repn=atoi(st2);
                    if (par->nreps<1) {
                        printf("Error! Nrep must be >=1");
                        exit(1);
                    }
                    if (par->reps<0) {
                        printf("Error! Nstart must be >=1");
                        exit(1);
                    }
                    if (par->repn<1) {
                        printf("Error! Nnum must be >=1");
                        exit(1);
                    }
                    
                }else{
                    usage();
                }                    
                break;

                case 'd':
                par->Mdev=atof(optarg);
                // strcpy(buf,optarg);
                // st1=buf;
                // st2=strsep(&st1,",");
                // if (st1!=NULL){
                //     par->Mtheta=atof(st2);
                //     st2=strsep(&st1,",");
                //     if (st1!=NULL){
                //         par->Mdev=atof(st1);
                //     }
                //     par->Mphi=atof(st2);
                // }else{
                //     usage();
                // }                    
                break;

                case 'c':
                par->scl=atof(optarg);
                break;

                case 'o':
                par->ord=1;
                break;

                case 'v':
                LOG_MAX=1;
                if (optarg) LOG_MAX=atoi(optarg);
                LOGMSG(1,"loglevel is set to %d",LOG_MAX);
                break;

                case '?':
                usage();
                break;

                case 'h':
                usage();
                break;

                case 'i':
                par->pin=1;
                strcpy(buf,optarg);
                st1=buf;
                if (st1!=NULL){
                    par->nlayer=atoi(st2);
                }
                break;

                case 'j':
                if (par->pin!=1) {
                  printf("you cannot use -j without -i");
                  exit(0);
                }
                par->pin=2;
                strcpy(buf,optarg);
                st1=buf;
                if (st1!=NULL){
                    par->nlayer2=atoi(st2);
                }
                break;
                
                default:
                usage();
                exit(0);
            }
            iproc=1;
        }
        if (!iproc) usage();
        if (fperlayer){
            par->plord=1;
        }
        printf("SuperCell size    : %dx%d\n",par->nx,par->ny);
        printf("Bulk layers to add: Left=%d,Right=%d\n",par->lbnl,par->rbnl);
        printf("SC Ordering       : %s\n",(par->ord==0)?"Default":"Alternative");
        printf("Fix concentration : %s\n",(par->plord==0)?"per slab":"per layer");            
        
        printf("Clean side-layers : %d\n",par->skip_dis);
    // printf("Scaling factor = %f, ordering = %d (0 - default)\n",par->scl,par->ord);
        if (par->Tdeb>=0.0) {
            printf("Thermal vibrations: Tdeb=%f, T=%f, %s\n",par->Tdeb,par->TK,(par->corr==1)?"correlated":"non-correlated");
            par->lbnl=MAX(par->lbnl,1);
            par->rbnl=MAX(par->rbnl,1);
            if (fvibmask){
                printf("                  : use mask from file \"ivmask\" \n");
                par->vibmask=1;
            }
        }else if(par->vib==0.0){
            printf("Thermal vibrations: none\n");        
        }
        if (par->Tc>=0.0){
            printf("Spin disorder with prefitted model: Tc=%f, T=%f\n",par->Tc,par->Tm);
            par->Mdev=getMdev(par->Tm/par->Tc);
        }

        if (par->Mdev>=0.0){
            printf("Magnetic config: given M2(Q)=%f, Mz=%f\n",par->Mdev,exp(-par->Mdev*par->Mdev*M_PI*M_PI/2.0));
            if (frotmask){
                par->rotmask=1;
                printf("                  : use mask from file \"irmask\" \n");
            }
        }
        if (par->pin==1) {
            printf("Configuration for current injection in y direction at layer=%i\n",par->nlayer);
        }
        if (par->pin==2) {
            printf("Configuration for third layer in y direction at layer=%i\n",par->nlayer2);
        }
        if(fvib==1){
            par->vib=1;
            }
    }

    int main(int argc, char *argv[])
    {
        t_scgeom scgeo;
        t_atomset atoms;
        t_opt par;
        t_trgeom mgeom,lg,rg;
        t_visdat vis;
        LOGMSG(2,"Supercell generator\n");
        initrnd();
        readopts(argc, argv, &par);

        LOGMSG(1,"Finished options processing");
        if(par.pin==1.0){
            readinput2(&scgeo, &atoms);
        }
        else if(par.pin==2.0){
            readinput3(&scgeo, &atoms);
        }
        else{
            readinput(&scgeo, &atoms);
        }
        LOGMSG(1,"Input information is loaded");


        loadmasks(&scgeo,par.rotmask,par.vibmask);
        LOGMSG(1,"Masks are loaded");

        
        if (par.dimcf>0.0) scgeo.dimcf=par.dimcf;

    //     if (mdef!=0.) printf("umsqr=%f, %f\n",devmeansqr(TK,Tdeb,mdef)/dimcf*sqrt(3.),devmeansqr(TK,Tdeb,mdef)/rawsr*sqrt(3.));
    // if (mdef!=0.) printf("umsqr=%f, %f\n",devmeansqr(TK,Tdeb,mdef)/dimcf*sqrt(3.),devmeansqr(TK,Tdeb,mdef)/rawsr*sqrt(3.));
        maketrgeom(&scgeo, &mgeom, &lg, &rg, par.nx,par.ny,par.lbnl,par.rbnl, par.ord,par.reps,par.repn,par.nreps,par.skip_dis,par.plord);
        LOGMSG(5,"WA1");

        writeatomlist(&atoms);
        LOGMSG(5,"WA2");
        if(par.pin==1.0){
 //           printf("sink:%s\n",par.sink);
//            genpingeom(&mgeom,par.nlayer,par.sink,1);
//            genpingeom(&lg,par.nlayer,par.sink,0);
//            genpingeom(&rg,par.nlayer,par.sink,0);
            genpingeom(&mgeom,&scgeo,par.nlayer,1);
            genpingeom(&lg,&scgeo,par.nlayer,0);
            genpingeom(&rg,&scgeo,par.nlayer,0);
        }
        else if(par.pin==2.0){
            genpingeom_trilayer(&mgeom,&scgeo,par.nlayer,par.nlayer2,1);
            genpingeom_trilayer(&lg,&scgeo,par.nlayer,par.nlayer2,0);
            genpingeom_trilayer(&rg,&scgeo,par.nlayer,par.nlayer2,0);
        }
        else{
            genpingeom(&mgeom,&scgeo,par.nlayer,0);
            genpingeom(&lg,&scgeo,par.nlayer,0);
            genpingeom(&rg,&scgeo,par.nlayer,0);
//            genpingeom(&mgeom,par.nlayer,par.sink,0);
//            genpingeom(&lg,par.nlayer,par.sink,0);
//            genpingeom(&rg,par.nlayer,par.sink,0);
        };
        if (par.Tdeb>=0.0) {
            if (par.corr==1){
                getvibrations(&mgeom,par.Tdeb, par.TK,par.mdef);
            }else {
                geteasyvibrations(&mgeom,par.Tdeb, par.TK,par.mdef);            
            }
        };
        if(par.vib==1) {
            getsepvibrations(&mgeom);
        };
        // geteasymagnons(&mgeom,par.Mtheta,par.Mphi,par.Mdev);        
        if (par.Mdev>=0.0){
            geteasymagnons(&mgeom,par.Mdev);  
            printf("Magnetic config: MAX(|Q|)=%f, M2(Q)=%f, M1(F)=%f\n",mgeom.mvals[2],mgeom.mvals[0],mgeom.mvals[1]);
            printf("Magnetic config: macrospin=(%f,%f), Mz=%f\n",mgeom.mvals[3],mgeom.mvals[4],mgeom.mvals[5]);
            writerotmask(&mgeom);
        }
        printf("Dimensional unit  : %f, %f\n",scgeo.dimcf, scgeo.dimcf*5.2917720859e-2);

        LOGMSG(5,"Main7");
        
        writegeom(&mgeom, "geom_m");

        LOGMSG(5,"Main8");
        
        writegeom(&lg, "geom_l");
        writegeom(&rg, "geom_r");

        if (par.mask!=0) writerotcf(&mgeom);
        if (par.vis!=0){
    	    LOGMSG(5,"Main9");
            prepvisdat(&mgeom,&vis);
    	    LOGMSG(5,"Main10");
            dumpvisdat(&vis);
        }
//     dumpvisdat(dimcf*5.291772108e-1,&vis);
        printf("Finished!\n");
    }

