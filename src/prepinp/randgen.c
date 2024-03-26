#include "libgen.h"

typedef struct{
    int ord, nx,ny,N,lbnl,rbnl;
    int vis,mask, corr;
    double scl;
    double dimcf,Tdeb,TK,mdef;    
    double Mtheta,Mphi,Mdev;
    } t_opt;


int genscdist(t_scsite *plane,t_trsite  *isite,double (*tr)[2],double *ofs,int nx,int ny,int nb, int ord, int pidx);
void maketrgeom(t_scgeom *inge, t_trgeom *geom, t_trgeom *lg, t_trgeom *rg, int nx,int ny,int nlb,int nrb, int ord);
void extrafields(t_trgeom *geom);
static void usage(void);
void readopts(int argc, char *argv[], t_opt *par);

inline int genscdist(t_scsite *plane,t_trsite  *isite,double (*tr)[2],double *ofs,int nx,int ny,int nb, int ord, int pidx){
    int N=nx*ny;
    int kind[nb][N];
    int j,at,nn,l,p;
    double pr;
    t_trsite  *site;
// generate atoms distribution    
    for (j=0;j<nb;++j){
        memset(kind[j],0,sizeof(int)*N);   
        for (at=1;at<((plane+j)->nc);at++) {
            nn=round(((plane+j)->conc[at])*N);
            l=0;
            while (l<nn) {
                pr=drand48()-1e-16;
                p=floor(N*pr);
                if (kind[j][p]==0){
                    l++;
                    kind[j][p]=at;
                };
            };
        };
    };

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
                site->at=(plane+at)->at[kind[at][nx*j+l]];
                site->pidx=pidx;
                site->mask=1;
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
                site->at=(plane+at)->at[kind[at][nx*j+l]];
                site->pidx=pidx;
                site->mask=1;
                site++;
            };
        };
    };


    return(N*nb);    
}


void maketrgeom(t_scgeom *inge, t_trgeom *geom, t_trgeom *lg, t_trgeom *rg, int nx,int ny,int nlb,int nrb, int ord){

    int na,N,nn,i,j,k,pl;
    double ofs[3]={0,0,0};
    double ofs0[3]={0,0,0};
    t_trsite *site;
    N=nx*ny;
    na=(inge->nat+inge->nb*nlb+inge->nb*nrb)*N;
    geom->site=(t_trsite  *)malloc(sizeof(t_trsite )*na);    
    nn=0;
    pl=0;
    for (i=0;i<nlb;i++){
        ofs[0]+=inge->ltrans[0];ofs[1]+=inge->ltrans[1];ofs[2]+=inge->ltrans[2];
        nn+=genscdist(inge->plane[0],geom->site+nn,inge->trv,ofs,nx,ny,inge->nb, ord,pl);
        pl++;
    }        

    for (i=0;i<(inge->np);i++){
        nn+=genscdist(inge->plane[i+1],geom->site+nn,inge->trv,ofs,nx,ny,inge->nb, ord,pl);
        pl++;
    }
    for (i=0;i<nrb;i++){
        nn+=genscdist(inge->plane[inge->np+1],geom->site+nn,inge->trv,ofs,nx,ny,inge->nb, ord,pl);
        ofs[0]+=inge->rtrans[0];ofs[1]+=inge->rtrans[1];ofs[2]+=inge->rtrans[2];
        pl++;
    }        
    for (site=geom->site;site<(geom->site+inge->nb*N);site++) site->mask=0;
    for (site=geom->site+nn-1;site>=(geom->site+nn-inge->nb*N);site--) site->mask=0;

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

    genscdist(inge->plane[0],lg->site,inge->trv,ofs0,1,1,inge->nb, ord,0);
    genscdist(inge->plane[inge->np+1],rg->site,inge->trv,ofs,1,1,inge->nb, ord,0);

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

    extrafields(geom);
    extrafields(lg);
    extrafields(rg);
}

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
        printf("randgen.x - input files generator for transport code\n");
        printf("Usage:\n"
            "\trandgen.x [options]\n"
            "\n");
        printf("Options:\n"
            "\t-s <size>    supercell size, number or NxM\n"
            "\t-b <num>     number of bulk layers to add on both sides\n"
            "\t-l <num>     number of bulk layers to add on LEFT side\n"
            "\t-r <num>     number of bulk layers to add on RIGHT side\n"
            "\t               '-r' and '-l' have priority over '-b'\n"            
            "\t-a <a0>      dimensional unit (lattice constant in most cases)\n"
            "\t               in atomic units. By default calculated from atomic Wsr's\n"
            "\t-c <scl>     scale resulting coordinates using <scl> as factor(not used)\n"
            "\t-o           use alternative ordering\n"
            "\t-m           write mask for damping calculations/temperature/complex energy calcs\n"
            "\t-f           do not generate visualisation files\n"        
            "\t-t <Tdeb,T>  used for temperature vibrations\n"
            "\t-d <Theta>   M2(Theta) for magnetic disorder\n"
            "\t               <Tdeb> - Debye temperature, <T> - temperature, in K\n"        
            "\t-z           use correlated Debye model for vibrations(not ready)\n"        
            "\t-h           Display help.\n"
            "\n");
        exit(0);

    }

    void readopts(int argc, char *argv[], t_opt *par){
        int opt, iproc=0;
        char *st1,*st2;
        char buf[400];

        par->ny=par->nx=1;
        par->lbnl=par->rbnl=0;
        par->ord=0;
        par->vis=1;
        par->mask=0;
        par->scl=1.0;
        par->Tdeb=par->TK=-1.;
        par->dimcf=-1.0;
        par->mdef=1.0;
        par->Mdev=0.0;
        par->Mtheta=0.0;
        par->Mphi=0.0;
        par->corr=0;
        while ((opt=getopt(argc,argv,"+s:b:l:r:a:t:d:c:ozmfh"))!=-1){ 
            switch(opt){
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

                case 'd':
                strcpy(buf,optarg);
                st1=buf;
                st2=strsep(&st1,",");
                if (st1!=NULL){
                    par->Mtheta=atof(st2);
                    st2=strsep(&st1,",");
                    if (st1!=NULL){
                        par->Mdev=atof(st1);
                    }
                    par->Mphi=atof(st2);
                }else{
                    usage();
                }                    
                break;

                case 'c':
                par->scl=atof(optarg);
                break;

                case 'o':
                par->ord=1;
                break;

                case '?':
                usage();
                break;

                case 'h':
                usage();
                break;

                default:
                usage();
                exit(0);
            }
            iproc=1;
        }
        if (!iproc) usage();
        printf("SuperCell size    : %dx%d\n",par->nx,par->ny);
        printf("Bulk layers to add: Left=%d,Right=%d\n",par->lbnl,par->rbnl);
        printf("Ordering          : %s\n",(par->ord==0)?"Default":"Alternative");
    // printf("Scaling factor = %f, ordering = %d (0 - default)\n",par->scl,par->ord);
        if (par->Tdeb>=0.0) {
            printf("Thermal vibrations: Tdeb=%f, T=%f, %s\n",par->Tdeb,par->TK,(par->corr==1)?"correlated":"non-correlated");
            par->lbnl=MAX(par->lbnl,1);
            par->rbnl=MAX(par->rbnl,1);
        }else{
            printf("Thermal vibrations: none\n");        
        }
        printf("Magnetic config: orientation=(%f ,%f), deviation=%f\n",par->Mtheta,par->Mphi,par->Mdev);
        
    }

    int main(int argc, char *argv[])
    {
        t_scgeom scgeo;
        t_atomset atoms;
        t_opt par;
        t_trgeom mgeom,lg,rg;
        t_visdat vis;
        initrnd();

        readopts(argc, argv, &par);

        readinput(&scgeo, &atoms);

        if (par.dimcf>0.0) scgeo.dimcf=par.dimcf;
        printf("Dimensional unit  : %f\n",scgeo.dimcf);

    //     if (mdef!=0.) printf("umsqr=%f, %f\n",devmeansqr(TK,Tdeb,mdef)/dimcf*sqrt(3.),devmeansqr(TK,Tdeb,mdef)/rawsr*sqrt(3.));
    // if (mdef!=0.) printf("umsqr=%f, %f\n",devmeansqr(TK,Tdeb,mdef)/dimcf*sqrt(3.),devmeansqr(TK,Tdeb,mdef)/rawsr*sqrt(3.));

        maketrgeom(&scgeo, &mgeom, &lg, &rg, par.nx,par.ny,par.lbnl,par.rbnl, par.ord);

        writeatomlist(&atoms);

        if (par.Tdeb>=0.0) {
            if (par.corr==1){
                getvibrations(&mgeom,par.Tdeb, par.TK,par.mdef);
            }else{
                geteasyvibrations(&mgeom,par.Tdeb, par.TK,par.mdef);            
            }
        };

        geteasymagnons(&mgeom,par.Mdev);        
        writegeom(&mgeom, "geom_m");
        writegeom(&lg, "geom_l");
        writegeom(&rg, "geom_r");

        if (par.mask!=0) writerotcf(&mgeom);
        writerotmask(&mgeom);
        if (par.vis!=0){
            prepvisdat(&mgeom,&vis);
            dumpvisdat(&vis);
        }
//     dumpvisdat(dimcf*5.291772108e-1,&vis);
        printf("Finished!\n");
    }

