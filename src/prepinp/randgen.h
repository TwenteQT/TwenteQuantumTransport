#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <getopt.h>
#include "masses.h"

#define kboltz 1.3806504e-23 //	J* K^−1
#define hbar 1.054571628e-34  // J*s
#define caem 1.660538782e-27 //kg
#define debcf  (hbar*hbar/(caem*kboltz))
#define bohrrad 0.52917720859e-1 // nm
#define aunit (5.291772108e-1) //a.u
#define ES_ii2 200
#define MAX_NUM_BOND 30
#define MAX_NC  150

#if !defined(round)
double round(double x);
#endif
#if !defined(MAX)
#define MAX(A,B) (A>B?A:B)
#endif

#if !defined(MIN)
#define MIN(A,B) (A<B?A:B)
#endif


#define atfmt " %-15s   %s   \t\t# %f\n"
#define fmt3s " %#15.11f %#15.11f %#15.11f   %s\n"
#define fmt2s " %#15.11g%#15.11g                  %s\n"

#define fpvec3(f,a,cm) \
    fprintf(f,fmt3s,a[0],a[1],a[2],cm)

#define fpvec2(f,a,cm) \
    fprintf(f,fmt2s,a[0],a[1],cm)
    
typedef struct{int n;double r;} cstr;
        
typedef struct{
    double *x,*y,*z,*wsr;
    int *ii1,*ii2,*plane;
    char (*info)[5];
    int num;
    int nuniq,nuniqa;
    char **ulabels,**ulats;
    int *bonds1,*bonds2,nbond;
} t_visdat;


typedef struct {
    double wsr,charge,mass;
    char *label;
    void *next;
    } t_atom;

typedef struct{
    int num;
    t_atom *ptr;
    } t_atomset;

typedef struct {
    double coord[3];
    int nc;
    double *conc;
    t_atom **at;
    } t_scsite;

typedef struct {
    int np,nb,nmtr;
    int nat, tnat;
    double cr;
    double dimcf;
    double scx[3],trv[2][2];
    double ltrans[3],rtrans[3];
    t_scsite *site;
    t_scsite **plane;
    } t_scgeom;

typedef struct {
    double coord[3];
    double mdir[2];
    int mask;
    int pidx;
    t_atom *at;
    } t_trsite;

typedef struct{
    int na,nmtr,tp;
    double cr;
    double dimcf;
    int sc[2];
    double scx[3],trv[2][2];
    t_trsite *site;
    double trans[3];
    double (*rcoord)[3];    
    } t_trgeom;

typedef struct{
    int ord, nx,ny,N,lbnl,rbnl;
    int vis,mask, corr;
    double scl;
    double dimcf,Tdeb,TK,mdef;    
    double Mtheta,Mphi,Mdev;
    } t_opt;


void scanv3(char *buf,FILE *f,double *a );
void scanv2(char *buf,FILE *f,double *a );
FILE *openread(char *nm);
double getmass(double charge);
t_atom *getatom(t_atomset *atoms,char *label);
int readinput(t_scgeom *geom, t_atomset *atoms);
double devmeansqr(double T0, double Tdeb);
void calcumsqr(t_atomset *atoms, double Tdeb, double TK, double dimcf);
void writeatomlist(t_atomset *atoms);
void writegeom(t_trgeom *geom, char *name);
void writerotcf(t_trgeom *geom);
double randn();
void readopts(int argc, char *argv[], t_opt *par);
void initrnd();
void maketrgeom(t_scgeom *inge, t_trgeom *geom, t_trgeom *lg, t_trgeom *rg, int nx,int ny,int nlb,int nrb, int ord);
int genscdist(t_scsite *plane,t_trsite  *isite,double (*tr)[2],double *ofs,int nx,int ny,int nb, int ord, int pidx);
void extrafields(t_trgeom *geom);
int cholesky(int N,int *Ap,int *Ai,double *Ax,int **Lp,int **Li,double **Lx, double **D);
int get_cf_dists(t_trgeom *geom, int** ApI,int **AiI, double **res, int **idx);
void geteasyvibrations(t_trgeom *geom,double Tdeb, double Tk,double mdef);
void getvibrations(t_trgeom *geom,double Tdeb, double Tk,double mdef);
void rotspheric (double *basei, double *anglei,double *nang);
void geteasymagnons(t_trgeom *geom,double Mtheta, double Mphi,double Mdev);
void writerotmask(t_trgeom *geom);