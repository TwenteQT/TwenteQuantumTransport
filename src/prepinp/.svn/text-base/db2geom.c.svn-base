#include "libgen.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <limits.h>

#define MAXSTRLEN 31

typedef struct{
    int vis, db, dry;
    char *dbloc;
} t_xopt;

// variables and parameters:

typedef struct t_varinput t_varinput;
struct t_varinput {
    char *name;
    char *value;
    t_varinput *next;
    int nread;
};

typedef union {
    int i;
    double d;
} t_value;

typedef struct t_varprops t_varprops;
struct t_varprops{
    char name[MAXSTRLEN];
    int type, n;
    t_value min;
    t_value max;
};

typedef struct t_par t_par;
typedef struct t_valptr t_valptr;
struct t_valptr{
    t_value *val;
    t_valptr *next;
    t_varprops *vp;
};
struct t_par {
    int nptr, nval;
    int type;
    t_valptr *ptr;
    t_value *val;
    char symbol;
    t_value min, max;
    t_par *next;
};

typedef struct t_var t_var;
struct t_var {
    int nset;
    t_varprops *vp;
    t_value *val;
};

// rest:

typedef struct {
    int nb, nmtr;
    double cr;
    double tr[3][3];
    t_scsite *site;
    double maxwsr;
} t_bulk;

typedef struct {
    int nb, nmtr, tp, ncell;
    double cr;
    double tr1[2], tr2[2], trperp[3];
    t_scsite *site;
    double maxwsr;
    double vol, volbelow;
} t_layer;

typedef struct {
    int num, nfix, nleft;
    double *initial;
    double *natexp;
    int npercell;
} t_basket;

typedef struct t_db_entry t_db_entry;
struct t_db_entry {
    int nvar;
    int kind; // see comment below
    char *name;
    t_varprops *vp;
    t_db_entry *next;
    t_bulk bulk;
    t_atomset atoms;
    t_basket basket;
};
// kind can be:
// 1 : turek bulk
// 2 : turek interface -> not implemented
// 3 : stuttgart supercell -> not implemented
// ... add and implement new kind of database entry

t_varprops dbe1_props[] = {
    {.name = "scvec1", .type = 1, .n = 3, .min = {.i = -INT_MAX}, .max = {.i = INT_MAX}},
    {.name = "scvec2", .type = 1, .n = 3, .min = {.i = -INT_MAX}, .max = {.i = INT_MAX}},
    {.name = "scvec3", .type = 1, .n = 3, .min = {.i = -INT_MAX}, .max = {.i = INT_MAX}},
    {.name = "shift", .type = 2, .n = 2, .min = {.d = -INFINITY}, .max = {.d = INFINITY}},
    {.name = "nlay", .type = 1, .n = 1, .min = {.i = 1}, .max = {.i = INT_MAX}},
    {.name = "nfixc", .type = 1, .n = 1, .min = {.i = 1}, .max = {.i = INT_MAX}},
    {.name = "stiff", .type = 2, .n = 1, .min = {.d = 0.0}, .max = {.d = INFINITY}},
    {.name = "rms", .type = 2, .n = 1, .min = {.d = 0.0}, .max = {.d = INFINITY}},
    {.name = "temp", .type = 2, .n = 1, .min = {.d = 0.0}, .max = {.d = INFINITY}},
    {.name = "tdeb", .type = 2, .n = 1, .min = {.d = 0.0}, .max = {.d = INFINITY}},
    {.name = "magdir", .type = 2, .n = 3, .min = {.d = -INFINITY}, .max = {.d = INFINITY}}
};

typedef struct t_slab t_slab;

struct t_slab{
    char *material;
    int nvarsin, nvar;
    t_varinput *varsin;
    t_var *var;
    t_db_entry *dbe;
    t_layer layer;
    double stiff;
};

typedef struct{
    int ns, npar, sc[2], nconf;
    t_slab *slab;
    t_varinput *parsin;
    t_par *par, *parptr;
    t_atomset atoms;
} t_multilay;

typedef struct{
    double (*coord)[3];
    t_trgeom lg, rg, mg;
    int ns, *na, *nfix;
    t_basket **basket;
    t_scsite **laysite;
    int *diskind, nint;
    double *rms, *zint;
    //int *sidx;
} t_config;

#define COMPLAIN_EXIT {LOGMSG(0,"Wrong input!");exit(1);};
#define EXITMSG(...) {printf("STOPPED\n");LOGMSG(0,__VA_ARGS__);exit(1);}; 

char const dbdirname[] = "potdb";
double const accuracy = 1.0e-15;

static void usage(void);
void readopts(int argc, char *argv[], t_xopt *xopt);
int read_input_ml(t_multilay *ml);
int count_char(char *buf, char c);
void store_var(char *buf, t_varinput *varsin);
FILE *cat_openread(char *name, char *path);
void get_db_dir(char *argv[], t_xopt *xopt);
int test_dir(char *dbloc1);
char * ask_dbloc();
int read_db_entry(t_db_entry *dbe, char *name, char *dbloc);
void strallcpy(char **str2, char *str1);
int read_db(t_db_entry **dbe, t_multilay *ml, char *dbloc);
int read_val(t_var *var, char *valstr, t_multilay *ml);


void alloc_basket(t_basket *a, int num) {
    a->num = num;
    a->initial = (double *)calloc(num, sizeof(double));
    a->natexp = (double *)calloc(num, sizeof(double));
    a->nleft = 0;
}    
    
int getatind(t_atom *at, t_atomset *atoms) {
    int i;
    t_atom *ptr;
    for(i=0, ptr=atoms->ptr; i<atoms->num; i++, ptr=ptr->next) {
        if (ptr == at) break;
    }
    return(i);
}

void reset_basket(t_basket *a) {
    int i;
    for(i=0; i<a->num; i++) {
        a->natexp[i] = a->initial[i]*a->nfix;
    }
    a->nleft = a->nfix*a->npercell;
}

int check_round(double num) {
    return fabs(num - round(num)) <= accuracy;
}

t_var *getvar(t_slab *slab, char *name, int opt) {
    int i;
    for(i=0; i<slab->nvar; i++) {
        if(strcmp(slab->var[i].vp->name, name) == 0) {
            if (slab->var[i].nset > 0) {
                return(&(slab->var[i]));
            }
            break;
        }
    }
    if (!opt) {
        EXITMSG("I need variable %s for %s", name, slab->material);
    }
    return(NULL);
}


void str2val(t_value *val, char *str, int type){
    int flag;
    switch(type) {
        case 1:
        val->i = scaninf(str)*INT_MAX + atoi(str);
        LOGMSG(5,"%d",val->i);
        break;
        
        case 2:
        val->d = atof(str);
        flag = scaninf(str);
        if (flag) val->d = (double)flag*INFINITY;
        LOGMSG(5,"%le",val->d);
        break;
        
        default:
        EXITMSG("str2val doesn't know that type");
        break;
    }
}

int test_minmax(t_value val, t_value min, t_value max, int type) {
    switch (type) {
        case 1:
        return val.i < min.i || val.i > max.i;
        break;
        case 2:
        return val.d < min.d || val.d > max.d;
        break;
        default:
        EXITMSG("test_minmax doesn't know that type");
        break;
    }
}

int scaninf(char *str) {
    int flag=0;
    double dummy;
    if (strstr(str, "inf") !=NULL) flag = 1;
    if (strstr(str, "Inf") !=NULL) flag = 1;
    if (strstr(str, "INF") !=NULL) flag = 1;
    if (strchr(str, '-') !=NULL) flag = -flag;
    if (flag) {
        return(flag);
    } else {
        if(sscanf(str,"%le",&dummy)!=1) COMPLAIN_EXIT;
        return(flag);
    }
}

char *posafter(char *str) {
    char *pos;
    pos = strchr(str, '{');
    if (pos != NULL) {
        pos++;
    } else {
        pos = str;
    }
    return(pos);
}
    
int count_char(char *buf, char c){
    int i=0, n=0;
    char *pos;
    while(NULL != (pos = strchr(buf+i, c))){
        i=pos-buf+1;
        n++;
    }
    return(n);
}

int count_str(char *buf, char *str){
    int i=0, n=0;
    char *pos;
    while(NULL != (pos = strstr(buf+i, str))){
        i=pos-buf+strlen(str);
        n++;
    }
    return(n);
}

void strallcpy(char **str2, char *str1) {
    *str2 = (char *)malloc(sizeof(char)*(strlen(str1)+1));
    strcpy(*str2, str1);
}

double getvolume(double *a, double *b, double *c) {
    double vol=0.0;
    vol += a[0]*(b[1]*c[2]-c[1]*b[2]);
    vol -= a[1]*(b[0]*c[2]-c[0]*b[2]);
    vol += a[2]*(b[0]*c[1]-c[0]*b[1]);
    return(fabs(vol));
}

double triple(double *a, double *b, double *c) {
    // a . b x c
    double vol=0.0;
    vol += a[0]*(b[1]*c[2]-c[1]*b[2]);
    vol -= a[1]*(b[0]*c[2]-c[0]*b[2]);
    vol += a[2]*(b[0]*c[1]-c[0]*b[1]);
    return(vol);
}

void inverse3(double res[3][3], double m[3][3]) {
    double vol = getvolume(&m[0][0],&m[1][0],&m[2][0]); //det of transpose = det
    res[0][0] = (m[1][1]*m[2][2]-m[1][2]*m[2][1])/vol;
    res[0][1] = (m[0][2]*m[2][1]-m[0][1]*m[2][2])/vol;
    res[0][2] = (m[0][1]*m[1][2]-m[0][2]*m[1][1])/vol;
    res[1][0] = (m[1][2]*m[2][0]-m[1][0]*m[2][2])/vol;
    res[1][1] = (m[0][0]*m[2][2]-m[0][2]*m[2][0])/vol;
    res[1][2] = (m[0][2]*m[1][0]-m[0][0]*m[1][2])/vol;
    res[2][0] = (m[1][0]*m[2][1]-m[1][1]*m[2][0])/vol;
    res[2][1] = (m[0][1]*m[2][0]-m[0][0]*m[2][1])/vol;
    res[2][2] = (m[0][0]*m[1][1]-m[0][1]*m[1][0])/vol;
}

void mmult3(double res[3][3], double m1[3][3], double m2[3][3]) {
    int i, j, k;
    for(i=0; i<3; i++) {
        for(j=0; j<3; j++) {
            res[i][j] = 0.0;
            for(k=0; k<3; k++) {
                res[i][j] += m1[i][k]*m2[k][j];
            }
        }
    }
}

void vmult3(double res[3], double m[3][3], double v[3]) {
    int i, k;
    for(i=0; i<3; i++) {
        res[i] = 0.0;
        for(k=0; k<3; k++) {
            res[i] += m[i][k]*v[k];
        }
    }
}

void cross3(double *res, double *a, double *b) {
    res[0] = a[1]*b[2]-a[2]*b[1];
    res[1] = a[2]*b[0]-a[0]*b[2];
    res[2] = a[0]*b[1]-a[1]*b[0];
}

double cross2(double *a, double *b) {
    return(a[0]*b[1]-a[1]*b[0]);
}

double inp3(double *vec1, double *vec2) {
    return(vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]);
}

double inp2(double *vec1, double *vec2) {
    return(vec1[0]*vec2[0] + vec1[1]*vec2[1]);
}

double abs3(double *vec) {
    return(sqrt(inp3(vec, vec)));
}

double abs2(double *vec) {
    return(sqrt(inp2(vec, vec)));
}

void zero3(double *vec) {
    vec[0] = 0.0;
    vec[1] = 0.0;
    vec[2] = 0.0;
}

void add3(double *res, double *v1, double *v2, double fact) {
    int i;
    for(i=0; i<3; i++) {
        res[i] = v1[i] + v2[i]*fact;
    }
}

double angle2(double *v1, double *v2) {
    double res;
    int sign;
    sign = v1[0]*v2[1]-v1[1]*v2[0] >= 0.0 ? 1 : -1;
    return( sign*acos(inp2(v1,v2)/abs2(v1)/abs2(v2)) );
}

int isinside(double *coord, double *v1, double *v2, double *v3) {
    double t, v;
    int i, s;
    v = triple(v1,v2,v3);
    s = v>0.0;
    v = fabs(v);
    
    t = triple(coord,v1,v2)*s;
    if (fabs(t/v)<=accuracy) t=0.0;
    if ( (t<0.0) || t/v >= 1.0-accuracy) return(0);
    
    t = triple(coord,v2,v3)*s;
    if (fabs(t/v)<=accuracy) t=0.0;
    if ( (t<0.0) || t/v >= 1.0-accuracy) return(0);
    
    t = triple(coord,v3,v1)*s;
    if (fabs(t/v)<=accuracy) t=0.0;
    if ( (t<0.0) || t/v >= 1.0-accuracy) return(0);
    
    return(1);
}

t_atom *getatom_db(t_atomset *atoms, char *label, char *path){
    t_atom *ptr;
    FILE *fl;
    int qq;
    for (ptr=atoms->ptr;ptr!=NULL; ptr=(ptr->next)){
        if (!strcmp(ptr->label,label)) break;
    }
    LOGMSG(2,"Loading atomic file for %s...",label);
    
    // Add new atom
    if (ptr==NULL) {
        char *full, buf[400];
        t_atom *nptr;
        
        full = (char *)malloc(sizeof(char)*(strlen(label)+strlen(path)+8));
        strcpy(full, path);
        strcat(full, "/atoms/");
        strcat(full, label);
        fl = cat_openread(full, "");
        nptr=(t_atom *)calloc(1, sizeof(t_atom));
        fgets(buf,400,fl);
        fgets(buf,400,fl);
        fgets(buf,400,fl);
        fgets(buf,400,fl);
        if(sscanf(buf,"%le %le",&(nptr->charge),&(nptr->wsr))!=2) COMPLAIN_EXIT;
        nptr->mass=getmass(nptr->charge);
        fclose(fl);
        nptr->label=(char *)malloc(sizeof(char)*(strlen(label)+1));
        strcpy(nptr->label,label);
        nptr->path = full;
        atoms->num++;
        nptr->next=atoms->ptr;
        atoms->ptr=nptr;
        ptr=nptr;
    }
    LOGMSG(2,"Done loading atomic file for %s!",label);    
    return(ptr);    
}

void read_bulk_turek(t_bulk *bulk, char *path, t_atomset *atoms) {
    char buf[400];
    FILE *FP;
    int i, j, k;
    double scaling[3];
    t_scsite *site;
    char label[32];
    double rawsr, dawsr, vol, dimcf, maxwsr;
    t_atom *ptr;
    
    LOGMSG(2,"Reading inpge");
    FP = cat_openread("inpge", path);
    fgets(buf,400,FP);
    fgets(buf,400,FP);
    if(sscanf(buf,"%d %d", &bulk->nb, &bulk->nmtr)!=2) COMPLAIN_EXIT;
    fgets(buf,400,FP);
    if(sscanf(buf,"%le", &bulk->cr)!=1) COMPLAIN_EXIT;
    
    scanv3(buf,FP,scaling);
    LOGMSG(4,"scaling: %le %le %le",scaling[0],scaling[1],scaling[2]);
    scanv3(buf,FP,&bulk->tr[0][0]);
    LOGMSG(4,"tr1: %le %le %le",bulk->tr[0][0], bulk->tr[0][1], bulk->tr[0][2]);
    scanv3(buf,FP,&bulk->tr[1][0]);
    LOGMSG(4,"tr2: %le %le %le",bulk->tr[1][0], bulk->tr[1][1], bulk->tr[1][2]);
    scanv3(buf,FP,&bulk->tr[2][0]);
    LOGMSG(4,"tr3: %le %le %le",bulk->tr[2][0], bulk->tr[2][1], bulk->tr[2][2]);
    for (i=0; i<3; i++) {
        bulk->tr[1][i] *= scaling[i];
        bulk->tr[2][i] *= scaling[i];
        bulk->tr[3][i] *= scaling[i];
    }
        LOGMSG(4,"tr1: %le %le %le",bulk->tr[0][0], bulk->tr[0][1], bulk->tr[0][2]);
        LOGMSG(4,"tr2: %le %le %le",bulk->tr[1][0], bulk->tr[1][1], bulk->tr[1][2]);
        LOGMSG(4,"tr3: %le %le %le",bulk->tr[2][0], bulk->tr[2][1], bulk->tr[2][2]);
    
    bulk->site = (t_scsite *)malloc(sizeof(t_scsite)*bulk->nb);
    for (i=0;i<bulk->nb;i++) {
        site = &bulk->site[i];
        scanv3(buf,FP,site->coord);
        for (j=0; j<3; j++) {
            site->coord[j] *= scaling[j];
        }
        site->vmask=1.0;
        site->rmask=1.0;
    }
    fclose(FP);
    
    LOGMSG(2,"Reading inpch");
    FP = cat_openread("inpch", path);
    fgets(buf,400,FP);
    
    for (j=0; j<bulk->nb; j++) {
        site = &bulk->site[j];
        LOGMSG(3,"Processing site %d",j+1);
        fgets(buf,400,FP);
        if(sscanf(buf,"%d",&(site->nc))!=1) COMPLAIN_EXIT;
        LOGMSG(3, "Number of chemical elements on site: %d", site->nc);
        site->conc=(double *)calloc(site->nc,sizeof(double));
        site->at=(t_atom **)calloc(site->nc,sizeof(t_atom *));
        for (i=0; i<site->nc; i++){
            LOGMSG(3,"Loading %d-th concentration...",i+1);
            fgets(buf,400,FP);
            if(sscanf(buf,"%s",label)!=1) COMPLAIN_EXIT;
            LOGMSG(3,"Label = %s",label);
            fgets(buf,400,FP);
            if(sscanf(buf,"%le",(site->conc)+i)!=1) COMPLAIN_EXIT;
            LOGMSG(3,"Concentration = %le",site->conc[i]);
            strcat(label, "_B");
            LOGMSG(3,"Path = %s",path);
            site->at[i]=getatom_db(atoms, label, path);
        }
    }
    
    
    rawsr=0.0;
    maxwsr = 0.0;
    for(i=0; i<bulk->nb; i++) {
        rawsr += pow(bulk->site[i].at[0]->wsr,3);
        if(bulk->site[i].at[0]->wsr > maxwsr) maxwsr = bulk->site[i].at[0]->wsr;
    }
    bulk->maxwsr = maxwsr;
    LOGMSG(3,"max wsr = %le",maxwsr);
    rawsr = pow((rawsr/bulk->nb), (1.0/3.0));
    vol = getvolume(&bulk->tr[0][0], &bulk->tr[1][0], &bulk->tr[2][0]);
    LOGMSG(4,"vol: %le",vol);
    dawsr = pow((0.75*vol/(bulk->nb*M_PI)),(1.0/3.0));
    dimcf = rawsr/dawsr;
    LOGMSG(4,"factor: %le",dimcf);
    for (i=0; i<bulk->nb; i++) {
        for(j=0; j<3; j++) {
            bulk->site[i].coord[j] *= dimcf;
        }
    }
    for (i=0; i<3; i++) {
        bulk->tr[0][i] *= dimcf;
        bulk->tr[1][i] *= dimcf;
        bulk->tr[2][i] *= dimcf;
    }
    LOGMSG(4,"tr1: %le %le %le",bulk->tr[0][0], bulk->tr[0][1], bulk->tr[0][2]);
    LOGMSG(4,"tr2: %le %le %le",bulk->tr[1][0], bulk->tr[1][1], bulk->tr[1][2]);
    LOGMSG(4,"tr3: %le %le %le",bulk->tr[2][0], bulk->tr[2][1], bulk->tr[2][2]);
}

void copybulk(t_bulk **bout, t_bulk *bin) {
    // Carefull: doesn't really copy everything inside site
    int i,j;
    free(*bout);
    *bout = (t_bulk *)malloc(sizeof(t_bulk));
    (*bout)->nb = bin->nb;
    (*bout)->nmtr = bin->nmtr;
    (*bout)->cr = bin->cr;
    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            (*bout)->tr[i][j] = bin->tr[i][j];
        }
    }
    (*bout)->site = (t_scsite *)malloc(sizeof(t_scsite)*bin->nb);
    for (i=0; i<bin->nb; i++) {
        (*bout)->site[i] = bin->site[i];
    }
}

static void usage(void){
    printf("db2geom.x - input files generator for transport code\n");
    printf("Usage:\n"
        "\tdb2geom.x [options]\n"
        "\n");
    printf("Options:\n"
        "\t-d <path>    specify path to database folder\n"
        "\t-f           do not generate visualisation files\n"
        "\t-h           Display help.\n"
        "\t-t           Test run, does not write any files\n"
        "\t-v <num>     verbose output with <num> level of detail\n"
        "\n");
    exit(0);

}

    void readopts(int argc, char *argv[], t_xopt *xopt){
        int opt;
        char *st1,*st2;
        char buf[400];
        int fvibmask=0,frotmask=0,fperlayer=0;
        struct option long_options[] =
          {
            {"vmask", no_argument,       &fvibmask, 1},
            {"rmask",   no_argument,     &frotmask, 1},
            {"perlayer",   no_argument,     &fperlayer, 1},
            {0, 0, 0, 0}
          };
        int option_index = 0;

        xopt->vis=1;
        xopt->db=0;
        xopt->dry=0;
        while ((opt=getopt_long(argc,argv,"+v:d:fht",long_options, &option_index))!=-1){ 
            switch(opt){
                
                case 'd':
                xopt->db=1;
                xopt->dbloc=(char *)malloc(sizeof(char)*(strlen(optarg)+1));
                strcpy(xopt->dbloc, optarg);
                break;
 
                case 'f':
                xopt->vis=0;
                break;
                
                case 't':
                xopt->dry=1;
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

                default:
                usage();
                exit(0);
            }
        }
        
        LOGMSG(1,"Finished options processing");
    }
    
    
    int read_input_ml(t_multilay *ml){    
        FILE *FP;
        char buf[400], label[30];
        t_slab *slab;
        char *pos;
        int nopen;
        t_varinput *varsin, *prev;

        LOGMSG(2,"Reading multilay_input...");

        FP=cat_openread("multilay_input", "");
        fgets(buf,400,FP);
        fgets(buf,400,FP);
        if(sscanf(buf,"%d",&(ml->ns))!=1) COMPLAIN_EXIT;
        LOGMSG(3, "Number of slabs: %d", ml->ns);
        fgets(buf,400,FP);
        if(sscanf(buf,"%dx%d", &(ml->sc[0]), &(ml->sc[1]))!=2) COMPLAIN_EXIT;
        LOGMSG(3, "Supercell: %d x %d", ml->sc[0], ml->sc[1]);
        fgets(buf,400,FP);
        if(sscanf(buf,"%d", &(ml->nconf))!=1) COMPLAIN_EXIT;
        LOGMSG(3,"Storing parameter input...");
        fgets(buf,400,FP);
        ml->parsin = (t_varinput *)malloc(sizeof(t_varinput));
        ml->npar = 0;
        varsin = ml->parsin;
        prev = varsin;
        while (count_char(buf, '{') == count_char(buf, '}')) {
            if (count_char(buf, '=') > 0 ) {
                store_var(buf, varsin);
                prev = varsin;
                varsin = varsin->next;
                ml->npar++;
            }
            fgets(buf,400,FP);
        }
        prev->next = NULL;
        free(varsin);
        
        ml->slab=(t_slab *)malloc(sizeof(t_slab)*ml->ns);
        LOGMSG(3,"Storing slab info...");
        for (slab=ml->slab;slab<(ml->slab+ml->ns);slab++){
            LOGMSG(3,"Slab: %d",slab-ml->slab+1);
            
            pos = strchr(buf, '{');
            nopen = 0;
            if (pos != NULL) {
                *pos = '\0';
                nopen = 1;
            } else {
                COMPLAIN_EXIT;
            };
            LOGMSG(5,"%s","Found {");
            if(sscanf(buf,"%s",label)!=1) COMPLAIN_EXIT;
            slab->material=(char *)malloc(sizeof(char)*(strlen(label)+1));
            strcpy(slab->material,label);
            LOGMSG(3,"Material: %s",slab->material);
            slab->varsin = (t_varinput *)malloc(sizeof(t_varinput));
            varsin = slab->varsin;
            prev = varsin;
            while ( nopen > 0 ){
                if ( fgets(buf,400,FP) == NULL) COMPLAIN_EXIT;
                nopen = nopen + count_char(buf, '{') - count_char(buf, '}');
                if (nopen > 0 || nopen == 0 && count_char(buf, '=') > 0) {
                    store_var(buf, varsin);
                    prev = varsin;
                    varsin = varsin->next;
                }
            };
            prev->next = NULL;
            free(varsin);
            fgets(buf,400,FP);
        };
        return(1);
    }

    
    void store_var(char *buf, t_varinput *varsin){
        int n, n2;
        char *pos = strchr(buf,'=');
        char str[MAXSTRLEN];
        
        if (pos == NULL) COMPLAIN_EXIT;
        n = pos-buf;
        if(sscanf(buf," %[^= ]",str)!= 1) COMPLAIN_EXIT;
        strallcpy(&((varsin)->name), str);
        LOGMSG(5, "Name: %s", (varsin)->name);
        n2 = strlen(buf)-n-1;
        strallcpy(&((varsin)->value), pos+1);
        LOGMSG(5, "Value string: %s", (varsin)->value);
        (varsin)->nread = 0;
        varsin->next = (t_varinput *)malloc(sizeof(t_varinput));
    }
    
    

    void get_db_dir(char *argv[], t_xopt *xopt){
        char *xloc, *dbloc, *pos, *p2;
        char path[400], dest[400];
        pid_t pid = getpid();
        struct stat s;
        int len, i;

        if (xopt->db) {
            pos = strstr(xopt->dbloc, dbdirname);
            if (pos == NULL) {
                dbloc = (char *)malloc(sizeof(char)*(strlen(xopt->dbloc)+strlen(dbdirname)+2));
                strcpy(dbloc, xopt->dbloc);
                strcat(dbloc, "/");
                strcat(dbloc, dbdirname);
            } else {
                dbloc = xopt->dbloc;
            };
            if(!test_dir(dbloc)) {
                dbloc = ask_dbloc();
            };
        } else {
            sprintf(path, "/proc/self/exe");
            if (readlink(path, dest, 400) != -1){
                //printf("proc works\n");
                xloc = (char *)malloc(sizeof(char)*(strlen(dest)+1));
                strcpy(xloc, dest);
            } else {    
                //printf("proc doesn't work...\n");
                xloc = (char *)malloc(sizeof(char)*(strlen(argv[0])+1));
                strcpy(xloc, argv[0]);
            };
            
            p2 = xloc;
            while (p2 != NULL) {
                pos = p2;
                p2 = strstr(pos+4, "/bin/");
            }
            len = pos - xloc + 1;
            dbloc = (char *)malloc(sizeof(char)*(len+strlen(dbdirname)+1));
            strncpy(dbloc, xloc, len);
            strcpy(dbloc+len, dbdirname);
            if (!test_dir(dbloc)) {
                free(dbloc);
                dbloc = ask_dbloc();
            }
            free(xloc);
        };
        if (!test_dir(dbloc)) {
            printf("%s = ", dbloc);
            printf("Still wrong.\n");
            exit(1);
        }
        xopt->dbloc = dbloc;
        LOGMSG(1,"Found potentials database");
    }

    int test_dir(char *dir){
        struct stat s;
        int flag;
        flag = 0;
        if( stat(dir, &s) >= 0) {
            if(S_ISDIR(s.st_mode)) {
                flag = 1;
            }
        }
        return(flag);
    }

    int test_file(char *file){
        struct stat s;
        int flag;
        flag = 0;
        if( stat(file, &s) >= 0) {
            flag = 1;
        }
        return(flag);
    }

    char * ask_dbloc() {
        char buf[400], *dbloc;
        char *pos;
        printf("%s %s %s","Cannot find the location of the",dbdirname,"folder.\n");
        printf("Enter path:  ");
        scanf("%s", buf);
        pos = strstr(buf, dbdirname);
        if (pos == NULL) {
            strcat(buf, "/");
            strcat(buf, dbdirname);
        };
        dbloc = (char *)malloc(sizeof(char)*(strlen(buf)+1));
        strcpy(dbloc, buf);
        return(dbloc);
    }
    
    FILE *cat_openread(char *name, char *path) {
        char *full;
        FILE *FP;
        full = (char *)calloc(strlen(name)+strlen(path)+2, sizeof(char));
        if (strlen(path)>0) {
            strcpy(full, path);
            strcat(full, "/");
        }
        strcat(full, name);
        LOGMSG(3,"Opening %s", full);
        FP = fopen(full, "r");
        if(FP==NULL) {
            EXITMSG("Error opening %s", full);
        }
        free(full);
        return(FP);
    }
    
    void init_dbe(t_db_entry *dbe) {
        int i;
        switch(dbe->kind) {
            case 1:
            dbe->nvar = sizeof(dbe1_props)/sizeof(dbe1_props[0]);
            dbe->vp = (t_varprops *)malloc(sizeof(t_varprops)*dbe->nvar);
            for (i=0; i<dbe->nvar; i++) {
                dbe->vp[i] = dbe1_props[i];
            }
            break;
            
            default:
            EXITMSG("Don't know database entry kind %d", dbe->kind);
            break;
        }
    }
    
    int read_db_entry(t_db_entry *dbe, char *name, char *dbloc) {
        FILE *f;
        char *path;
        char buf[400], str[MAXSTRLEN], str1[MAXSTRLEN], str2[MAXSTRLEN];
        int kind;
        t_varprops *vp;
        
        LOGMSG(3,"%s %s","Reading db entry for", name);
        path = (char *)malloc(sizeof(char)*(strlen(name)+strlen(dbloc)+2));
        strcpy(path, dbloc);
        strcat(path, "/");
        strcat(path, name);
        
        f = cat_openread("minmax", path);
        
        strallcpy(&(dbe->name), name);
        fgets(buf,400,f);
        if(sscanf(buf, "%d", &dbe->kind) != 1) COMPLAIN_EXIT;
        init_dbe(dbe);
        while (fgets(buf,400,f) != NULL) {
            if(sscanf(buf,"%s %s %s", str, str1, str2) != 3) COMPLAIN_EXIT;
            for (vp=dbe->vp; vp<dbe->vp+dbe->nvar; vp++) {
                if (strcmp(str, vp->name) == 0) {
                    LOGMSG(4,"new min and max for %s", vp->name);
                    str2val(&(vp->min), str1, vp->type);
                    str2val(&(vp->max), str2, vp->type);
                }
            }
        }
        fclose(f);
        switch(dbe->kind) {
            case 1:
            read_bulk_turek(&dbe->bulk, path, &dbe->atoms);
            break;
            
            default:
            EXITMSG("Database entry kind undefined");
            break;
        }
        free(path);
        return(1);
    }
    
    int read_db (t_db_entry **dbe, t_multilay *ml, char *dbloc){
        t_slab *slab;
        t_db_entry **p, *dbe0, **dbep;
        int nb, i, j, k;
        t_atom *ptr, **pp;
        t_scsite *site;
        
        LOGMSG(1,"Reading database entries");
        dbep = dbe;
        pp=&ml->atoms.ptr;
        for (slab=ml->slab; slab<(ml->slab+ml->ns); slab++) {
            LOGMSG(3,"\n");
            LOGMSG(3,"Checking slab: %d",slab-ml->slab+1);
            p = dbe;
            while (*p != NULL) {
                LOGMSG(4,"%s = %s ?", slab->material, (*p)->name);
                if (strcmp(slab->material, (*p)->name) == 0) {
                    LOGMSG(4,"Yes: no new entry.");
                    slab->dbe = *p;
                    break;
                }
                p = &(*p)->next;
            }
            if (*p == NULL) {
                *dbep = dbe0 = (t_db_entry *)calloc(1, sizeof(t_db_entry));
                read_db_entry(dbe0, slab->material, dbloc);
                slab->dbe = dbe0;
                for( *pp = dbe0->atoms.ptr; *pp!=NULL ;pp=(t_atom **)&(*pp)->next);
                ml->atoms.num += dbe0->atoms.num;
                LOGMSG(4,"Increased number of atoms to %d", ml->atoms.num);
                dbep = &dbe0->next;
            }
        }
        // determine the expectations
        for (dbe0=*dbe; dbe0!=NULL; dbe0=dbe0->next) {
            switch(dbe0->kind) {
                case 1:
                site=dbe0->bulk.site;
                nb=dbe0->bulk.nb;
                break;
            }
            //dbe0->natexp = (double *)calloc(ml->atoms.num, sizeof(double));
            alloc_basket(&dbe0->basket, ml->atoms.num);
            for (j=0; j<nb; j++, site++) {
                for (i=0; i<site->nc; i++) {
                    //for (ptr=ml->atoms.ptr, k=0; k<ml->atoms.num; ptr=ptr->next, k++) {
                    //    if(site->at[i] == ptr) break;
                    //}
                    k = getatind(site->at[i], &ml->atoms);
                    dbe0->basket.initial[k] += site->conc[i];
                }
            }
            for (k=0; k<ml->atoms.num; k++) {
                LOGMSG(4,"Expectation atom %d: %g", k, dbe0->basket.initial[k]);
            }
        }
        return(1);
    }
    
    int read_vars(t_multilay *ml) {
        t_slab *slab;
        t_varinput *vin, *pin;
        t_varprops *vp;
        t_var *var;
        t_par *par, *p2;
        t_par **prevnext=&ml->parptr;
        int i, j, n, ndots, i1,i2,i3;
        double d1, d2, d3;
        char *pos, str[MAXSTRLEN];
        
        LOGMSG(1,"\n");
        LOGMSG(1,"Initializing parameters");
        pin = ml->parsin;
        ml->par = (t_par *)calloc(ml->npar, sizeof(t_par));
        for(par = ml->par; par-ml->par < ml->npar; par++) {
            LOGMSG(3,"%s", pin->name);
            if (pin->name[1] == '\0' && isalpha(pin->name[0])) {
                par->symbol = pin->name[0];
                for (p2 = ml->par; p2<par; p2++) {
                    if(p2->symbol == par->symbol) {
                        EXITMSG("%s","Found parameters with identical name!");
                    }
                }
            } else {
                EXITMSG("%s","Use 1 letter for parameters!");
            }
            pin = pin->next;
        }
        
        LOGMSG(1,"Reading input variables");
        for (slab=ml->slab; slab<(ml->slab+ml->ns); slab++) {
            LOGMSG(2,"Slab: %d",slab-ml->slab+1);
            slab->nvar = slab->dbe->nvar;
            slab->var = (t_var *)malloc(sizeof(t_var)*slab->nvar);
            vin = slab->varsin;
            vp = slab->dbe->vp;
            var = slab->var;
            while( vp - slab->dbe->vp < slab->dbe->nvar ) {
                var->vp = vp;
                var->nset = 0;
                vp++;
                var++;
            }
            while(vin != NULL) {
                var = slab->var;
                while( var - slab->var < slab->nvar ) {
                    if (strcmp(var->vp->name, vin->name) == 0) {
                        LOGMSG(3,"%s",var->vp->name);
                        vin->nread = read_val(var, vin->value, ml);
                        break;
                    }
                    var++;
                }
                
                vin = vin->next;
            }
            vin = slab->varsin;
            while(vin != NULL) {
                if(!vin->nread) {
                    LOGMSG(0,"Don't know how to read; ignored \'%s\'",vin->name);
                }
                vin = vin->next;
            }
        }
        
        LOGMSG(1,"Reading parameter values");
        pin = ml->parsin;
        for(i=0; i < ml->npar; i++) {
            par = &ml->par[i];
            if (par->type == 0) {
                *prevnext = i+1==ml->npar ? NULL : par+1;
            } else {
                *prevnext = par;
                prevnext = &par->next;
                while (pin != NULL) {
                    LOGMSG(5, "Comparing %c to %s", par->symbol, pin->name);
                    if(par->symbol == pin->name[0]) {
                        LOGMSG(2, "Reading %c...", par->symbol);
                        ndots = count_str(pin->value, "..");
                        pos = posafter(pin->value);
                        switch(ndots) {
                        
                            case 0:
                            n = count_char(pin->value, ',') + 1;
                            par->val = (t_value *)malloc(sizeof(t_value)*n);
                            par->nval = n;
                            for(j=0; j<n; j++) {
                                if(sscanf(pos," %[^ ,\t\n}]",str)!= 1) COMPLAIN_EXIT;
                                str2val(&(par->val[j]), str, par->type);
                                pos = strchr(pos, ',')+1;
                            }
                            break;
                        
                            case 1:
                            if(count_char(pin->value, ',') > 0) COMPLAIN_EXIT;
                            if(par->type != 1) {
                                EXITMSG("%c: Cannot process a..b input for non-integers, use a .. b .. c", par->symbol);
                            }
                            if(sscanf(pos," %d .. %d",&i1,&i2)!= 2) COMPLAIN_EXIT;
                            n = i2-i1+1;
                            par->val = (t_value *)malloc(sizeof(t_value)*n);
                            par->nval = n;
                            for(j=0;j<n;j++) {
                                par->val[j].i = i1+j;
                                LOGMSG(5,"%d", par->val[j].i);
                            }
                            break;
                            
                            case 2:
                            switch(par->type) {
                                case 1:
                                if(sscanf(pos," %d .. %d .. %d",&i1,&i2,&i3)!= 3) COMPLAIN_EXIT;
                                n = (i3-i1)/i2 + 1;
                                par->val = (t_value *)malloc(sizeof(t_value)*n);
                                par->nval = n;
                                for(j=0;j<n;j++) {
                                    par->val[j].i = i1+j*i2;
                                    LOGMSG(5,"%d", par->val[j].i);
                                }
                                break;
                                
                                case 2:
                                if(sscanf(pos," %le .. %le .. %le",&d1,&d2,&d3)!= 3) COMPLAIN_EXIT;
                                n = floor((d3-d1)/d2) + 1;
                                par->val = (t_value *)malloc(sizeof(t_value)*n);
                                par->nval = n;
                                for(j=0;j<n;j++) {
                                    par->val[j].d = d1+j*d2;
                                    LOGMSG(5,"%lg", par->val[j].d);
                                }
                                break;
                            }
                            break;
                            
                            default:
                            EXITMSG("I'm confused about all the double dots");
                            break;
                        }
                        for(j=0;j<n;j++) {
                            LOGMSG(5,"VAL MIN MAX : %g %g %g",par->val[j].d, par->min.d, par->max.d);
                            if(test_minmax(par->val[j],par->min,par->max,par->type)) {
                                EXITMSG("min or max of %c violated", par->symbol);
                            }
                        }
                        break;
                    }
                    pin = pin->next;
                }
            }
        }
        return(1);
    }
    
    int read_val(t_var *var, char *valstr, t_multilay *ml) {
        int i, j, n, flag;
        char str[MAXSTRLEN];
        char *pos;
        t_valptr **ptr;
        t_par *par;
        
        LOGMSG(3,"Reading %s...",var->vp->name);
        n = count_char(valstr, ',') + 1;
        //allocval(&(var->val), var->vp->type, n);
        var->val = (t_value *)malloc(sizeof(t_value)*n);
        pos = strchr(valstr, '{');
        if (pos != NULL) {
            pos++;
        } else {
            pos = valstr;
        }
        for(i=0; i<n; i++) {
            if(sscanf(pos," %[^ ,\t\n}]",str)!= 1) COMPLAIN_EXIT;
            if(isalpha(str[0])) {
                if (str[1] == '\0') {
                    LOGMSG(3,"Linking parameter %s",str);
                   //check parameter
                    flag = 1;
                    for (j=0;j<ml->npar;j++) {
                        par = &ml->par[j];
                        if(par->symbol == str[0]) {
                            flag = 0;
                            par->nptr++;
                            ptr = &(par->ptr);
                            while (*ptr != NULL) {
                                if(var->vp->type != (*ptr)->vp->type) {
                                    EXITMSG("Data type of parameter %s does not match in all occurrences", str);
                                }
                                ptr = &(*ptr)->next;
                            }
                            *ptr = (t_valptr *)calloc(1, sizeof(t_valptr));
                            (*ptr)->val = &var->val[i];
                            (*ptr)->vp = var->vp;
                            switch(var->vp->type) {
                                case 1:
                                if (par->type == 0) {
                                    par->type = var->vp->type;
                                    par->min.i = var->vp->min.i;
                                    par->max.i = var->vp->max.i;
                                } else {
                                    if(par->min.i < var->vp->min.i) par->min.i = var->vp->min.i;
                                    if(par->max.i > var->vp->max.i) par->max.i = var->vp->max.i;
                                }
                                break;
                                case 2:
                                if (par->type == 0) {
                                    par->type = var->vp->type;
                                    par->min.d = var->vp->min.d;
                                    par->max.d = var->vp->max.d;
                                } else {
                                    if(par->min.d < var->vp->min.d) par->min.d = var->vp->min.d;
                                    if(par->max.d > var->vp->max.d) par->max.d = var->vp->max.d;
                                }
                                break;
                            }
                        }
                    }
                    if (flag) EXITMSG("Parameter %s not defined", str);
                    
                } else {
                    EXITMSG("%s","Use 1 letter for parameters!");
                }
                
            } else {
                str2val(&(var->val[i]), str, var->vp->type);
                if(test_minmax(var->val[i],var->vp->min,var->vp->max,var->vp->type)) {
                    EXITMSG("min or max of %s violated", var->vp->name);
                }
            }
            pos = strchr(pos, ',')+1;
        }
        if (n != var->vp->n) EXITMSG("Number of elements in %s not correct.", var->vp->name);
        var->nset = n;
        return(n);
    }
    
    void bulk2layer(t_slab *slab) {
        t_var *var;
        t_bulk *bulk, *bulk2=NULL;
        t_layer *lay;
        t_atom *ptr;
        int v1[3], v2[3], v3[3];
        double cross[3], v1a, crabs;
        double vec1[3] = {0,0,0};
        double vec2[3] = {0,0,0};
        double vec3[3] = {0,0,0};
        double m1[3][3], m2[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        double m1inv[3][3], trans[3][3];
        double coord[3], orig[3];
        double vol, rr, cut, dist;
        int i, j, k, n, r, ib, ii;
        
        LOGMSG(2,"Making layer geometry for %s", slab->material);
        bulk = &slab->dbe->bulk;
        lay = &slab->layer;
        
        rr=0.0;
        for(i=0; i<bulk->nb; i++) {
            rr += pow(bulk->site[i].at[0]->wsr,3);
        }
        rr *= M_PI/0.75;
        
        var = getvar(slab, "scvec1", 0);
        for (i=0; i<3; i++) {
            v1[i] = var->val[i].i;
            for (j=0; j<3; j++) {
                vec1[j] += v1[i]*bulk->tr[i][j];
            }
        }
        var = getvar(slab, "scvec2", 0);
        for (i=0; i<3; i++) {
            v2[i] = var->val[i].i;
            for (j=0; j<3; j++) {
                vec2[j] += v2[i]*bulk->tr[i][j];
            }
        }
        var = getvar(slab, "scvec3", 0);
        for (i=0; i<3; i++) {
            v3[i] = var->val[i].i;
            for (j=0; j<3; j++) {
                vec3[j] += v3[i]*bulk->tr[i][j];
            }
        }
        
        v1a = abs3(vec1);
        cross3(cross, vec1, vec2);
        crabs = abs3(cross);
        vol = inp3(cross, vec3);
        if ( vol > 0 ) {
            k = 1;
        } else if (vol < 0) {
            k = -1;
            } else {
            EXITMSG("wrong vectors!");
        }
        lay->vol = fabs(vol);
        //m2 = { {v1a, inp3(vec1,vec2)/v1a, 0}, {0, crabs/v1a ,0} , {0,0, crabs } };
        lay->tr1[0] = m2[0][0] = v1a;
        lay->tr1[1] = 0;
        lay->tr2[0] = m2[0][1] = inp3(vec1,vec2)/v1a;
        lay->tr2[1] = m2[1][1] = k*crabs/v1a;
        m2[2][2] = k*crabs;
        for (i=0; i<3; i++) {
            m1[i][0] = vec1[i];
            m1[i][1] = vec2[i];
            m1[i][2] = cross[i];
        };
        
        inverse3(m1inv, m1);
        mmult3(trans, m2, m1inv);
        LOGMSG(4,"Found transformation:");
        LOGMSG(4,"%le %le %le", trans[0][0], trans[0][1], trans[0][2]);
        LOGMSG(4,"%le %le %le", trans[1][0], trans[1][1], trans[1][2]);
        LOGMSG(4,"%le %le %le", trans[2][0], trans[2][1], trans[2][2]);
        
        LOGMSG(4,"layer tr1: %le %le", lay->tr1[0], lay->tr1[1]);
        LOGMSG(4,"layer tr2: %le %le", lay->tr2[0], lay->tr2[1]);
        vmult3(lay->trperp, trans, vec3);
        LOGMSG(3,"layer trperp: %le %le %le", lay->trperp[0], lay->trperp[1], lay->trperp[2]);
        
        copybulk(&bulk2, bulk);
        
        vmult3(&bulk2->tr[0][0], trans, &bulk->tr[0][0]);
        LOGMSG(4,"transformed bulk tr1: %le %le %le", bulk2->tr[0][0], bulk2->tr[0][1], bulk2->tr[0][2]);
        vmult3(&bulk2->tr[1][0], trans, &bulk->tr[1][0]);
        LOGMSG(4,"transformed bulk tr2: %le %le %le", bulk2->tr[1][0], bulk2->tr[1][1], bulk2->tr[1][2]);
        vmult3(&bulk2->tr[2][0], trans, &bulk->tr[2][0]);
        LOGMSG(4,"transformed bulk tr3: %le %le %le", bulk2->tr[2][0], bulk2->tr[2][1], bulk2->tr[2][2]);
        for (i=0; i<bulk->nb; i++) {
            vmult3(coord, trans, bulk->site[i].coord);
            bulk2->site[i].coord[0]=coord[0];
            bulk2->site[i].coord[1]=coord[1];
            bulk2->site[i].coord[2]=coord[2];
            LOGMSG(4,"transformed bulk coord: %le %le %le", coord[0], coord[1], coord[2]);
        }
        
        n = round(fabs(vol)/rr);
        LOGMSG(3,"%d bulk cells in 1 layer", n);
        lay->ncell = n;
        lay->nb = n*bulk2->nb;
        
        lay->site = (t_scsite *)realloc(lay->site, sizeof(t_scsite)*lay->nb);
        
        LOGMSG(3,"Checking which sites are inside the layer cell");
        zero3(vec1); zero3(vec2); zero3(vec3);
        vec1[0] = lay->tr1[0];
        vec2[0] = lay->tr2[0];
        vec2[1] = lay->tr2[1];
        vec3[0] = lay->trperp[0];
        vec3[1] = lay->trperp[1];
        vec3[2] = lay->trperp[2];
        lay->volbelow = 0.0;
        // Simple way to loop over translations outwards from the origin:
        for(n=0, r=0; r<5/*n<lay->nb*/; r++) {
            for(i=-r; i<=r; i++) {
                for(j=-r; j<=r; j++) {
                    for(k=-r; k<=r; k++) {
                        if (abs(i)+abs(j)+abs(k) == r) {
                            //LOGMSG(5,"Checking translation: %d, %d, %d",i,j,k);
                            orig[0] = i*bulk2->tr[0][0] + j*bulk2->tr[1][0] + k*bulk2->tr[2][0];
                            orig[1] = i*bulk2->tr[0][1] + j*bulk2->tr[1][1] + k*bulk2->tr[2][1];
                            orig[2] = i*bulk2->tr[0][2] + j*bulk2->tr[1][2] + k*bulk2->tr[2][2];
                            // look for atoms inside layer cell
                            for(ib=0; ib<bulk2->nb; ib++) {
                                for(ii=0; ii<3; ii++) {
                                    coord[ii] = orig[ii] + bulk2->site[ib].coord[ii];
                                }
                                //LOGMSG(5,"coord: %g, %g, %g",coord[0],coord[1],coord[2]);
                                
                                if (isinside(coord, vec1, vec2, vec3)) {
                                    LOGMSG(4,"inside: %d, %d, %d",i,j,k);
                                    LOGMSG(4,"coord: %le %le %le", coord[0], coord[1], coord[2]);
                                    // copy site
                                    memcpy(&lay->site[n], &bulk2->site[ib], sizeof(t_scsite));
                                    for(ii=0; ii<3; ii++) {
                                        lay->site[n].coord[ii] = coord[ii];
                                    }
                                    rr = lay->site[n].at[0]->wsr;
                                    dist = rr - coord[2];
                                    if (dist > 0) {
                                        lay->volbelow += M_PI*dist*dist*(rr-dist/3.0);
                                    }
                                    n++;
                                    //LOGMSG(0,"HERE %d", n);
                                }
                            }
                        }
                    }
                }
            }
        }
        LOGMSG(3,"found %d atoms inside", n);
        if(n>lay->nb) EXITMSG("Too many atoms inside; something wrong with coordinates.");
        lay->cr = bulk->cr;
        lay->maxwsr = bulk->maxwsr;
        //set nmtr
        cut = lay->cr*lay->maxwsr;
        cross3(cross, vec1, vec3);
        crabs = abs3(cross);
        dist = fabs(inp3(cross, vec2)/crabs);
        lay->nmtr = ceil(cut/dist);
        cross3(cross, vec2, vec3);
        crabs = abs3(cross);
        dist = fabs(inp3(cross, vec1)/crabs);
        n = ceil(cut/dist);
        if (n > lay->nmtr) lay->nmtr = n;
        LOGMSG(4,"nmtr set to: %d", lay->nmtr);
        lay->tp = 1; //this one will be ignored anyway
        
        free(bulk2->site);
        free(bulk2);
        LOGMSG(2,"Done\n");
    }
    
    void stretch(t_layer *lay, double *vec1, double *vec2) {
        double angle, area, x, z, r, d;
        int ib;
        
        angle = angle2(vec1, vec2);
        area = cross2(vec1, vec2);
        x = vec1[0]/lay->tr1[0];
        printf("Stretching by %g %%\n", (x-1)*100);
        LOGMSG(5,"angle diff: %le", angle-angle2(lay->tr1, lay->tr2));
        LOGMSG(5,"area diff: %le", area-x*x*cross2(lay->tr1, lay->tr2));
        if (fabs(1-angle2(lay->tr1, lay->tr2)/angle) > accuracy 
                    || fabs(1-x*x*cross2(lay->tr1, lay->tr2)/area) > accuracy ) {
            EXITMSG("Shapes of cells do not match (in xy-plane).");
        }
        z = 1.0/(x*x);
        lay->tr1[0] = vec1[0];
        lay->tr1[1] = vec1[1];
        lay->tr2[0] = vec2[0];
        lay->tr2[1] = vec2[1];
        lay->trperp[0] *= x;
        lay->trperp[1] *= x;
        lay->trperp[2] *= z;
        lay->volbelow = 0.0;
        for(ib=0; ib<lay->nb; ib++) {
            lay->site[ib].coord[0] *= x;
            lay->site[ib].coord[1] *= x;
            lay->site[ib].coord[2] *= z;
            // recalculate volume of sphere caps outside the cell:
            r = lay->site[ib].at[0]->wsr;
            d = r - lay->site[ib].coord[2];
            if (d > 0) {
                lay->volbelow += M_PI*d*d*(r-d/3.0);
            }
        }
        LOGMSG(4, "volbelow/vol = %g", lay->volbelow/lay->vol);
    }
    
    t_atom *pick_atom(t_scsite *site, t_basket *basket, t_atomset *atoms) {
        double p1, thres[site->nc+1];
        int ic, ia[site->nc], nc=site->nc;
        double number=drand48();
        if (basket->nleft == 0) reset_basket(basket);
        thres[0] = 0.0;
        for (ic=0; ic<nc; ic++) {
            ia[ic] = getatind(site->at[ic], atoms);
            p1 = basket->natexp[ia[ic]]/(basket->initial[ia[ic]]*basket->nfix);
            thres[ic+1] = thres[ic] + site->conc[ic]*p1;
        }
        for (ic=0; ic<nc; ic++) {
            thres[ic+1] /= thres[nc];
            //LOGMSG(6,"thres(%d) = %g",ic,thres[ic+1]);
            if (number < thres[ic+1]) {
                basket->nleft--;
                basket->natexp[ia[ic]] -= 1.0;
                basket->natexp[ia[ic]] *= basket->natexp[ia[ic]] > 0.0;
                //LOGMSG(6,"Picked: %d",ia[ic]);
                return(site->at[ic]);
            }
        }
        return(NULL);
    }
    
    void extrafields(t_trgeom *geom){
        int i,j;
        double d=geom->dimcf*5.2917721092e-1;
        
        for(i=0;i<geom->na;++i){
            for (j=0;j<3;++j){
                geom->rcoord[i][j]=geom->site[i].coord[j]*geom->scx[j]*d;            
            }
        }    
    }
    
    void prep_conf(t_config *conf, t_multilay *ml) {
        t_trgeom *rg, *lg, *mgeom;
        t_slab *slab;
        t_db_entry *dbe;
        t_var *var;
        t_layer *lay;
        t_scsite *laysite;
        t_trsite *site;
        double stot=0.0, area, volin=0.0;
        double vec1[2] = {0,0};
        double vec2[2] = {0,0};
        double off1[3]={0,0,0}, off2[2], shift[ml->ns][2];
        int nlay[ml->ns], *nfix;
        int i, j, k, l, n, m, ncell, t, newnmtr, p;
        t_basket *basket;
        double minvav, rms;
        //int sidx;
        
        rg = &conf->rg;
        lg = &conf->lg;
        mgeom = &conf->mg;
        
        LOGMSG(2,"Preparing configurations...\n");
        
        conf->na = (int *)malloc(sizeof(int)*ml->ns);
        conf->basket = (t_basket **)malloc(sizeof(t_basket*)*ml->ns);
        conf->ns = ml->ns;
        
        conf->nfix = (int *)malloc(sizeof(int)*ml->ns);
        conf->diskind = (int *)malloc(sizeof(int)*ml->ns);
        nfix = conf->nfix;
        
        mgeom->na = 0;
        mgeom->tp = 1;
        mgeom->nmtr = 0;
        mgeom->sc[0]=mgeom->sc[1]=1;
        mgeom->cr = 0.0;
        mgeom->dimcf = 1.0;
        mgeom->scx[0]=mgeom->scx[1]=mgeom->scx[2]=1.0;
        
        conf->nint= -1;
        for (i=0; i<ml->ns; i++) {
            slab = &ml->slab[i];
            switch (slab->dbe->kind) {
                case 1:
                bulk2layer(slab);
                slab->dbe->basket.npercell = slab->dbe->bulk.nb;
                conf->nint += 1; // nint has to be done smarter for self-consistent interfaces
                break;
            }
            // get stiffness
            var = getvar(slab, "stiff", 1);
            slab->stiff = (var == NULL ? 0.0 : var->val[0].d);
            stot += slab->stiff;
            vec1[0] += slab->layer.tr1[0]*slab->stiff;
            vec1[1] += slab->layer.tr1[1]*slab->stiff;
            vec2[0] += slab->layer.tr2[0]*slab->stiff;
            vec2[1] += slab->layer.tr2[1]*slab->stiff;
            
            // get number of layers
            var = getvar(slab, "nlay", 0);
            nlay[i] = var->val[0].i;
            
            mgeom->na += slab->layer.nb*nlay[i];
            if (slab->layer.cr > mgeom->cr) mgeom->cr = slab->layer.cr;
            if (slab->layer.nmtr > mgeom->nmtr) mgeom->nmtr = slab->layer.nmtr;
        }
        mgeom->na *= ncell = ml->sc[0]*ml->sc[1];
        conf->coord = (double (*)[3])malloc(sizeof(double[3])*mgeom->na);
        conf->laysite = (t_scsite **)malloc(sizeof(t_scsite*)*mgeom->na);
        conf->rms = (double *)malloc(sizeof(double)*mgeom->na);
        conf->zint = (double *)malloc(sizeof(double)*conf->nint);
        //conf->sidx = (int *)malloc(sizeof(int)*mgeom->na);
        LOGMSG(2,"%d atoms in scattering region", mgeom->na);
        
        t = ml->sc[0] < ml->sc[1] ? ml->sc[0] : ml->sc[1];
        
        if (stot == 0.0) {
            EXITMSG("Sum of stiffnesses = 0, cannot handle that.");
        }
        mgeom->trv[0][0] = ml->sc[0]*(vec1[0] /= stot);
        mgeom->trv[0][1] = ml->sc[0]*(vec1[1] /= stot);
        mgeom->trv[1][0] = ml->sc[1]*(vec2[0] /= stot);
        mgeom->trv[1][1] = ml->sc[1]*(vec2[1] /= stot);
        LOGMSG(3,"Calculated average tr1: %g %g", vec1[0], vec1[1]);
        LOGMSG(3,"Calculated average tr2: %g %g\n", vec2[0], vec2[1]);
        
        mgeom->site = (t_trsite *)malloc(sizeof(t_trsite)*mgeom->na);
        area = cross2(vec1, vec2);
        
        // preparing middle region (atoms are not yet determined)
        n=0;
        l=0;
        p=0;
        //sidx=0;
        while (1) {
            slab = &ml->slab[l];
            //sidx++;
            lay = &slab->layer;
            dbe = slab->dbe;
            basket = &dbe->basket;
            conf->basket[l] = basket;
            LOGMSG(2,"Preparing slab %d ...", l+1);
            
            if (lay->cr > mgeom->cr) mgeom->cr = lay->cr;
            
            newnmtr = ceil(lay->nmtr/(double)t);
            if (newnmtr > mgeom->nmtr) mgeom->nmtr = newnmtr;
            
            var = getvar(slab, "nfixc", 1);
            if (var != NULL) {
                nfix[l] = var->val[0].i;
            } else {
                nfix[l] = nlay[l]*ncell*lay->ncell;
            }
            basket->nfix = nfix[l];
            reset_basket(basket);
            for (i=0; i<basket->num; i++) {
                if(!check_round(basket->natexp[i])) {
                    LOGMSG(0,"WARNING: nfixc not consistent with concentrations; not fixing exactly");
                    break;
                }
            }
            
            var = getvar(slab, "shift", 1);
            if (var != NULL) {
                shift[l][0] = var->val[0].d*vec1[0] + var->val[1].d*vec2[0];
                shift[l][1] = var->val[0].d*vec1[1] + var->val[1].d*vec2[1];
            } else {
                shift[l][0]=shift[l][1]=0.0;
            }
            
            var = getvar(slab, "rms", 1);
            if (var != NULL) {
                minvav = 0.0;
                for(i=0; i<lay->nb; i++) {
                    for(j=0; j<lay->site[i].nc; j++) {
                        minvav += lay->site[i].conc[j]/lay->site[i].at[j]->mass;
                    }
                }
                minvav /= (double)lay->nb;
                // rms per inverse square root atomic mass unit:
                rms = var->val[0].d/sqrt(3.0*minvav);
                conf->diskind[l] = 1;
            } else {
                rms = 0.0;
                conf->diskind[l] = 0;
            }

            stretch(lay, vec1, vec2);
            
            //volin += lay->volbelow;
            //if(l!=0) off1[2] += volin/area;
            if(l!=0) off1[2] += lay->volbelow/area;
            // fill mgeom
            k=0;
            conf->na[l]=0;
            while (1) {
                for (i=0; i<ml->sc[0]; i++) {
                    for (j=0; j<ml->sc[1]; j++) {
                        off2[0] = i*lay->tr1[0] + j*lay->tr2[0] + shift[l][0];
                        off2[1] = i*lay->tr1[1] + j*lay->tr2[1] + shift[l][1];
                        for (m=0; m<lay->nb; m++) {
                            conf->laysite[n] = &lay->site[m];
                            conf->coord[n][0] = lay->site[m].coord[0]+off1[0]+off2[0];
                            conf->coord[n][1] = lay->site[m].coord[1]+off1[1]+off2[1];
                            conf->coord[n][2] = lay->site[m].coord[2]+off1[2];
                            conf->rms[n] = rms;
                            mgeom->site[n].pidx = p;
                            mgeom->site[n].vmask = 1.0;
                            //conf->sidx[n] = sidx;
                            n++;
                            conf->na[l]++;
                        }
                    }
                }
                off1[0] += lay->trperp[0];
                off1[1] += lay->trperp[1];
                if(k>=nlay[l]-1) break;
                off1[2] += lay->trperp[2];
                k++;
                p++;
            }
            LOGMSG(2,"done\n");
            if(l>=ml->ns-1) break;
            //volin = lay->vol - lay->volbelow;
            off1[2] += (lay->vol - lay->volbelow)/area;
            conf->zint[l] = off1[2];
            l++;
        }
        off1[2] += lay->trperp[2];
        for (i=0; i<3; i++) {
            mgeom->trans[i] = off1[i];
        }
        
        // left lead:
        for (i=0; i<3; i++) {
            rg->trans[i] = lay->trperp[i];
        }
        rg->na = lay->nb;
        rg->tp = 1;
        rg->nmtr = lay->nmtr;
        rg->sc[0]= ml->sc[0];
        rg->sc[1]= ml->sc[1];
        rg->cr = mgeom->cr;
        rg->dimcf = 1.0;
        rg->scx[0]=rg->scx[1]=rg->scx[2]=1.0;
        rg->trv[0][0] = vec1[0];
        rg->trv[0][1] = vec1[1];
        rg->trv[1][0] = vec2[0];
        rg->trv[1][1] = vec2[1];
        rg->site = (t_trsite *)malloc(sizeof(t_trsite)*rg->na);
        for (i=0; i<rg->na; i++) {
            if (lay->site[i].nc > 1) EXITMSG("Do not take alloys as leads.");
            for (j=0; j<3; j++) {
                rg->site[i].coord[j] = lay->site[i].coord[j]+off1[j];
            }
            rg->site[i].at = lay->site[i].at[0];
        }
        
        // right lead:
        lay = &ml->slab[0].layer;
        for (i=0; i<3; i++) {
            lg->trans[i] = lay->trperp[i];
        }
        lg->na = lay->nb;
        lg->tp = 1;
        lg->nmtr = lay->nmtr;
        lg->sc[0]= ml->sc[0];
        lg->sc[1]= ml->sc[1];
        lg->cr = mgeom->cr;
        lg->dimcf = 1.0;
        lg->scx[0]=lg->scx[1]=lg->scx[2]=1.0;
        lg->trv[0][0] = vec1[0];
        lg->trv[0][1] = vec1[1];
        lg->trv[1][0] = vec2[0];
        lg->trv[1][1] = vec2[1];
        lg->site = (t_trsite *)malloc(sizeof(t_trsite)*lg->na);
        for (i=0; i<lg->na; i++) {
            if (lay->site[i].nc > 1) EXITMSG("Do not take alloys as leads.");
            for (j=0; j<3; j++) {
                lg->site[i].coord[j] = lay->site[i].coord[j]-lg->trans[j];
            }
            lg->site[i].at = lay->site[i].at[0];
        }
        
        lg->rcoord=(double (*)[3])malloc(sizeof(double[3])*lg->na);
        rg->rcoord=(double (*)[3])malloc(sizeof(double[3])*rg->na);
        mgeom->rcoord=(double (*)[3])malloc(sizeof(double[3])*mgeom->na);
        extrafields(lg);
        extrafields(rg);
    }
    
    void fill_conf(t_config *conf, t_atomset *atoms) {
        int l, n=0, i, j;
        t_trsite *site;
        t_basket *basket;
        double c;
        
        for(l=0; l<conf->ns; l++) {
            basket = conf->basket[l];
            basket->nfix = conf->nfix[l];
            reset_basket(basket);
            
            for(i=0; i<conf->na[l]; i++, n++) {
                site = &conf->mg.site[n];
                // determine disorder:
                site->at = pick_atom(conf->laysite[n], basket, atoms);
                c = conf->rms[n]/sqrt(site->at->mass);
                for (j=0; j<3; j++) {
                    site->coord[j] = conf->coord[n][j] + c*randn();
                }
            }
        }
        extrafields(&conf->mg);
    }
    
    void makedir(char *path, int dry) {
        if (!dry) {
            struct stat st = {0};
            if(stat(path, &st) == -1) {
                if (mkdir(path, 0777) != 0) {
                    EXITMSG("Failed to create directory %s", path);
                }
            }
        }
    }
    
    void copy_atoms(t_atomset *atoms, int dry) {
        int i, j, flag[atoms->num];
        t_atom *ptr, *p2;
        
        // checking for identical labels: (still has to be tested)
        for(i=0, ptr=atoms->ptr; i<atoms->num; i++, ptr=ptr->next) {
            flag[i]=0;
            for(j=0, p2=atoms->ptr; j<i; j++, p2=p2->next) {
                if(strcmp(p2->label, ptr->label)==0) {
                    flag[i]++;
                }
            }
        }
        for(i=0, ptr=atoms->ptr; i<atoms->num; i++, ptr=ptr->next) {
            if(flag[i]) {
                char str[MAXSTRLEN+1], *l=ptr->label;
                sprintf(str,"_%d",flag[i]);
                ptr->label = (char *)realloc(l, sizeof(char)*(strlen(str)+strlen(l)+1));
                strcat(ptr->label, str);
            }
        }
        // the above still has to be tested !
        
        makedir("./atoms0", dry);
        for(i=0, ptr=atoms->ptr; i<atoms->num; i++, ptr=ptr->next) {
            char buf[100+strlen(ptr->path)];
            sprintf(buf,"cp -f %s atoms0/%s", ptr->path, ptr->label);
            LOGMSG(3,"Copied %s to atoms0/%s",ptr->path, ptr->label);
            system(buf);
        }
        LOGMSG(2,"writing atomlist");
        if(!dry) writeatomlist(atoms);
        LOGMSG(2,"done with atoms");
    }
    
    void free_conf(t_config *conf) {
        free(conf->na);
        free(conf->basket);
        free(conf->nfix);
        free(conf->laysite);
        free(conf->coord);
        free(conf->mg.site);
        free(conf->lg.site);
        free(conf->rg.site);
        free(conf->mg.rcoord);
        free(conf->lg.rcoord);
        free(conf->rg.rcoord);
        free(conf->diskind);
        free(conf->rms);
        //free(conf->sidx);
        free(conf->zint);
        conf->ns = 0;
    }
    
    void write_pdb(t_config *conf, char *fname) {
        int i, j;
        t_trgeom *mg=&conf->mg;
        char *lbl, name[5];
        double x, y, z;
        FILE *fp;
        
        fp = fopen(fname,"w");
        for (i=0; i<mg->na; i++) {
            x = mg->rcoord[i][0];
            y = mg->rcoord[i][1];
            z = mg->rcoord[i][2];
            lbl = mg->site[i].at->label;
            j=strcspn(lbl, "0123456789_");
            if (j > 4) j=4;
            strncpy(name, lbl, j);
            name[j]=0;
            fprintf(fp,"ATOM  %5i %-4.4s              %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s\n",
                i+1,name,x,y,z,1.0,0.0,name);
        }
        fprintf(fp,"END\n");
        fclose(fp);
    }
    
/*    writesidx(t_config *conf, char *fname) {
        int i;
        FILE *fp;
        
        fp = fopen(fname, "w");
        for(i=0; i<conf->mg.na; i++) {
            fprintf(fp," %d\n", conf->sidx[i]);
        }
        fclose(fp);
    }*/
    
    write_geom_i(t_config *conf, char *fname) {
        int i;
        FILE *fp;

        fp = fopen(fname, "w");
        fprintf(fp,"----- POSITION OF INTERFACES -----\n");
        fprintf(fp," %d\n", conf->nint);
        for(i=0; i<conf->nint; i++) {
            fprintf(fp," %#15.11f\n", conf->zint[i]);
        }
        fclose(fp);
    }
    
    int make_link(char *target, char* lname) {
        if(test_file(lname)) {
            remove(lname);
        }
        if(symlink(target, lname) == 0) {
            return(0);
        } else {
            return(1);
        }
    }
    
    void setpars_dostuff(t_multilay *ml, t_par *par, char *path, int level, t_xopt xopt) {
        int i;
        t_valptr *ptr;
        char buf[strlen(path)+2*MAXSTRLEN];
        char buf2[MAXSTRLEN];
        FILE *fp;

        if (par!=NULL) {
            sprintf(buf,"%s/dirlist",path);
            fp = fopen(buf,"w");
            fprintf(fp," %d\n",ml->npar-level);
            fprintf(fp," %d\n",par->nval);
            LOGMSG(3, "Looping over values of %c ...", par->symbol);
            for (i=0; i<par->nval; i++) {
                for (ptr=par->ptr; ptr!=NULL; ptr=ptr->next) {
                    *(ptr->val) = par->val[i];
                }
                sprintf(buf2,"%c%g", par->symbol, (par->type-1?par->val[i].d:par->val[i].i));
                fprintf(fp," %s\n", buf2);
                sprintf(buf,"%s/%s", path, buf2);
                makedir(buf, xopt.dry);
                setpars_dostuff(ml, par->next, buf, level+1, xopt);
            }
            fclose(fp);
        } else {
            t_config conf;
            int icf, j;
            char *pos, lname[strlen(path)+2*MAXSTRLEN];
            
            printf("\ndoing stuff in folder %s\n", path);
            sprintf(buf,"%s/dirlist",path);
            fp = fopen(buf,"w");
            fprintf(fp," 0\n");
            fprintf(fp," %d\n", ml->nconf);
            
            prep_conf(&conf, ml);
            
            for (icf=0; icf<ml->nconf; icf++) {
                fprintf(fp," cf-%d\n", icf+1);
                LOGMSG(2,"Making configuration %d", icf+1);
                sprintf(buf,"%s/cf-%d",path,icf+1);
                makedir(buf, xopt.dry);
                
                //make the actual configuration:
                fill_conf(&conf, &ml->atoms);
                
                // output stuff:
                if (!xopt.dry){
                    //link atoms folder:
                    buf[0] ='\0';
                    for(j=0; j<level+1; j++) {
                        strcat(buf, "../");
                    }
                    pos = &buf[strlen(buf)];
                    
                    sprintf(pos, "atoms0");
                    sprintf(lname, "%s/cf-%d/atoms", path, icf+1);
                    if (make_link(buf, lname)) EXITMSG("Could not make link %s", lname);
                    //link atomlist:
                    sprintf(pos,"atomlist");
                    sprintf(lname,"%s/cf-%d/atomlist", path, icf+1);
                    if (make_link(buf, lname)) EXITMSG("Could not make link %s", lname);
                    //write geom files:
                    sprintf(buf,"%s/cf-%d/geom_m",path,icf+1);
                    writegeom(&conf.mg, buf);
                    sprintf(buf,"%s/cf-%d/geom_l",path,icf+1);
                    writegeom(&conf.lg, buf);
                    sprintf(buf,"%s/cf-%d/geom_r",path,icf+1);
                    writegeom(&conf.rg, buf);
                    //write slabidx file:
                    //sprintf(buf,"%s/cf-%d/slabidx",path,icf+1);
                    //writesidx(&conf, buf);
                    
                    //write interface positions:
                    sprintf(buf,"%s/cf-%d/geom_i",path,icf+1);
                    write_geom_i(&conf, buf);
                    //write pdb file:
                    if (xopt.vis!=0){
                        sprintf(buf,"%s/cf-%d/vis.pdb",path,icf+1);
                        write_pdb(&conf, buf);
                    }
                }
            }
            fclose(fp);
            free_conf(&conf);
        }
    }
    
    void write_dirloop() {
        char white[100];
        char str1[MAXSTRLEN], str[MAXSTRLEN];
        char path[400]=".", p[101]=".";
        FILE *fp, *fl;
        int lev, ndir, i, l, len=1;
        
        fl = fopen("./dirloop", "w");
        fchmod(fileno(fl), S_IRWXU | S_IRWXG | S_IRWXO);
        fprintf(fl,"#!/bin/sh \n\n");
        
        fp = fopen("./dirlist", "r");
        if(fscanf(fp," %d", &lev)!=1) EXITMSG("Error reading dirlist.");
        fclose(fp);
        
        for (l=0; l<=lev; l++) {
            fprintf(fl,"%sfor d%d in",white,l);
            fp = cat_openread("dirlist", path);
            if(fscanf(fp," %*d %d", &ndir)!=1) EXITMSG("Error reading dirlist.");
            
            for (i=0; i<ndir; i++) {
                fscanf(fp, " %s", str);
                fprintf(fl," %s",str);
            }
            fclose(fp);
            fprintf(fl,"\n%sdo\n",white);
            strcat(white,"    ");
            strcat(path, "/");
            strcat(path, str);
            len += sprintf(p+len, "/$d%d", l);
        }
        fprintf(fl,"%s#cp ./trans.conf %s/trans.conf\n", white, p);
        fprintf(fl,"%s#cp ./runscript %s/runscript\n", white, p);
        fprintf(fl,"%s#cd %s\n", white, p);
        fprintf(fl,"%s#sbatch runscript\n", white);
        fprintf(fl,"%s#cd ", white);
        for (l=0; l<=lev; l++) {
            fprintf(fl,"../");
        }
        fprintf(fl,"\n");
        for (l=lev; l>=0; l--) {
            white[l*4]='\0';
            fprintf(fl,"%sdone\n",white);
        }
        fclose(fl);
    }
    
    int main(int argc, char *argv[])
    {
        t_multilay ml = {.atoms = {.ptr = NULL, .num=0}};
        t_xopt xopt;
        t_db_entry *dbe=NULL;
        
        LOGMSG(1,"Geometry generator using database\n");
        
        initrnd();

        readopts(argc, argv, &xopt);
        
        get_db_dir(argv, &xopt);

        read_input_ml(&ml);
        
        read_db(&dbe, &ml, xopt.dbloc);
        
        read_vars(&ml);
        
        copy_atoms(&ml.atoms, xopt.dry);
        
        setpars_dostuff(&ml, ml.parptr, ".", 0, xopt);
        
        write_dirloop();
        
        printf("Finished!\n");
    }
