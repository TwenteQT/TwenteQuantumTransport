#!/usr/bin/env qoctave
cfr=2.581280756014811e-06*1e9;
opt.Robust='off';
opt.confint=0.7;

# Read parameters from geometry
tic
try
    fls=ls('-d','L*/cf*'); 
    ifile=[strtrim(fls(1,:)),'/transout.h5'];
    a0=h5read(ifile,'/geom/lgeo/anm');
    base=h5read(ifile,'/geom/lgeo/base');
    trperp=h5read(ifile,'/geom/lgeo/perp_trans');
    sc=h5read(ifile,'/geom/lgeo/sc_size');
    ar=det(base);
    rcfa0=a0^2*ar*cfr/100;
    lscl=a0*trperp(3);
    scsz=sc(1)*sc(2);
    rcfa=scsz*rcfa0;
    #cf=1/pi/4*gfact;  
    #dcfa=cf/scz/lch*lscl;
    lsh=h5read(ifile,'/cond/l_shc')/rcfa;
    rsh=h5read(ifile,'/cond/r_shc')/rcfa;
    shres=.5*(1./lsh+1./rsh);
catch
    disp('no results yet');
    exit
end
toc

dirs2=dir('L*');
mconds=[];
conds=[];
mres=[];
res=[];
for i2=1:length(dirs2)
    ldir=dirs2(i2).name;
    dirs3=dir([ldir,'/cf*']);
    cnd=[];
    res0=[];
    for i3=1:length(dirs3)
        cfdir=[ldir,'/',dirs3(i3).name];
        if (exist([cfdir,'/transout.h5'],'file') && exist([cfdir,'/andout'])) 
            ifile=[cfdir,'/transout.h5'];
            cnd0=h5read(ifile,'/cond/lr');
            cnd=[cnd;[cnd0(1,1),cnd0(2,2)]];
        end
    end;
    if (size(cnd,2)>0)
        cnd=cnd/rcfa;
        l0=str2num(ldir(2:end))*lscl;
        res0=[1./cnd(:,1)-shres(1),1./cnd(:,2)-shres(2)];
        mcnd=mean(cnd,1);
        dcnd=std(cnd,0,1);
        mres0=mean(res0,1);
        dres0=std(res0,0,1);
        nconf=size(cnd,1);
        ml0=repmat(l0,nconf,1);
        mconds=[mconds;[l0,mcnd,dcnd]];
        conds=[conds;[ml0,cnd]];
        mres=[mres;[l0,mres0,dres0]];
        res=[res;[ml0,res0]];
    end
end;

if (size(mconds,1)>1) 
    mconds=sortrows(mconds,1);
    conds=sortrows(conds,1);
    res=sortrows(res,1);
    mres=sortrows(mres,1);
    [rfu,rdu]=bisqpolyfit(res(:,1),res(:,2),1,opt);
    [rfd,rdd]=bisqpolyfit(res(:,1),res(:,3),1,opt);
    rhos=[rfu(1),rfd(1),rdu(1),rdd(1)]*1e2;
end;        

save -v7 result.mat mconds conds res mres rhos
disp(rhos);
