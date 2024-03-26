#!/usr/bin/env qoctave
cfr=2.581280756014811e-06*1e9;
opt.Robust='on';
opt.confint=0.7;
opts.MaxIter=2e4;
gfact=2.1;
rncut=0;
# To use effective magnetic moment, set mgi=2
mgi=2;
rwp=0;


try
    lbase=load('nlofs');
catch
    lbase=0;
end;
disp(['NLadd=' num2str(lbase)]);
# Read parameters from geometry
tic
try
    fls=ls('-d','L*/cf*/transout.h5');
    ifile=strtrim(fls(round(end/2),:));
%    fls=ls('-1 --sort=time','L*/cf*/transout.h5');
%    ifile=strtrim(fls(1,:));
    a0=h5read(ifile,'/geom/lgeo/anm');
    base=h5read(ifile,'/geom/lgeo/base');
    trperp=h5read(ifile,'/geom/lgeo/perp_trans');
    sc=h5read(ifile,'/geom/lgeo/sc_size');
    ar=det(base);
    rcfa0=a0^2*ar*cfr/100;
    lscl=a0*trperp(3);
    scsz=sc(1)*sc(2);
    rcfa=scsz*rcfa0;
    cf=1/pi/4*gfact;  
    lsh=h5read(ifile,'/cond/l_shc')/rcfa;
    rsh=h5read(ifile,'/cond/r_shc')/rcfa;
    shres=.5*(1./sum(lsh)+1./sum(rsh));
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
damps=[];
mdamps=[];
mvals=[];
momst=[];
for i2=1:length(dirs2)
    ldir=dirs2(i2).name;
    dirs3=dir([ldir,'/cf*']);
    cnd=[];
    res0=[];
    dmp=[];
    mgm=[];
    mgm1=[];
    mgv0=[];
    moms=[];
    for i3=1:length(dirs3)
        ifile=[ldir,'/',dirs3(i3).name,'/transout.h5'];
        try
            flags=h5read(ifile,'/compat');
            rdy=1;
        catch
            rdy=0;
        end
        if (rdy==1) 
            if (flags(5)==1), cnd=[cnd;h5read(ifile,'/cond/totlr')];end;
            if (flags(6)==1), cnd=[cnd;h5read(ifile,'/cond/totrl')];end;
            mgv=h5read(ifile,'/damp/magvals');
            moms=[moms;[mgv(1),mgv(2)]];
            if (flags(3)==1), 
                dmp=[dmp;diag(real(h5read(ifile,'/damp/total')))];
                mgm=[mgm;[mgv(1);mgv(1)]];
                mgm1=[mgm1;[mgv(2);mgv(2)]];
                mgv0=[mgv0;mgv(3)];
            end
        end
    end;
    l0=(str2num(ldir(2:end))+lbase)*lscl;
    if (size(cnd,1)>0)
        cnd=cnd/rcfa;
        res0=[1./cnd-shres];
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
        momst=[momst;[l0,mean(moms,1)/scsz/str2num(ldir(2:end))]];
    end
    if (size(dmp,1)>0)
        dmp=dmp*cf;
        ml0=repmat(l0,size(dmp,1),1);
        m0=mean(mgm);
        m1=mean(mgm1);
        mdmp=mean(dmp,1);
        ddmp=std(dmp,0,1);
        mvals=[mvals;[l0,mean(mgv0)]];
        mdamps=[mdamps;[m0,m1,l0,mdmp,ddmp]];
        damps=[damps;[mgm,mgm1,ml0,dmp]];
    end;
    
end;

if (size(mconds,1)>1) 
    mconds=sortrows(mconds,1);
    conds=sortrows(conds,1);
    res=sortrows(res,1);
    mres=sortrows(mres,1);
    momst=sortrows(momst,1);
    moms=mean(momst(:,2:3),1);
end;        

disp('----------');
if (nnz(mres(:,1)>rncut)>1)

    pp=nnz(mres(:,1)>=rncut);
    ind=mres(:,1)>=rncut;

    q1=polyfit(mres(ind,1),mres(ind,2),1);
    pin=[q1(1)   q1(2)   .5    -.5];
    
    
    wt=1./sqrt(mres(ind,1));
#disp('----------1');
    [cfs2,dev2,pout2]=resfit(mres(ind,1),mres(ind,2),pin,wt,opt);
    disp('res1');
    pin=cfs2;
#    cfs2
#    pin(3:4)=pin(3:4)*.9;
    ind1=res(:,1)>=rncut;
#disp('----------2');
    [cfs1,dev1,pout1]=resfit(res(ind1,1),res(ind1,2),pin,res(ind1,1).^rwp,opt);
    disp('res2');

    rhos=[cfs1(1),dev1(1);cfs2(1),dev2(1)]*1e2;
    Rs=[cfs1(2),dev1(2);cfs2(2),dev2(2)];

    disp(rhos);
    
    [rc1,rd1]=bisqpolyfit(mres(ind,1),mres(ind,2),1,opt);
    [rc2,rd2]=bisqpolyfit(res(ind1,1),res(ind1,2),1,opt);

    rhosL=[rc1(1),rd1(1);rc2(1),rd2(1)]*1e2;
    RsL=[rc1(2),rd1(2);rc2(2),rd2(2)];

    disp(rhosL);

else
    rhos=zeros(2,2);
end
mdamps
if (size(mdamps,1)>1) 
    damps=sortrows(damps,1);
    mdamps=sortrows(mdamps,1);
    [d0,dd]=bisqpolyfit(damps(:,1),damps(:,4),1,opt);
    [d1,dd1]=bisqpolyfit(mdamps(:,1),mdamps(:,4),1,opt);
    disp('damp1');
    [d0q,ddq]=bisqpolyfit(damps(:,2),damps(:,4),1,opt);
    [d1q,dd1q]=bisqpolyfit(mdamps(:,2),mdamps(:,4),1,opt);
    disp('damp2');
    
    alph=[d0(1),dd(1);d1(1),dd1(1)]*1e3;
    dEnch=[d0(2),dd(2);d1(2),dd1(2)]*1e3*mean(mdamps(:,3)./mdamps(:,1));
    alph1=[d0q(1),ddq(1);d1q(1),dd1q(1)]*1e3;
    dEnch1=[d0q(2),ddq(2);d1q(2),dd1q(2)]*1e3*mean(mdamps(:,3)./mdamps(:,2));
    disp(alph);
    
    mfact=mean(mvals(:,2));
end;


save -v7 result.mat mconds conds res mres rhos damps mdamps alph a0 sc scsz opt ar lscl gfact rcfa lsh rsh shres Rs dEnch alph1 dEnch1 mvals mfact momst moms
#disp(mres);


fid=fopen('base.txt','w');
if (nnz(mres(:,1)>rncut)>1)
fprintf(fid,'-------------Resitivity (rho) and interface term (Rb)--------------\n');
fprintf(fid,'                             rho       d(rho)      Rb        d(Rb)\n');
fprintf(fid,'nonlin(all set)       : %10.5f %10.5f %10.5f %10.5f\n',rhos(1,:),Rs(1,:));
fprintf(fid,'nonlin(mean vals)     : %10.5f %10.5f %10.5f %10.5f\n',rhos(2,:),Rs(2,:));
fprintf(fid,'   lin(all set)       : %10.5f %10.5f %10.5f %10.5f\n',rhosL(1,:),RsL(1,:));
fprintf(fid,'   lin(mean vals)     : %10.5f %10.5f %10.5f %10.5f\n',rhosL(2,:),RsL(2,:));
fprintf(fid,'-------------------------------------------------------------------\n');
fprintf(fid,'\n');
end

if (size(mdamps,1)>1)
fprintf(fid,'---------------------------Gilbert damping-------------------------\n');
fprintf(fid,'                            alpha     d(alpha)   G_enh     d(G_enh)\n');
fprintf(fid,'  all dataset (Meff)  : %10.5f %10.5f %10.5f %10.5f\n',alph1(1,:),dEnch1(1,:));
fprintf(fid,'    mean.vals.(Meff)  : %10.5f %10.5f %10.5f %10.5f\n',alph1(2,:),dEnch1(2,:));
fprintf(fid,'  all dataset (Msat)  : %10.5f %10.5f %10.5f %10.5f\n',alph(1,:),dEnch(1,:));
fprintf(fid,'    mean.vals.(Msat)  : %10.5f %10.5f %10.5f %10.5f\n',alph(2,:),dEnch(2,:));
fprintf(fid,'-------------------------------------------------------------------\n');
fprintf(fid,'\n');
end

fclose(fid);
