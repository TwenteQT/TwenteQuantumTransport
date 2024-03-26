#!/usr/bin/env qoctave 
rescut=20;
rscut=2;
rncut=2;
dfdcut=5;
needsdf=1;
rwp=.0;
printmeans=1;
global confi=0.7;
global verb=0;
global resiter=20000;
global sdliter=20000;



function [cfs,dev,cout]=resfit(x1,rr,pin,wt)
    global confi;
    global resiter;
    opts.MaxIter=resiter;
    opts.Tol=1e-8;
    opts.Bounds=[0,Inf;0,Inf;-Inf,0;0,Inf];
    opts.confint=confi;
    opts.Robust='on';
    if (exist('wt'))
      opts.Weights=wt;
    end;

    mf1=@(x,p) p(1)+p(2)*x+p(3)./(p(4)+x*p(2));    
    [cout,dev,opi]=bisqnonlin(mf1,x1,rr,pin,opts);

    a=cout(1);
    b=cout(2);
    c=cout(3);
    d=cout(4);

    rs1=[b,a+c/d,(-a + d)/sqrt(-4*c + (a - d)^2),(-2*c + d*(-a + d))/(sqrt(-4*c + (a - d)^2)*d)];
    rs2=[b,a+c/d,(a - d)/sqrt(-4*c + (a - d)^2),(2*c + (a - d)*d)/(sqrt(-4*c + (a - d)^2)*d)];

    betadev=2*abs((-a + d)/(-4*c + (a - d)^2)^1.5)*abs(dev(3))+4*abs(c/(-4*c + (a - d)^2)^1.5)*(abs(dev(1))+abs(dev(4)));
    gammadev=2*abs((c*(a + d))/((-4*c + (a - d)^2)^1.5*d))*abs(dev(1))+2*abs((2*c + a*(-a + d))/((-4*c + (a - d)^2)^1.5*d))*abs(dev(3))+2*abs((c*(a^2 - 4*c - 3*a*d))/((-4*c + (a - d)^2)^1.5*d^2))*abs(dev(4));
    rbdev=abs(dev(1))+1/abs(d)*abs(dev(3))+abs(c/d^2)*abs(dev(4));

    cfs=rs2;
    dev=[dev(2),rbdev,betadev,gammadev];
    
endfunction;


% Read parameters
isoct=length(ver('Octave'));
pf=fopen('params.txt');
v=fread(pf,Inf,'uchar=>uchar')';
fclose(pf);
pars=strsplit(strrep(v,'#','%'),char(10));
for ii=1:length(pars)
    if strfind(pars{ii},'=')
        eval([pars{ii},';']);
    end
end;
format long;

cf=1/pi/4*gfact;
ar1=sc*a0^2*ar;
lw=a0*d;
ladd=bmom/lch*lw;
lar=round(bmom/lch);
ladd=lar*lw;
nladd=bmom/lch;
cfr=2.581280756014811e-06*1e9;

ldr=dir('L*');
cond=[];
damp=[];
sharv=[];
condbs=[];
res=[];
resbs=[];
diffpol=[];
resbase=[];
dampbase=[];
diffbase=[];
meanval=[];
meandev=[];


for id=1:length(ldr)
    if (ldr(id).isdir) 
        dn=ldr(id).name;
        nl0=sscanf(dn,'L%d');
        x0=nl0*lw;
        l0=x0+ladd;
        nl0t=nl0+nladd;
        m0=nl0*lch+bmom;
        xco=[l0,nl0t,x0,nl0,m0];
        rcfa=ar1*cfr/100;
        dcfa=cf/sc/lch*lw;
        mus=lch/lw;
        ldr1=dir([dn,'/cf*']);
        dampt=[];
        condbst=[];
        sharvt=[];
        condt=[]; 
        rest=[];
        resbst=[];
        hasdata=0;
        for id=1:length(ldr1)
            if (ldr1(id).isdir) 
                cdn=[dn,'/',ldr1(id).name];
                andf=[cdn,'/andout'];
                af=fopen(andf);
                if (af~=-1)
                    ao=fread(af,Inf,'uchar=>uchar')';
                    fclose(af);
                    td=length(strfind(ao,'carefull with decomposition'));
                    td=1;
                    olines=deblank(strsplit(ao,char(10)));
                    for ind=1:length(olines)
                        if (strfind(olines{ind},'Damping constants, Re part'))
                            tp=[str2num(olines{ind+1}(10:end)),str2num(olines{ind+2}(10:end))];
                            dampt=[dampt;tp];
                            break;
                        end;
                    end;
                    sh=[];
                    for ind=1:length(olines)
                        if (strfind(olines{ind},'Sharvin'))
                            if (strfind(olines{ind+1},'spin'))
                                sh=[sh;str2num(olines{ind+1}(strfind(olines{ind+1},':')+1:end))];
                            else
                                sh=[sh;str2num(olines{ind+1}(1:end))];
                            end;
                            if (size(sh,1)==2); 
                                break;
                            end;
                        end;
                    end;
                    tl=sh(1,1)+sh(1,2);
                    tr=sh(2,1)+sh(2,2);
                    sharvt=[sharvt;[sh(1,:),sh(2,:),tl,tr]];
                    for ind=1:length(olines)
                        if (strfind(olines{ind},'Grand total'))
                            tp=str2num(olines{ind}(strfind(olines{ind},':')+1:end));
                            condt=[condt;tp];
                            rest=[rest;[1./tp,1./tp-.5*(1/tl+1/tr)]];
                            break;
                        end;
                    end;
                    % tp=zeros(1,4);
                    % tpi=0;
                    if (td~=0) 
                        for ind=1:length(olines)
                            if (strfind(olines{ind},'L->R') ||strfind(olines{ind},'R->L'))
                                tp=[str2num(olines{ind+2}(strfind(olines{ind+2},':')+1:end)),str2num(olines{ind+3}(strfind(olines{ind+3},':')+1:end))];
                                % tpi=tpi+1;
                                condbst=[condbst;tp];
                                resbst=[resbst;[1./tp]];                    
                            end;
                        end;
                    end;
%                    if (nl0 == 50) || (nl0==40), disp([nl0,tp]);end; 
                    % tp=tp/tpi;
                    % condbst=[condbst;tp];
                    % resbst=[resbst;[1./tp]];                    
                end;
            end;
        end;
        if (size(rest,1)>0)
            condt=condt/rcfa;
            rest=rest*rcfa;
            dampt=dampt*dcfa;
            sharvt=sharvt/rcfa;
            res=[res;[repmat(xco,size(rest,1),1),rest]];
            cond=[cond;[repmat(xco,size(condt,1),1),condt]];

            if (size(resbst,1)>0)
                resbst=resbst*rcfa;
                condbst=condbst/rcfa;
                resbs=[resbs;[repmat(xco,size(resbst,1),1),resbst]];
                condbs=[condbs;[repmat(xco,size(condbst,1),1),condbst]];
            end

            if (size(dampt,1)>0)
                damp=[damp;[repmat(xco,size(dampt,1),1),dampt]];
            end;
            sharv=[sharv;[repmat(xco,size(sharvt,1),1),sharvt]];

            if (size(condt,1)>1); conddev=std(condt,[],1);
            else conddev=zeros(size(condt));
            end;

            if (size(resbst,1)>0)
                if (size(condbst,1)>1); condbstdev=std(condbst,[],1);
                else condbstdev=zeros(size(condbst));
                end;
                rbstdev=condbstdev./mean(condbst,1).^2;
            end 

            if (size(sharvt,1)>1); sharvdev=std(sharvt,[],1);
            else sharvdev=zeros(size(sharvt));
            end;
            rdev=conddev./mean(condt,1).^2;

            resbase=[resbase;[l0,1./mean(condt,1),rdev]];

            if (isempty(dampt)); 
                dampt=[nan,nan,nan,nan];
            else
                dampbase=[dampbase;l0,mean([dampt(:,1);dampt(:,4)],1),std([dampt(:,1);dampt(:,4)],[],1)];
            end;
            if (size(dampt,1)>1); dampdev=std(dampt,[],1);
            else dampdev=zeros(size(dampt));
            end;

            if (size(resbst,1)>0)
                gm=.5*(condbst(:,2)+condbst(:,3));
                gu=condbst(:,1);
                gd=condbst(:,4);

                df=[gu,gd,gm];
                format long;
                %disp([sh(1,1),sh(1,2)]/rcfa);
                diffpol=[diffpol;[l0*ones(size(df,1),1),df]];
                diffbase=[diffbase;[l0,mean(df,1),std(df,[],1)]];
                mcb=mean(condbst,1);
                mrb=mean(resbst,1);
            else
                mcb=nan*ones(1,4);
                mrb=nan*ones(1,4);
                condbstdev=nan*ones(1,4);
                rbstdev=nan*ones(1,4);
            end
            % disp(condbst)   1-5p    6p            7p                           8p                                 9-12 13-16 17-22         23-26                    
            meanval=[meanval;[xco,mean(condt,1),1./mean(condt,1),1./mean(condt,1)-.5*(sum(1./mean(sharvt(:,5:6),1))),mcb,mrb,mean(sharvt,1),mean(dampt,1)]];
            meandev=[meandev;[xco,conddev,rdev,rdev,condbstdev,rbstdev,sharvdev,dampdev]];
        end;
    end;
end;
meanval=sortrows(meanval,1);
meandev=sortrows(meandev,1);
cond=sortrows(cond,1);
res=sortrows(res,1);
if (size(resbs,1)>0)
    resbs=sortrows(resbs,1);
    condbs=sortrows(condbs,1);
end;
if (size(damp,1)>0)
    damp=sortrows(damp,1);
    dampbase=sortrows(dampbase,1);
end

if (size(diffbase,1)>0)
    diffbase=sortrows(diffbase,1);
    diffpol=sortrows(diffpol,1);
end;

sharv=sortrows(sharv,1);
resbase=sortrows(resbase,1);
save -v7 results.mat cond damp sharv condbs res resbs resbase dampbase meanval meandev diffpol
save -ascii dampout.txt  dampbase
save -ascii resout.txt  resbase
save -ascii diffout.txt  diffbase
save -ascii diffall.txt diffpol
save -ascii dampall.txt  damp
save -ascii resall.txt  res


opt.Robust='on';
opt.confint=confi;
%%%%%%%  find resistivity
if (nnz(resbase(:,1)>=rescut)>1)
    ind1=res(:,1)>=rescut;
    [rho,rhodev]=bisqpolyfit(res(ind1,1),res(ind1,7),1,opt);
    rr=[rho,rhodev];
    

    ind=resbase(:,1)>=rescut;
    pp=nnz(resbase(:,1)>=rescut);
    [rho,S]=wpolyfit(resbase(ind,1),resbase(ind,2)-rcfa*.5*(1/tl+1/tr),resbase(ind,3),1);
    rhodev=(sqrt(sumsq(inv(S.R'))'/S.df)*S.normr)'*tinv(1-(1-confi)/2,length(ind)-2);
    rr=[rr,[rho',rhodev']];
    rr(1,:)=rr(1,:)*1e2;
else
    rr=zeros(2,4);
end


%%%%%%%%%% Find spin-resistivities (assuming independent), be aware it is wrong with SOC
if (nnz(diffbase(:,1)>rscut)>1)

    ind1=diffpol(:,1)>rscut;    

    [rho,rhodev]=bisqpolyfit(diffpol(ind1,1),1./(diffpol(ind1,2)+diffpol(ind1,4))-rcfa*.5*(1/sh(1,1)+1/sh(1,2)),1,opt);
    rho(1,1)=rho(1,1)*1e2;rhodev(1,1)=rhodev(1,1)*1e2;
    ru=[rho,rhodev];

    [rho,rhodev]=bisqpolyfit(diffpol(ind1,1),1./(diffpol(ind1,3)+diffpol(ind1,4))-rcfa*.5*(1/sh(1,1)+1/sh(1,2)),1,opt);
    rho(1,1)=rho(1,1)*1e2;rhodev(1,1)=rhodev(1,1)*1e2;
    rd=[rho,rhodev];


    ind=diffbase(:,1)>rscut;    
    pp=length(ind);

    dat=1./(diffbase(:,2)+diffbase(:,4))-rcfa*.5*(1/sh(1,1)+1/sh(1,2));
    dev=(diffbase(:,5)+diffbase(:,7)).*dat.^2;

    [rho,S]=wpolyfit(diffbase(ind,1),dat(ind),dev(ind),1);
    rhodev=(sqrt(sumsq(inv(S.R'))'/S.df)*S.normr)'*tinv(1-(1-confi)/2,pp-2);
    rho(1,1)=rho(1,1)*1e2;rhodev(1,1)=rhodev(1,1)*1e2;
    ru=[ru,[rho',rhodev']];

    dat=1./(diffbase(:,3)+diffbase(:,4))-rcfa*.5*(1/sh(2,1)+1/sh(2,2));
    dev=(diffbase(:,6)+diffbase(:,7)).*dat.^2;
    
    [rho,S]=wpolyfit(diffbase(ind,1),dat(ind),dev(ind),1);
    rhodev=(sqrt(sumsq(inv(S.R'))'/S.df)*S.normr)'*tinv(1-(1-confi)/2,length(diffbase(ind,1))-2);
    rho(1,1)=rho(1,1)*1e2;rhodev(1,1)=rhodev(1,1)*1e2;
    rd=[rd,[rho',rhodev']];

    rt=zeros(2,4);
    rt(:,[1,3])=1./(1./ru(:,[1,3])+1./rd(:,[1,3]));
    rt(:,[2,4])=1./ru(:,[1,3]).^2.*ru(:,[2,4])+1./rd(:,[1,3]).^2.*rd(:,[2,4]);
    rt(:,[2,4])=rt(:,[1,3]).^2.*rt(:,[2,4]);
else
    ru=zeros(2,4);
    rd=zeros(2,4);
    rt=zeros(2,4);
end;




%resb=[diffbase(:,1),1./(diffbase(:,3)+diffbase(:,2))];
%resp=[diffpol(:,1),1./(diffpol(:,3)+diffpol(:,2))];


if (nnz(diffbase(:,1)>rncut)>1)

pp=nnz(resbase(:,1)>=rncut);
ind=resbase(:,1)>=rncut;

pin=[1.5313    0.0053   -0.0299    0.0272];
yy=resbase(ind,2)-rcfa*.5*(1/tl+1/tr);
wt=1./sqrt(resbase(ind,3));
[cfs2,dev2,pout2]=resfit(resbase(ind,1),yy,pin,wt);

pin=pout2;
ind1=res(:,1)>=rncut;
[cfs1,dev1,pout1]=resfit(res(ind1,1),res(ind1,7),pin,res(ind1,1).^rwp);

ra=[[cfs1(1),dev1(1),cfs2(1),dev2(1)]*1e2;cfs1(2),dev1(2),cfs2(2),dev2(2)];
else
  ra=zeros(2,4);
end


if (size(dampbase,1)>0)

    %%%%%%%  find damping constant
    [dc,dcdev]=bisqpolyfit([damp(:,1);damp(:,1)],[damp(:,6);damp(:,9)],1,opt);
    dp=[dc,dcdev]*1e3;


    [dc,dcdev]=bisqpolyfit(damp(:,1),damp(:,6),1,opt);
    dpxx=[dc,dcdev]*1e3;

    [dc,dcdev]=bisqpolyfit(damp(:,1),damp(:,9),1,opt);
    dpyy=[dc,dcdev]*1e3;


    [dc,S]=wpolyfit(dampbase(:,1),dampbase(:,2),dampbase(:,3),1);
    dcdev=(sqrt(sumsq(inv(S.R'))'/S.df)*S.normr)'*tinv(1-(1-confi)/2,size(dampbase,1)-2);
    dp=[dp,[dc',dcdev']*1e3];
end

if (needsdf==1), 

    if (size(diffbase,1)>0)&&0
        pp=size(diffbase,1);

        %%%%%%%  find spin-diffusion length
        ind=find(diffbase(:,1)>=dfdcut);
        [cfs2,dev2,mfu2,mfd2]=ndiffit2(diffbase(:,1),diffbase(:,2),diffbase(:,3),diffbase(:,4),pp,.5);
        ind=find(diffpol(:,1)>=dfdcut);
        [cfs1,dev1,mfu1,mfd1]=ndiffit2(diffpol(:,1),diffpol(:,2),diffpol(:,3),diffpol(:,4),pp,.5);

        for ii=1:3
            dfl(ii,1:4)=[cfs1(ii),dev1(ii),cfs2(ii),dev2(ii)];
        end;
    else
	dfl=zeros(2,4);
    end
end;

if (verb==1)
    disp('rho:');
    disp(rr);
    disp('ru:');
    disp(ru);
    disp('rd:');
    disp(rd);
    disp('rt:');
    disp(rt);
    disp('ra:');
    disp(ra);
    disp('damp:')
%    disp(dp)
    disp('sdf:')
    disp(dfl)
end;
%dp

fid=fopen('baseresults.txt','w');
fprintf(fid,'-------------Resitivity (rho) and interface term (Rb)--------------\n');
fprintf(fid,'                             rho       d(rho)      Rb        d(Rb)\n');
fprintf(fid,'FULL: nonlin total    : %10.5f %10.5f %10.5f %10.5f\n',ra(1,1:2),ra(2,1:2));
fprintf(fid,'FULL: lin (L>%6.2f)  : %10.5f %10.5f %10.5f %10.5f\n',rescut,rr(1,1:2),rr(2,1:2));
fprintf(fid,'FULL: parallel resist.: %10.5f %10.5f %10.5f %10.5f\n',rt(1,1:2),rt(2,1:2));
fprintf(fid,'FULL: up channel      : %10.5f %10.5f %10.5f %10.5f\n',ru(1,1:2),ru(2,1:2));
fprintf(fid,'FULL: down channel    : %10.5f %10.5f %10.5f %10.5f\n',rd(1,1:2),rd(2,1:2));
fprintf(fid,'-------------------------------------------------------------------\n');
if (printmeans==1)
fprintf(fid,'-------------Resistivity: fit of averaged mean values--------------\n');
fprintf(fid,'                             rho       d(rho)      Rb        d(Rb)\n');
fprintf(fid,'MEAN: nonlin total    : %10.5f %10.5f %10.5f %10.5f\n',ra(1,3:4),ra(2,3:4));
fprintf(fid,'MEAN: lin (L>%6.2f)  : %10.5f %10.5f %10.5f %10.5f\n',rescut,rr(1,3:4),rr(2,3:4));
fprintf(fid,'MEAN: parallel resist.: %10.5f %10.5f %10.5f %10.5f\n',rt(1,3:4),rt(2,3:4));
fprintf(fid,'MEAN: up channel      : %10.5f %10.5f %10.5f %10.5f\n',ru(1,3:4),ru(2,3:4));
fprintf(fid,'MEAN: down channel    : %10.5f %10.5f %10.5f %10.5f\n',rd(1,3:4),rd(2,3:4));
fprintf(fid,'-------------------------------------------------------------------\n');
fprintf(fid,'\n');
end;
fprintf(fid,'\n');


if (size(dampbase,1)>0)
fprintf(fid,'---------------------------Gilbert damping-------------------------\n');
fprintf(fid,'                            alpha     d(alpha)   G_enh     d(G_enh)\n');
fprintf(fid,'FULL: averaged        : %10.5f %10.5f %10.5f %10.5f\n',dp(1,1:2),dp(2,1:2));
fprintf(fid,'FULL: xx component    : %10.5f %10.5f %10.5f %10.5f\n',dpxx(1,1:2),dpxx(2,1:2));
fprintf(fid,'FULL: yy component    : %10.5f %10.5f %10.5f %10.5f\n',dpyy(1,1:2),dpyy(2,1:2));
fprintf(fid,'-------------------------------------------------------------------\n');
if (printmeans==1)
fprintf(fid,'-------------Damping: fit of averaged mean values------------------\n');
fprintf(fid,'MEAN: averaged       : %10.5f %10.5f %10.5f %10.5f\n',dp(1,3:4),dp(2,3:4));
fprintf(fid,'-------------------------------------------------------------------\n');
end;
fprintf(fid,'\n');
end

if (needsdf==1), 

    if (size(diffbase,1)>0)
%        fprintf(fid,'sdl: %10.5f %10.5f %10.5f %10.5f\n',dfl(1,:));
%        fprintf(fid,'beta: %10.5f %10.5f %10.5f %10.5f\n',dfl(2,:));
%        fprintf(fid,'offs : %10.5f %10.5f %10.5f %10.5f\n',dfl(3,:));
        %fprintf(fid,'sdl1: %10.5f %10.5f %10.5f %10.5f\n',dfl(4,:));
        %fprintf(fid,'sdl2: %10.5f %10.5f %10.5f %10.5f\n',dfl(5,:));
        %fprintf(fid,'bet1: %10.5f %10.5f %10.5f %10.5f\n',dfl(6,:));
        %fprintf(fid,'bet2: %10.5f %10.5f %10.5f %10.5f\n',dfl(7,:));
    end
end
fprintf(fid,'\n\n');
fprintf(fid,'Sharvin: %10.5f %10.5f \n',[sh(1,1),sh(1,2)]/rcfa);
fprintf(fid,'a0     : %10.5f\n',a0);
fprintf(fid,'mus    : %10.5f\n',mus);
fprintf(fid,'mus3   : %10.5f\n',mus/(a0^2*ar));

fclose(fid);
