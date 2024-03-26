function [cfs,dev,mfu,mfd]=diffit3(x1,y1,x2,y2,prec,pwr)

confi=.7;
maxit=200;
cfs1(2)=y1(end)*2-1;
cfs1(1)=5;
cfs1(3)=0;
cfs2(2)=-(y2(end)*2-1);
cfs2(1)=10;
cfs2(3)=0;
pp1=size(x1,1);
pp2=size(x2,1);

mf1=@(x,p) .5*(1+p(2))*(1+(1-p(2))/(1+p(2))*exp(-x/p(1)+p(3)));
[f,cfs1,kvg,iter,corp,covp,covr,stdresid,Z,r2]=leasqr(x1,y1,cfs1,mf1,1e-16,1000);
dev1=sqrt(diag(covp)*(pp1-3)/pp1)*tinv(1-(1-confi)/2,pp1-3); 

mf2=@(x,p) .5*(1-p(2))*(1+(1+p(2))/(1-p(2))*exp(-x/p(1)+p(3)));
[f,cfs2,kvg,iter,corp,covp,covr,stdresid,Z,r2]=leasqr(x2,y2,cfs2,mf2,1e-16,1e3);
dev2=sqrt(diag(covp)*(pp2-3)/pp2)*tinv(1-(1-confi)/2,pp2-3); 
%exit
ncfs1(1)=cfs1(1);
ncfs1(2)=cfs1(3);
ncfs2(1)=cfs2(1);
ncfs2(2)=cfs2(3);

beta0=0;
sdl0=0;
xx1=exp(-x1./ncfs1(1)+ncfs1(2))-1;
xx2=exp(-x2./ncfs2(1)+ncfs2(2))-1;
yy1=(1-2*y1+(xx1+1));
yy2=(-1+2*y2-(xx2+1));

cc=polyfit([xx1;xx2],[yy1;yy2],1);
beta=cc(1);
cfs1(1)=[];
cfs2(1)=[];


for ii=1:maxit


% ncfs1(1)=cfs1(1);
% ncfs1(2)=cfs1(2);

mf1=@(x,p) .5*(1+beta)*(1+(1-beta)/(1+beta)*exp(-x/p(1)+p(2)));
[f,ncfs1,kvg,iter,corp,covp,covr,stdresid,Z,r2]=leasqr(x1,y1,ncfs1,mf1,1e-16,1e3);
ndev1=sqrt(diag(covp)*(pp1-3)/pp1)*tinv(1-(1-confi)/2,pp1-3); 

% ncfs2(1)=cfs2(1);
% ncfs2(2)=cfs2(2);
mf2=@(x,p) .5*(1-beta)*(1+(1+beta)/(1-beta)*exp(-x/p(1)+p(2)));
[f,ncfs2,kvg,iter,corp,covp,covr,stdresid,Z,r2]=leasqr(x2,y2,ncfs2,mf2,1e-16,1e3);
ndev2=sqrt(diag(covp)*(pp2-3)/pp2)*tinv(1-(1-confi)/2,pp2-3); 


l1=ncfs1(1);
l2=ncfs2(1);

sdl=1/sqrt(1/l1^2+1/l2^2);

sdldev=sdl^3*(abs(ndev1(1))/l1^3+abs(ndev2(1))/l2^3);

cfs1(1)=beta;
cfs2(1)=beta;

mf1=@(x,p) .5*(1+p(1))*(1+(1-p(1))/(1+p(1))*exp(-x/l1+p(2)));
[f,cfs1,kvg,iter,corp,covp,covr,stdresid,Z,r2]=leasqr(x1,y1,cfs1,mf1,1e-16,1e3);
dev1=sqrt(diag(covp)*(pp1-3)/pp1)*tinv(1-(1-confi)/2,pp1-3); 

mf2=@(x,p) .5*(1-p(1))*(1+(1+p(1))/(1-p(1))*exp(-x/l2+p(2)));
[f,cfs2,kvg,iter,corp,covp,covr,stdresid,Z,r2]=leasqr(x2,y2,cfs2,mf2,1e-16,1e3);
dev2=sqrt(diag(covp)*(pp2-3)/pp2)*tinv(1-(1-confi)/2,pp2-3); 



bw1=(1/(dev1(1)/cfs1(1)))^pwr;
bw2=(1/(dev2(1)/cfs2(1)))^pwr;
tw=bw1+bw2;
bw1=bw1/(tw);
bw2=bw2/(tw);

beta=cfs1(1)*bw1+cfs2(1)*bw2;
bdev=dev1(1)*bw1+dev2(1)*bw2;

%disp([beta,sdl]);
if ((max(abs(beta0-beta))<prec) && (max(abs(sdl0-sdl))<prec)), break;end;
beta0=beta;
sdl0=sdl;
end
% disp(ncfs1);
% disp(ncfs2);
if (ii==maxit)
disp('Fit for Spin-Diffusion length did not converge!');
exit;
end
cfs=[beta,sdl,ncfs1(1),ncfs2(1),ncfs1(2),ncfs2(2)];

dev=[bdev*2,sdldev*2,ndev1(1),ndev2(1),ndev1(2),ndev2(2)];
bstr=num2str(beta,'%22.15e');
l1str=num2str(ncfs1(1),'%22.15e');
o1str=num2str(ncfs1(2),'%22.15e');
l2str=num2str(ncfs2(1),'%22.15e');
o2str=num2str(ncfs2(2),'%22.15e');


mfustr=['@(x) ', '.5*(1+' , bstr, ')*(1+(1-' , bstr , ')./(1+', bstr,')*exp(-(x/',l1str,')+',o1str,'))'];
#mfu=str2func(qq);
mfdstr=['@(x) ', '.5*(1-' , bstr, ')*(1+(1+' , bstr , ')./(1-', bstr,')*exp(-(x/',l2str,')+',o2str,'))'];
#mfd=str2func(qq);

mfu=eval(mfustr);
mfd=eval(mfdstr);
#mfu=@(x,p) .5*(1+beta)*(1+(1-beta)/(1+beta)*exp(-x/ncfs1(1)+ncfs1(2)));
#mfd=@(x,p) .5*(1-beta)*(1+(1+beta)/(1-beta)*exp(-x/ncfs2(1)+ncfs2(2)));
