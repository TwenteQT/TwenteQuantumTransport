function [cfs,dev,mfu,mfd]=diffit2(x1,y1,x2,y2,prec,pwr)

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


l1=cfs1(1);
l2=cfs2(1);

sdl=1/sqrt(1/l1^2+1/l2^2);

sdldev=sdl^3*(abs(dev1(1))/l1^3+abs(dev2(1))/l2^3);

beta1=cfs1(2);
beta2=cfs2(2);

bw1=(1/(dev1(2)/beta1))^pwr;
bw2=(1/(dev2(2)/beta2))^pwr;
tw=bw1+bw2;
bw1=bw1/(tw);
bw2=bw2/(tw);

beta=beta1*bw1+beta2*bw2;
bdev=dev1(2)*bw1+dev2(2)*bw2;


cfs=[beta,sdl,l1,l2,beta1,beta2];

dev=[bdev*2,sdldev*2,dev1(1),dev2(1),dev1(2),dev2(2)];
bstr1=num2str(beta1,'(%22.15e)');
bstr2=num2str(beta2,'(%22.15e)');

l1str=num2str(cfs1(1),'(%22.15e)');
o1str=num2str(cfs1(3),'(%22.15e)');
l2str=num2str(cfs2(1),'(%22.15e)');
o2str=num2str(cfs2(3),'(%22.15e)');


mfustr=['@(x) ', '.5*(1+' , bstr1, ')*(1+(1-' , bstr1 , ')./(1+', bstr1,')*exp(-(x/',l1str,')+',o1str,'))'];
#mfu=str2func(qq);
mfdstr=['@(x) ', '.5*(1-' , bstr2, ')*(1+(1+' , bstr2 , ')./(1-', bstr2,')*exp(-(x/',l2str,')+',o2str,'))'];
#mfd=str2func(qq);

mfu=eval(mfustr);
mfd=eval(mfdstr);
#mfu=@(x,p) .5*(1+beta)*(1+(1-beta)/(1+beta)*exp(-x/ncfs1(1)+ncfs1(2)));
#mfd=@(x,p) .5*(1-beta)*(1+(1+beta)/(1-beta)*exp(-x/ncfs2(1)+ncfs2(2)));
