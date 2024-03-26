function [cfs,dev,cout]=resfit(x1,rr,pp,wt,cf0)

opts.abstol=1e-8;
confi=.7;
maxit=2000;
y1=rr;
dp=.001*ones(size(cf0));
opt.bounds=[0,Inf;0,Inf;-Inf,0;0,Inf];
pp1=size(x1,1);
mf1=@(x,p) p(1)+p(2)*x+p(3)./(p(4)+x*p(2));
[f,cfs,kvg,iter,corp,covp,covr,stdresid,Z,r2]=leasqr(x1,y1,cf0,mf1,1e-16,maxit,wt,dp,'dfdp');
cout=cfs;
dev=sqrt(diag(covp)*(pp1-3)/pp1)*tinv(1-(1-confi)/2,pp-3); 

a=cfs(1);
b=cfs(2);
c=cfs(3);
d=cfs(4);


rs1=[b,a+c/d,(-a + d)/sqrt(-4*c + (a - d)^2),(-2*c + d*(-a + d))/(sqrt(-4*c + (a - d)^2)*d)];
rs2=[b,a+c/d,(a - d)/sqrt(-4*c + (a - d)^2),(2*c + (a - d)*d)/(sqrt(-4*c + (a - d)^2)*d)];

betadev=2*abs((-a + d)/(-4*c + (a - d)^2)^1.5)*abs(dev(3))+4*abs(c/(-4*c + (a - d)^2)^1.5)*(abs(dev(1))+abs(dev(4)));
gammadev=2*abs((c*(a + d))/((-4*c + (a - d)^2)^1.5*d))*abs(dev(1))+2*abs((2*c + a*(-a + d))/((-4*c + (a - d)^2)^1.5*d))*abs(dev(3))+2*abs((c*(a^2 - 4*c - 3*a*d))/((-4*c + (a - d)^2)^1.5*d^2))*abs(dev(4));
rbdev=abs(dev(1))+1/abs(d)*abs(dev(3))+abs(c/d^2)*abs(dev(4));

cfs=rs2;
dev=[dev(2),rbdev,betadev,gammadev];
endfunction
