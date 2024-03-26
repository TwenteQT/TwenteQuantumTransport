function [cfs,dev,mfu,mfd]=ndiffit2(x1,gu,gd,gm,pp,pwr)
    opts.abstol=1e-6;
    gt=gu+gd+2*gm;
    gut=gu+gm;
    gdt=gd+gm;
    y1=gut./gt;
    wt=x1.^(pwr/2);
    confi=.7;
    maxit=2000;
    cfs(2)=y1(end)*2-1;
    cfs(1)=5;
    cfs(3)=0;
    pp1=size(x1,1);
    mf1=@(x,p) 0.5+p(2)/2*(1-exp(-x/p(1)+p(3)));
    [f,cfs,kvg,iter,corp,covp,covr,stdresid,Z,r2]=leasqr(x1,y1,cfs,mf1,1e-16,1000,wt);
    dev=sqrt(diag(covp)*(pp1-3)/pp1)*tinv(1-(1-confi)/2,pp-3); 

    bstr1=num2str(cfs(2),'(%22.15e)');
    l1str=num2str(cfs(1),'(%22.15e)');
    o1str=num2str(cfs(3),'(%22.15e)');


    mfustr=['@(x) ', '0.5+0.5*' , bstr1, '*(1-exp(-(x/',l1str,')+',o1str,'))'];
    mfdstr=['@(x) ', '0.5-0.5*' , bstr1, '*(1-exp(-(x/',l1str,')+',o1str,'))'];
    mfu=eval(mfustr);
    mfd=eval(mfdstr);
end;