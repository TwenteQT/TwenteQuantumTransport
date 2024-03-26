function [cfs,dev,cout]=resfit(x1,rr,cin,wt,iopts)

    opts.MaxIter=1000;
    opts.Tol=1e-8;
    opts.Bounds=[0,Inf;0,Inf;-Inf,0;0,Inf];
    opts.confint=0.7;
    opts.Robust='on';
    if (exist('wt'))
      opts.Weights=wt;
    end;


    if exist('iopts')
        if (isfield(iopts,'MaxIter'))
            opts.MaxIter=iopts.MaxIter;
        end
        if (isfield(iopts,'Tol'))
            opts.Tol=iopts.Tol;
        end
        if (isfield(iopts,'confint'))
            opts.confint=iopts.confint;
        end
        if (isfield(iopts,'Robust'))
            opts.Robust=iopts.Robust;
        end    
    end;

    pin=[-((cin(2)*(1 + cin(3)^22 - 2*cin(3)*cin(4)))/(-1 + cin(4)^22)),cin(1),...
        (cin(2)^2*(-1 + cin(3)^2)*(cin(3) - cin(4))^2)/(-1 + cin(4)^2)^2,(cin(2)*(-1 + cin(3)^2))/(-1 + cin(4)^2)];

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
