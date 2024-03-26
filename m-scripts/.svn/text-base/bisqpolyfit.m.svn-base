function [p,dp,opi]=bisqpolyfit(x1,y1,pwr,opts)

    MaxIter=500;
    Tol=1e-16;
%    pin=pin(:);
    dp=0.001*ones(pwr+1,1);
    isRobust=0;
%    Bounds=Inf*[-ones(length(pin(:)),1),ones(length(pin(:)),1)];
    iteroutput=0;
    wts=ones(size(y1));
    confint=0.7;
    exitflag=1;
    if exist('opts')
        if (isfield(opts,'MaxIter'))
            MaxIter=opts.MaxIter;
        end
        if (isfield(opts,'Tol'))
            Tol=opts.Tol;
        end
        if (isfield(opts,'dp'))
            dp=opts.dp;
        end
        if (isfield(opts,'Bounds'))
            Bounds=opts.Bounds;
        end
        if (isfield(opts,'Weights')) 
            wts=opts.Weights;
        end
        if (isfield(opts,'confint'))
            confint=opts.confint;
        end
        if (isfield(opts,'Robust') && strcmpi(opts.Robust,'on'))
            isRobust=1;
        end    
        if (isfield(opts,'Display') && strcmpi(opts.Display,'iter'))
            iteroutput=1;
        end    
    end;

    sqrtwt = sqrt(wts);
    wtdy=y1.*sqrtwt;


    N = length(x1);
    
    % Number of unique X values, to define correctly numbers of freedom degrees (FD)
    Nu=length(unique(x1));

    [p0,S]=wpolyfit(x1,y1,sqrtwt,pwr);f=polyval(p0,x1);
	    
    p = p0(:);
    totaliter = 1;
    res=y1-f;

    NP=length(p);
    resin=res;
    % Formal degrees of freedom
    dfe=N-NP;
    % Real degrees of freedom    
    dfeu=Nu-NP;
    jacob=polyjacob(x1,pwr);
    
    [Q,ignore]=qr(jacob,0);
    h = min(.9999, sum(Q.*Q,2));
    adjfactor = 1 ./ sqrt(1-h);

    % Still unsure which FD one has to use here.
%    sigma = norm(resin) / sqrt(dfeu);
    sigma = norm(resin) / sqrt(dfe);

    tiny_s = 1e-6 * std(wtdy) / sum(wts);
    if tiny_s==0
        tiny_s = 1;
    end

    D = 1e-7;
    robiter = 0;

    origMaxIter=MaxIter;
    iterlim=origMaxIter-1;
    res=resin;
    if (isRobust==1)
        iter1=1;
        while true
            robiter = robiter + 1;

            if iteroutput
                fprintf( '\nRobust fitting iteration %i:\n---------------------------', robiter) ;
            end

            % Compute residuals from previous fit, then compute scale estimate
            radj = res .* adjfactor;
            rs = sort(abs(radj));
            sigma = median(rs(NP:end)) / 0.6745;
   
            % Compute new weights from these residuals, then re-fit
            tune = 4.685;
            bw = cfrobwts(radj/(max(tiny_s,sigma)*tune));
#            bw(1:10);
            p0 = p(:);
            MaxIter = max( 0, iterlim );

	    [p,S]=wpolyfit(x1,y1,1.1./(wts.*bw),pwr);
	    p=p(:);
	    f=polyval(p,x1);
%	    p=p(:);
%	    p1(:)-p0
%            [f,p,exitflag,iter1,corp,covp,covr,stdresid,Z,r2]=...
%            leasqr(x1,y1,p0,mf1,Tol,MaxIter,sqrt(wts.*bw),dp,'dfdp',opt);

            % Reduce the iteration limit by the number of nonlinear fit iterations
            % used for this set of robust weights.  Also keep track of the total
            % number of iterations, as we would like to report that total.
            iterlim = iterlim - iter1;
            totaliter = totaliter + iter1;

            res = wtdy - f.*sqrtwt;

            % Check WHILE loop exit conditions
            if cfInterrupt( 'get' )
                % Fitting stopped by user.
                exitflag = -1;
                break
            elseif (iterlim <= 0)
                % Number of iterations or function evaluations exceeded
                %
                % Set Exit Flag to zero as the limit on number of iterations or
                % function evaluations may have been exceeded without the flag
                % already been set to zero.
                exitflag = 0;
                break
            elseif all( abs( p-p0 ) <= D*max( abs( p ), abs( p0 ) ) )
                % Convergence
                break
            end
        end

        res=wtdy-f.*sqrtwt;
  
        if all(bw<D | bw>1-D)
            % All weights 0 or 1, this amounts to ols using a subset of the data
            included = (bw>1-D);
            % FIXME: This probably has to be corrected for real FD number
            warning('error estimation can be wrong! Correction has to be implemented for this case');
            robust_s = norm(res(included)) / sqrt(sum(included) - NP); 
        else
            % Compute robust mse according to DuMouchel & O'Brien (1989)
            radj = res .* adjfactor;
            robust_s = cfrobsigma(max(tiny_s,sigma),radj, NP, tune, h);
            % Again, corect for "right" number of FDs, not totally sure that it is correct.
%            robust_s=robust_s*sqrt(dfe/dfeu);
        end
        sigma = max(robust_s, sqrt((sigma^2 * NP^2 + robust_s^2 * N) / (NP^2 + N)));
    end;
    % again unclear,dfe or dfeu
%    resnorm = dfeu * sigma^2; % new resnorm based on sigma
    resnorm = dfe * sigma^2; % new resnorm based on sigma
    jacob=polyjacob(x1,pwr);
    [~,R]=qr(jacob,0);
    Rinv = R \ eye(length(R));
    v= sum(Rinv.^2,2)* (sigma^2);
    % Coefficient of student distribution for given confidence interval and NFDs
%    dfeu
%    confint
    t=-tinv((1-confint)/2,dfeu);
%    t=1;
    dp=sqrt(sum(v,2)')*t;
    dp=dp(:);

    opi.isConv=exitflag;
    opi.iter=totaliter;
    opi.robiter=robiter;
    opi.res=res;
    opi.resnorm=resnorm;
    opi.f=f;
endfunction;

function jacob=polyjacob(x1,pwr)
jacob=[];
for ii=0:pwr
jacob=[x1(:).^ii,jacob];
end;
endfunction;

function w = cfrobwts(r)
    w = (abs(r)<1) .* (1 - r.^2).^2;
    if all(w(:)==0)
        w(:) = 1;
    end
endfunction;

function rob_s = cfrobsigma(s,r,p,t,h)
    % Include tuning constant in sigma value
    st = s*t;

    % Get standardized residuals
    n = length(r);
    u = r ./ st;

    % Compute derivative of phi function
    wfun = @(u) cfrobwts(u);
    phi = u .* feval(wfun,u);
    delta = 0.0001;
    u1 = u - delta;
    phi0 = u1 .* feval(wfun,u1);
    u1 = u + delta;
    phi1 = u1 .* feval(wfun,u1);
    dphi = (phi1 - phi0) ./ (2*delta);

    % Compute means of dphi and phi^2; called a and b by Street.  Note that we
    % are including the leverage value here as recommended by O'Brien.
    m1 = mean(dphi);
    m2 = sum((1-h).*phi.^2)/(n-p);

    % Compute factor that is called K by Huber and O'Brien, and lambda by
    % Street.  Note that O'Brien uses a different expression, but we are using
    % the expression that both other sources use.
    K = 1 + (p/n) * (1-m1) / m1;

    % Compute final sigma estimate.  Note that Street uses sqrt(K) in place of
    % K, and that some Huber expressions do not show the st term here.
    rob_s = K*sqrt(m2) * st /(m1);

endfunction;

function stop = cfInterrupt(action, stop)

    persistent PERSISTENT_STOP
    if isempty( PERSISTENT_STOP )
        PERSISTENT_STOP = false;
    end

    switch action
    case 'get'
        stop = PERSISTENT_STOP;
    case 'set'
        PERSISTENT_STOP = stop;
    otherwise
        error(message('curvefit:cfInterrupt:InvalidAction'));
    end

end
