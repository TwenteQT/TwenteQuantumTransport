#!/usr/bin/env qoctave
cfr=2.581280756014811e-06*1e9;
opt.Robust='on';
opt.confint=0.7;
opts.MaxIter=2e3;
gfact=2.1;
rncut=1;
# To use effective magnetic moment, set mgi=2
mgi=2;
rwp=0;



function [z,ii,jj] = uniquetol(x,tol)
    if size(x,1) == 1, x = x(:); end
        [y,ii,jj] = unique(x);
        if size(x,2) > 1
            [~,ord] = sort(sum(x.^2,1),2,'descend');
            [y,io] = sortrows(y,ord);
            [~,jo] = sort(io);
            ii = ii(io);
            jj = jo(jj);
        end
        d = sum(abs(diff(y,1,1)),2);
        isTol = [true;d > tol];
        z = y(isTol,:);
        bin = cumsum(isTol); % [n,bin] = histc(y,z);
        jj = bin(jj);
        ii = ii(isTol);

    end % UNIQUETOL;



    try
        lbase=load('nlofs');
    catch
        lbase=0;
    end;
    disp(['NLadd=' num2str(lbase)]);
    # Read parameters from geometry
    tic
    qq=ls('-d1','L*');
    l1=sort(str2num(qq(:,2:end)));
    l2=-1;
    for ll=length(l1):-1:1
        try
            fls=ls('-d1',['L',num2str(l1(ll)),'/cf*/transout.h5']);
            if (size(fls,1)>0),
                l2=l1(ll);
                break;
            end
        catch
        end
    end
    if (l2<0),
        disp('no results yet');
        exit
    end

    ifile=strtrim(fls(1,:));
    a0=h5read(ifile,'/geom/lgeo/anm');
    base=h5read(ifile,'/geom/lgeo/base');
    trperp=h5read(ifile,'/geom/lgeo/perp_trans');
    sc=h5read(ifile,'/geom/lgeo/sc_size');
    nl0=h5read(ifile,'/geom/lgeo/num');
    nl1=nl0*sc(1)*sc(2);
    ar=det(base);
    rcfa0=a0^2*ar*cfr/100;
    lscl=a0*trperp(3);
    scsz=sc(1)*sc(2);
    rcfa=scsz*rcfa0;
    cf=1/pi/4*gfact;
    lsh=h5read(ifile,'/cond/l_shc')/rcfa;
    rsh=h5read(ifile,'/cond/r_shc')/rcfa;
    shres=.5*(1./sum(lsh)+1./sum(rsh));
    toc

    rzu=[];
    rzd=[];
    lzu=[];
    lzd=[];
    for i2=1:size(fls,1)
        ifile=strtrim(fls(i2,:));
        try
            flags=h5read(ifile,'/compat');
            rdy=1;
        catch
            rdy=0;
        end
        if (rdy==1)
            if ((flags(5)==1) && (flags(6)==1)),
                Lz=h5read(ifile,'/current/Lz');
                Rz=h5read(ifile,'/current/Rz');
                aucoord=h5read(ifile,'/geom/mgeo/aucoord');
                cz=aucoord(:,3);
                [zz,q1,q2]=uniquetol(cz,1);
                nz=length(zz);
                lz1=zeros(nz,2);
                rz1=lz1;
                for ii=1:nz
                    lz1(ii,:)=sum(Lz(q2==ii,:),1);
                    rz1(ii,:)=sum(Rz(q2==ii,:),1);
                end
                rzu=[rzu,rz1(:,1)];
                rzd=[rzd,rz1(:,2)];
                lzu=[lzu,lz1(:,1)];
                lzd=[lzd,lz1(:,2)];
            end
        end
    end
    mlzu=mean(lzu,2);
    mlzd=mean(lzd,2);    
    mrzu=mean(rzu,2);
    mrzd=mean(rzd,2);        