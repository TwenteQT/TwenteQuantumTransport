#!/usr/bin/env qoctave

fdir='full';
spotdir='spot';

MS='PA';
fdn=dir([fdir,'/L*']);
for ii=1:length(fdn)
    dn=fdn(ii).name;
    L=str2num(dn(2:end));
    
    for im=1:2
	disp([dn,'/',MS(im)]);
	bfile=[fdir,'/',dn,'/',MS(im),'/tran.mat'];
	if (exist(bfile,'file'))

        load(bfile);    
        kx0=kx1;ky0=ky1;
        nb=size(kx0,1)*size(kx0,2);
        base=zeros(size(kx0));
        cnd(ii,1)=L;
        for ii1=1:4
    	    qq1=tr1{ii1};
    	    qq1(isnan(qq1))=0;
    	    base=base+qq1;
    	    qq1=tr1{ii1};
    	    qq1(isnan(qq1))=0;
    	    base=base+qq1;
    	end
    	base=base/2;
    	sfile=[spotdir,'/',dn,'/',MS(im),'/tran.mat'];
    	if (exist(sfile,'file'))
    	    load(sfile);
    	    if (~isempty(BZ_AREA))
    		lx1=BZ_AREA(1);
    		lx2=BZ_AREA(2);
    		ly1=BZ_AREA(3);
    		ly2=BZ_AREA(4);
    		inds=(kx0>lx1)&(kx0<lx2)&(ky0>ly1)&(ky0<ly2);
    		nr=nnz(inds);
    		base(inds)=0;
    		inds=(kx0>-lx2)&(kx0<-lx1)&(ky0>-ly2)&(ky0<-ly1);
    		base(inds)=0;
    		nr=nr+nnz(inds);
#    		disp([nr,nnz(inds)]);
    		tr=sum(sum(base))/nb;

    		spt=0;
	        for ii1=1:4
    		    qq1=tr1{ii1};
        	    qq1(isnan(qq1))=0;
        	    spt=spt+sum(sum(qq1));
		    qq1=tr1{ii1};
		    qq1(isnan(qq1))=0;
        	    spt=spt+sum(sum(qq1));
    		end
    		dv=(1+(BZ_AREA_SYM==0));
    		spt=spt/2*dv;
    		cfs=nr/(size(kx1,1)*size(kx1,2)*dv);
    		tr=tr+spt*cfs/nb;
    		
    		
    	    end
    	else
    	    tr=sum(sum(base))/nb;
    	end
    	cnd(ii,im+1)=tr;
    end
    end;
%    disp(cnd(ii,:));
end

%    disp(cnd);
    cnd(:,4)=(cnd(:,2)-cnd(:,3))./cnd(:,2);
    cnd(:,5)=(cnd(:,2)-cnd(:,3))./cnd(:,3);
    cnd=sortrows(cnd,1);
    printf('%d  %e  %e  %e  %e\n',cnd');    
    
    save -ascii reslt.txt cnd