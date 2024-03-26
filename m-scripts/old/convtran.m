#!/usr/bin/env qoctave

function [kx,ky,tr]=loadtran(fn,sym)
df=0;
if (strcmp(fn(end-2:end),'.gz')==1)
%    confirm_recursive_rmdir(0)
%    td=tempname;
    q=gunzip(fn);
    fn=q{1};
%    df=1;
end;

fp=fopen(fn); 
if (fp~=-1)
    fgetl(fp);
      [dat,q2]=fscanf(fp,'%*s %*s %le %le %le %le %le %le %le %le %le %le', [10,Inf]);
%      disp(dat(3,:));
      dv=1+(sym==1);
      nn=sqrt(q2/dv/10);
      
      dat(3,dat(3,:)<0)=nan;
      dat(4,dat(4,:)<0)=nan;
      dat(5,dat(5,:)<0)=nan;
      dat(6,dat(6,:)<0)=nan;
      kx=reshape(dat(1,:)',nn,nn*dv);
      ky=reshape(dat(2,:)',nn,nn*dv);
      tr{1}=reshape(dat(3,:)',nn,nn*dv);
      tr{2}=reshape(dat(4,:)',nn,nn*dv);
      tr{3}=reshape(dat(5,:)',nn,nn*dv);
      tr{4}=reshape(dat(6,:)',nn,nn*dv);
      clear dat;
      
else
    disp('error!');
    return;
end

%if (df==1)
%confirm_recursive_rmdir(0);
%rmdir(tempname,'s');
%end

endfunction

try

if (~exist('andout','file')), 
    disp('no ando');
    exit(1); 
end;

fl=dir('tran_*');
if (length(fl) <1)
disp('too few files!');
exit(1);
end;


    BZ_AREA=[];
    BZ_AREA_SYM=0;
    if (exist('output','file')) 
	[ret,res]=system('head -n 120  output | grep "^BZ_AREA"');
    else
	if (exist('output.gz','file')) 
	    [ret,res]=system('zcat output.gz | head -n 120 | grep "^BZ_AREA"');
	    else
	    [ret,res]=system('bzcat output.bz2 | head -n 120 | grep "^BZ_AREA"');
	end
    end	

    if (length(res)>0)
	vals=strread (res,"%s",'delimiter','\n');
	for ii1=1:1:length(vals)
	    res=eval(vals{ii1});
	end
    end;

    
    if (exist('tran_lr','file'))
    [kx1,ky1,tr1]=loadtran('tran_lr',BZ_AREA_SYM);
    else
    if (exist('tran_lr.gz','file'))
    [kx1,ky1,tr1]=loadtran('tran_lr.gz',BZ_AREA_SYM);
    end;
    end;
    
    
    if (exist('tran_rl','file'))
    [kx2,ky2,tr2]=loadtran('tran_rl',BZ_AREA_SYM);
    else
    if (exist('tran_rl.gz','file'))
    [kx2,ky2,tr2]=loadtran('tran_rl.gz',BZ_AREA_SYM);
    end;
    end;
    save('tran.mat','-7','kx1','ky1','kx2','ky2','tr1','tr2','BZ_AREA','BZ_AREA_SYM');
catch 
exit(1);
end 
exit(0);
            
