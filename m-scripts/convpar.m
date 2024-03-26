#!/usr/bin/env qoctave

pf=fopen('params.txt');
v=fread(pf,Inf,'uchar=>uchar')';
fclose(pf);
pars=strsplit(strrep(v,'#','%'),char(10));
for ii=1:length(pars)
    if strfind(pars{ii},'=')
        eval([pars{ii},';']);
    end
end;

clear pf v pars ii
save -hdf5 params.h5 
