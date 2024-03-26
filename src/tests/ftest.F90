program qq
implicit none
integer*8 :: fid
integer*8,external::c_fopen
real(kind=8):: vv(4)=(/1.,2.,3.,4./), vv1(4)
    fid=c_fopen('this-file',len('this-file') ,'w')
    call c_fwrite(vv,8*4,1,fid)
    call c_fclose(fid)
    
    fid=c_fopen('this-file',len('this-file') ,'r')
    call c_fread(vv1,8*4,1,fid)    
    write(*,*) vv1
    call c_fclose(fid)
    write(*,'(A)') 'wwww'
end program qq