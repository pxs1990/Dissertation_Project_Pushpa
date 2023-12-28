! ======================================================================
!    FileName: exportData.f90
!     Project: exportData
!      Author: Yong Yang, Department of Mathematics, UTA
!       Email: yong.yang@mavs.uta.edu
!     Created: 2017-04-06 16:19:11
! Modified by: Pushpa Shrestha, Department of Mathematics, UTA
!       Email: pushpa.shrestha@mavs.uta.edu
! Last Change: 2019-09-05       
! ======================================================================
program main
    implicit none
	integer :: i,j,k,m,imax,jmax,kmax,mmax,n,dm
    integer :: istart,iend,jstart,jend,kstart,kend,nstart,nend,mstart,mend
	integer :: alloc_err
	real 	,dimension(:,:,:,:), allocatable :: Q
	logical :: fileExist
	character(30) :: fname,qname,cn
    character(15) :: string
    integer :: ios
    character(100) :: msg
    integer :: n1, n2, n3

    mmax=5
    open(11,file='input.dat')
    read(11,*) nstart !=1751
    read(11,*) nend !=1950
    read(11,*) istart !=101
    read(11,*) iend !=200
    read(11,*) jstart !=1
    read(11,*) jend !=128
    read(11,*) kstart !=1
    read(11,*) kend !=200
    read(11,*) mstart
	read(11,*) mend
    close(11)
   fname='plot3d296500.q'
	inquire(file=trim(adjustl(fname)),exist=fileExist)
	  if (.not. fileExist) then
		   write(*,*) 'File '//trim(adjustl(fname))//' is not ready! Please check it!'
		   write(*,*) 'Exiting...'
		   stop
      end if

    open(12,file=trim(adjustl(fname)),status='old',form='unformatted', action='read',       &
       iostat=ios, iomsg=msg)
	    if(ios /= 0) then
           write(*,*) 'ERROR: ', msg
           stop
        end if	
	read(12, iostat=ios, iomsg=msg) imax,jmax,kmax
	    if(ios /= 0) then
           write(*,*) 'ERROR: ', msg
           stop
        end if
	read(12, iostat=ios, iomsg=msg) n1, n2, n3
	    if(ios /= 0) then
           write(*,*) 'ERROR: ', msg
           stop
        end if
	close(12)

    write(*,*) 'Allocate arrays...'
	allocate(Q(imax,jmax,kmax,mmax),STAT=alloc_err)
	if (alloc_err/=0) then
		write(*,*) 'Q allocate error'
		stop
    end if
	write(*,*) 'Q Allocate successful'  
	
	write(cn,'(I6)') nstart
	    
		qname='plot3d'//trim(adjustl(cn))//'.q'
		inquire(file=trim(adjustl(qname)),exist=fileExist)
		if (.not. fileExist) then
			write(*,*) 'File '//trim(adjustl(qname))//' is not ready! Please check it!'
			write(*,*) 'Exiting...'
			stop
		end if
		write(*,*) 'reading'
		open(13,file=trim(adjustl(qname)),form='unformatted', action='read',       &
       iostat=ios, iomsg=msg)
		  if(ios /= 0) then
             write(*,*) 'ERROR: ', msg
             stop
          end if
		read(13, iostat=ios, iomsg=msg) imax,jmax,kmax
		  if(ios /= 0) then
             write(*,*) 'ERROR: ', msg
             stop
          end if
		read(13, iostat=ios, iomsg=msg) n1,n2,n3
		  if(ios /= 0) then
             write(*,*) 'ERROR: ', msg
             stop
          end if   
		read(13, iostat=ios, iomsg=msg) ((((Q(i,j,k,m),i=1,imax),j=1,jmax),k=1,kmax),m=1,mmax)
		  if(ios /= 0) then
             write(*,*) 'ERROR: ', msg
             stop
          end if
		close(13)
        
		write(*,*) 'writting'
        fname='export'//'_whole_'//trim(adjustl(cn))//'.dat'
        open(14,file=trim(adjustl(fname)),form='unformatted')
		write(14) iend-istart+1,jend-jstart+1,kend-kstart+1,mend-mstart+1
	    write(14) ((((Q(i,j,k,m),i=istart,iend),j=jstart,jend),k=kstart,kend),m=mstart,mend)
		close(14)
end program main
