! ======================================================================
!    FileName: exportGrid.f90
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
	integer :: i,j,k,imax,jmax,kmax,ns,nstart,nend
    integer :: istart,iend,jstart,jend,kstart,kend
	integer :: alloc_err
	!real 	:: rho,u,v,w
	real 	,dimension(:,:,:), allocatable :: X,Y,Z

	logical :: fileExist
	character(30) :: qname,grdname,cn

    open(11,file='input.dat')
    read(11,*) nstart 
    read(11,*) nend !=1950
    read(11,*) istart !=101
    read(11,*) iend !=200
    read(11,*) jstart !=1
    read(11,*) jend !=128
    read(11,*) kstart !=1
    read(11,*) kend !=200

    close(11)
    
    grdname="3d.grd"
	inquire(file=trim(adjustl(grdname)),exist=fileExist)
	if (.not. fileExist) then
		write(*,*) 'File '//trim(adjustl(grdname))//' is not ready! Please check it!'
		write(*,*) 'Exiting...'
		stop
    end if

    open(11,file=trim(adjustl(grdname)),status='old',form='unformatted')
	read(11) 
	read(11) imax,jmax,kmax
    write(*,*) 'Allocate arrays...'
	allocate(X(imax,jmax,kmax),STAT=alloc_err)
	if (alloc_err/=0) then
		write(*,*) 'X allocate error'
		stop
	end if
	allocate(Y(imax,jmax,kmax),STAT=alloc_err)
	if (alloc_err/=0) then
		write(*,*) 'Y allocate error'
		stop
	end if
	allocate(Z(imax,jmax,kmax),STAT=alloc_err)
	if (alloc_err/=0) then
		write(*,*) 'Z allocate error'
		stop
	end if
	write(*,*) 'Allocate successful'
	write(*,*) 'Reading grid file...'

	read(11) &
	& ((( X(i,j,k), i=1,imax), j=1,jmax), k=1,kmax),&
	& ((( Y(i,j,k), i=1,imax), j=1,jmax), k=1,kmax),&
	& ((( Z(i,j,k), i=1,imax), j=1,jmax), k=1,kmax)
	close(11)

	open(12,file=trim(adjustl(grdname//'gridA2')),form='unformatted')
	write(12) iend-istart+1,jend-jstart+1,kend-kstart+1
	write(12) &
	& ((( X(i,j,k), i=istart,iend), j=jstart,jend), k=kstart,kend),&
	& ((( Y(i,j,k), i=istart,iend), j=jstart,jend), k=kstart,kend),&
	& ((( Z(i,j,k), i=istart,iend), j=jstart,jend), k=kstart,kend)
	close(12)




end program main
