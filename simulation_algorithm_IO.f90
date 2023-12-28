Module IO
use CommVar
implicit none
contains

subroutine ReadGrid
	
  grdfile='grid/D-scheme-I.grd'
  open(11,file=grdfile,form='unformatted',status='unknown')
  read(11) nblocks
  print*,  nblocks
  read(11) nx,ny,nz
  print*, nx,ny,nz
  allocate(x(nx,ny,nz), &
  &		   y(nx,ny,nz), &
  &		   z(nx,ny,nz)  )
  read(11) (((x(i,j,k),i=1,nx),j=1,ny),k=1,nz), &
  & 	   (((y(i,j,k),i=1,nx),j=1,ny),k=1,nz), &
  & 	   (((z(i,j,k),i=1,nx),j=1,ny),k=1,nz)
  close(11)
  print*, x(nx,ny,nz),y(nx,ny,nz),z(nx,ny,nz)

end subroutine

subroutine WriteSubGrid
	real*8,allocatable,dimension(:,:,:):: subx,suby,subz
	allocate(subx(1:137,1:120,1:601), &
  &		     suby(1:137,1:120,1:601), &
  &		     subz(1:137,1:120,1:601)  )
	subx(1:137,1:120,1:601) = x(1:137,1:120,360:960)
	suby(1:137,1:120,1:601) = y(1:137,1:120,360:960)
	subz(1:137,1:120,1:601) = z(1:137,1:120,360:960)
	print*, subx(137,120,601),suby(137,120,601),subz(137,120,601)
  open(12,file='SubGrid.grd',form='unformatted')
  write(12) 1
  write(12) 137,120,601
  write(12) (((subx(i,j,k),i=1,137),j=1,120),k=1,601), &
  & 	    (((suby(i,j,k),i=1,137),j=1,120),k=1,601), &
  & 	    (((subz(i,j,k),i=1,137),j=1,120),k=1,601)
  close(12)	
end subroutine	

subroutine ReadField(ts)
  integer ts
  
  write (cts, "(I6)") ts
  fldfile='data/plot3d'//trim(cts)//'.q'
  print*, fldfile
  open(15,file=fldfile,  form='unformatted',status='unknown')
  rewind(15)
  read(15) nx,ny,nz
  print*,  nx,ny,nz
  read(15) mach, alpha, reyn, time
  print*,  mach, alpha, reyn, time
   
  read(15) ((((U(i,j,k,n),i=1,nx),j=1,ny),k=1,nz),n=1,5)	
end subroutine

subroutine ReadLamd2Field(ts)
  integer ts 
  write (cts, "(I6)") ts
  fldfile='lamd2/plot3dLamdTwo'//trim(cts)//'.q'
  print*, fldfile
  open(15,file=fldfile,  form='unformatted',status='unknown')
  rewind(15)
  read(15) nx,ny,nz
  print*,  nx,ny,nz
  ! read(15) mach, alpha, reyn, time
  ! print*,  mach, alpha, reyn, time
   
   allocate(Lamd2(nx,ny,nz,4))
  read(15) ((((Lamd2(i,j,k,n),i=1,nx),j=1,ny),k=1,nz),n=1,4)
  print*, Lamd2(nx,ny,nz,1:4)	
end subroutine

subroutine ReadMeanField
	fldfile='MeanU.uft'
  	print*, "reading ", fldfile
  	open(15,file=fldfile,  form='unformatted',status='unknown')
  	rewind(15)
	read(15) ((((TM_U(i,j,k,n),i=1,nx),j=1,ny),k=1,nz),n=1,5)
	close(15)
	print*, "redading ", fldfile, " finished"

	fldfile='MeanPU.uft'
  	print*, "reading ", fldfile
  	open(15,file=fldfile,  form='unformatted',status='unknown')
  	rewind(15)
	read(15) ((((TM_PU(i,j,k,n),i=1,nx),j=1,ny),k=1,nz),n=1,5)
	close(15)
	print*, "redading ", fldfile, " finished"
end subroutine

subroutine WriteMeanUFT
	open(11,file="MeanU.uft",form='unformatted')
	write(11) TM_U
	close(11)

	open(11,file="MeanPU.uft",form='unformatted')
	write(11) TM_PU
	close(11)
end subroutine
	
subroutine WriteMeanSlice
	 
	open(12,file="Mean_Side.plt")
	write(12,*) 'VARIABLES=y,z,ro,ru,rv,rw,rE,u,v,w,E'
	write(12,*) 'ZONE J=', ny, 'K=', nz, 'F=POINT'
	do k=1,nz
		do j=1,ny
			write(12,"(11(1X,E15.7E3))") y(1,j,k),z(1,j,k), &
	&		TM_U(1,j,k,1:5), TM_PU(1,j,k,2:5)
		end do
	end do
	close(12)

        open(12,file="Mean_Mid.plt")
	write(12,*) 'VARIABLES=y,z,ro,ru,rv,rw,rE,u,v,w,E'
	write(12,*) 'ZONE J=', ny, 'K=', nz, 'F=POINT'
	do k=1,nz
		do j=1,ny
			write(12,"(11(1X,E15.7E3))") y(nx/2,j,k),z(nx/2,j,k), &
	&		TM_U(nx/2,j,k,1:5), TM_PU(nx/2,j,k,2:5)
		end do
	end do
	close(12)
end subroutine

! subroutine WriteFlucSlice

! 	open(12,file="Fluc_Mid.plt")
! 	write(12,*) 'VARIABLES=y,z, u,v,w'
! 	write(12,*) 'ZONE J=', ny, 'K=', nz, 'F=POINT'
! 	do k=1,nz
! 		do j=1,ny
! 			write(12,"(5(1X,E15.7E3))") y(nx/2,j,k),z(nx/2,j,k), &
! 	&		fluc_u(nx/2,j,k,2:4)	
! 		end do
! 	end do
! 	close(12)

! end subroutine

subroutine WriteRSUFT
	open(11,file="RS.uft",form='unformatted')
	write(11) RS
	close(11)
end subroutine

subroutine ReadRS
	fldfile='RS.uft'
  	print*, "reading ", fldfile
  	open(15,file=fldfile,  form='unformatted',status='unknown')
  	rewind(15)
	read(15) ((((RS(i,j,k,n),i=1,nx),j=1,ny),k=1,nz),n=1,4)
	close(15)
	print*, "redading ", fldfile, " finished"
end subroutine

subroutine WriteRSSlice
	open(12,file="RS_Mid.plt")
	write(12,*) 'VARIABLES=y,z,uu,vv,ww,vw'
	write(12,*) 'ZONE J=', ny, 'K=', nz, 'F=POINT'
	do k=1,nz
		do j=1,ny
			write(12,"(6(1X,E15.7E3))") y(nx/2,j,k),z(nx/2,j,k), &
	&		RS(nx/2,j,k,1:4)
		end do
	end do
	close(12)
end subroutine

subroutine WriteSnapshots

end subroutine

subroutine ReadSnapshots

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine ExtractMVG
! integer subnx,subny,subnz
! integer nx1,nx2,nxt
! real*8:: MAXLZ,LOCLZ,Coe
! real*8:: LX1,LX2
! real*8,allocatable,dimension(:,:,:):: subx,suby,subz, xx,yy,zz

! 	subnx = 301
! 	subny = 192
! 	subnz = 137
	
!   allocate(subx(subnx,subny,subnz), &
!   &				 suby(subnx,subny,subnz), &
!   &				 subz(subnx,subny,subnz)  )
  
! 	do k=1,subnz
! 		do j=1,subny
! 			do i=1,subnx
! 				subx(i,j,k) = z(k,j,i+307)
! 				suby(i,j,k) = y(k,j,i+307)
! 				subz(i,j,k) = x(k,j,i+307)
! 			end do
! 		end do
! 	end do
	
! 	!!! Rescale spanwise length
! 	MAXLZ = subz(100,1,subnz)-subz(100,1,1)
! 	print*, MAXLZ
! 	do j=1,subny			
! 		do i=1,subnx
! 			LOCLZ = subz(i,j,subnz)-subz(i,j,1)
! 			do k=1,subnz
! 				subz(i,j,k) = subz(i,j,k)*MAXLZ/LOCLZ
! 			end do
! 		end do
! 	end do
! 	!!! Linearly scale 
! 	Coe=5.13
! 	subx = Coe*subx
! 	suby = Coe*suby
! 	subz = Coe*subz
! 	subx = subx - subx(1,1,1)
! 	!!!!  Add two computational domains
! 	LX1 = 50.0
! 	LX2 = 140.0
! 	nx1 = 150
! 	nx2 = 350
! 	nxt = 801
! 	allocate(xx(1:nxt,1:subny,1:subnz),&
! 					 yy(1:nxt,1:subny,1:subnz),&
! 					 zz(1:nxt,1:subny,1:subnz))
! 	!!! the block before MVG				 
!   do k=1,subnz
!   	do j=1,subny
!   		do i=1,nx1
!   			xx(i,j,k) = LX1*real(i-1)/real(nx1)
!   			yy(i,j,k) = suby(1,j,k)
!   			zz(i,j,k) = subz(1,j,k)
!   		end do
!   	end do
!   end do
!   print*, "block1 finished"
!   !!! MVG block
!   do k=1,subnz
!   	do j=1,subny
!   		do i=1,subnx
! 			  xx(150+i,j,k) = LX1 + subx(i,j,k)
! 			  yy(150+i,j,k) = suby(i,j,k)
! 			  zz(150+i,j,k) = subz(i,j,k) 
! 			end do
! 		end do
! 	end do
!   print*, "MVG block finished"
!   !!! the block after MVG
!   do k=1,subnz
!   	do j=1,subny
!   		do i=1,nx2
!   			xx(451+i,j,k) = xx(451,1,1) + LX2*real(i)/real(nx2)
!   			yy(451+i,j,k) = suby(subnx,j,k)
!   			zz(451+i,j,k) = subz(subnx,j,k)
!   		end do
!   	end do
!   end do  
!   print*, "block2 finished"

!   open(11,file='MVG.grd',form='unformatted')
!   write(11) 1
!   write(11) nxt,subny,subnz
!   write(11) xx,yy,zz
!   close(11)	
  	
! ! 	open(11,file='MVG.grd',form='unformatted')
! !   write(11) 1
! !   write(11) subnx,subny,subnz
! !   write(11) subx,suby,subz
! !   close(11)	
!  end subroutine 
end Module

