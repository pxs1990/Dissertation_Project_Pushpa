Module CommCal
use CommVar
use IO
implicit none

contains

subroutine TimeMeanCal
	integer l

	nstart = 380000
	nend   = 490000
 	do ts = nstart,nend,500
 		call ReadField(ts)
		print*, ts, "  finished"
		!!!!!! Mean of Conservative Var 
		TM_U = TM_U + U

		!!!!!! Mean of Primitive Var
		do k=1,nz
			do j=1,ny
				do i=1,nx
					do l=2,4
						TM_PU(i,j,k,l) = TM_PU(i,j,k,l) + U(i,j,k,l)/U(i,j,k,1)
					end do
				end do
			end do
		end do

 	end do
 	nsample = (nend-nstart)/500 + 1
	TM_U = TM_U/nsample
	TM_PU = TM_PU/nsample			

end Subroutine

! subroutine FindMaxMinPos
! 	integer a(3), b(3)
! 	a = Maxloc(U(1:nx,1:ny,1:nz,1))
! 	b = Minloc(U(1:nx,1:ny,1:nz,1))
! 	
! 	print*, a, U(a(1),a(2),a(3),1)
! 	print*, b, U(b(1),b(2),b(3),1)
! end subroutine

subroutine FlucCal(ts,Fluc_U)

	integer ts
	integer l
	real*8,allocatable,dimension(:,:,:,:):: Fluc_U

	! allocate(Fluc_U(nx,ny,nz,5))
	Fluc_U(:,:,:,:) = 0.0d0

	call ReadField(ts)

	! Fluc_U(:,:,:,1) = U(:,:,:,1) - TM_U(:,:,:,1)

	do l=2,5
		Fluc_U(:,:,:,l) = U(:,:,:,l)/U(:,:,:,1) - TM_PU(:,:,:,l)
	end do

	print*, "Fluctuation field of ", ts , "finished ..."
end subroutine

!!!!!! 2 order statistics !!!!!!!!!
subroutine RSCal
	real*8,allocatable,dimension(:,:,:,:):: Fluc_U
	allocate(Fluc_U(nx,ny,nz,5))
	nstart = 380000
	nend   = 490000
 	do ts = nstart,nend,500
 		call FlucCal(ts,Fluc_U)
 		RS(:,:,:,1) = RS(:,:,:,1) + Fluc_U(:,:,:,2)**2
 		RS(:,:,:,2) = RS(:,:,:,2) + Fluc_U(:,:,:,3)**2
 		RS(:,:,:,3) = RS(:,:,:,3) + Fluc_U(:,:,:,4)**2
 		RS(:,:,:,4) = RS(:,:,:,4) + Fluc_U(:,:,:,3)*Fluc_U(:,:,:,4)
 	end do
 	nsample = (nend-nstart)/500 + 1
	RS = RS/nsample
end subroutine

!!!!!!!!!!!!!!!!!!!!! POD  &  DMD  !!!!!!!!!!!!!!!!
subroutine CorvMCal

	integer NS, snx,sny,snz,snp
	integer ts, m
	real*8,allocatable,dimension(:,:):: CorvM
	real*8,allocatable,dimension(:,:,:,:):: fluc_vel,samp_fluc
	! real*8,allocatable,dimension(:,:) :: Snapshots

	NS = 60

	snx = 137
	sny = 120
	snz = 601
	snp = snx*sny*snz*3

	allocate(fluc_vel(nx,ny,nz,5))
	allocate(samp_fluc(snx,sny,snz,3))
	allocate(Snapshots(snp,NS))
	allocate(CorvM(NS,NS))
	CorvM = 0.0


	do m=1,NS
		ts	= 380000 + 1500*(m-1)			
		call FlucCal(ts, fluc_vel)
		samp_fluc = fluc_vel(1:137, 1:120, 360:960, 2:4)
		Snapshots(:,m) = reshape(samp_fluc,(/snp/))
	end do

	open(11,file="Snapshots.uft",form='unformatted')
		write(11) Snapshots
	close(11)

	do j=1,NS
		do i=1,j
			CorvM(i,j) = dot_product(Snapshots(:,i),Snapshots(:,j))
			CorvM(j,i) = CorvM(i,j)
		end do
	end do		

	CorvM = CorvM/snp

	print*, "the Corvariance Matrix is : "
	print*, CorvM

	open(11,file="CorvM.txt")
	do j=1,NS
		do i=1,NS
			write(11,*) CorvM(i,j)
		end do
	end do

end subroutine

subroutine PODModesCal
	integer NS, snx,sny,snz,snp
	integer l, m, ii
	real*8,allocatable,dimension(:) :: EigValue
	real*8,allocatable,dimension(:,:) :: EigVector
	real*8,allocatable,dimension(:,:) :: POD
	real*8,allocatable,dimension(:) :: Tmp_Snap
	real*8,allocatable,dimension(:,:,:,:):: Tmp_U
	character*80 tmp_char

	NS = 60

	allocate(EigValue(NS),EigVector(NS,NS))	

	snx = 137
	sny = 120
	snz = 601
	snp = snx*sny*snz*3

	allocate(Snapshots(snp,NS),POD(snp,NS),Tmp_Snap(snp))
	allocate(Tmp_U(snx,sny,snz,3))


	!!!!!!!!! Read Snapshots of Fluc field
	fldfile='Snapshots.uft'
  	print*, "reading ", fldfile
  	open(15,file=fldfile,  form='unformatted',status='unknown')
  	rewind(15)
	read(15) ((Snapshots(l,m),l=1,snp),m=1,NS)
	close(15)
	print*, "redading ", fldfile, " finished"

	! Snapshots = 1.0
	! print*, size(Snapshots(:,1)), Snapshots(29641320,1)
	! Tmp_Snap = Snapshots(:,1)
	! Tmp_U = Reshape(Tmp_Snap, (/snx,sny,snz,3/))
	! ii=1
	! do n=1,3
	! 	do k=1,snz
	! 		do j=1,sny
	! 			do i=1,snx
	! 				Tmp_U(i,j,k,n) = Snapshots(ii,1)
	! 				ii = ii+1
	! 			end do
	! 		end do
	! 	end do
	! end do
	! print*, ii


	!!!!!!!!! Read Eig System
	open(15,file='EigValue.txt')
  	rewind(15)
  	do i=1,NS
		read(15,*) EigValue(i)
	end do
	close(15)
	! print*, EigValue

	open(15,file='EigVector.txt')
  	rewind(15)
  	do j=1,NS
	  	do i=1,NS
			read(15,*) EigVector(i,j)
		end do
	end do
	close(15)
	! print*, EigVector(:,1)

	!!!!!!!!! Calculate POD Modes  
	! POD = MATMUL(Snapshots,EigVector)

	POD = 0.0	
	do m = 3,6
		do j = 1,NS
			POD(:,m) = POD(:,m) + Snapshots(:,j)*EigVector(j,m)
		end do
		print*, norm2(POD(:,m))
		POD(:,m) = POD(:,m)/norm2(POD(:,m))
	end do 	


	!!!!!!!  Write POD modes
	do m=3,6
		write (tmp_char, "(I1)") m
		fldfile = 'POD_'//trim(tmp_char)//'.plt'
		! open(15,file=fldfile,  form='unformatted')
		! write(15) snx,sny,snz
		open(12,file=fldfile)
		write(12,*) 'VARIABLES=x,y,z,u,v,w'
		write(12,*) 'ZONE I=', snx , 'J=', sny, 'K=', snz, 'F=POINT'
		! write(15) mach, alpha, reyn, time
		!!!!!!! Reshape the POD Fields
		Tmp_Snap = POD(:,m)
		Tmp_U = Reshape(Tmp_Snap,(/snx,sny,snz,3/))
		print*, Tmp_U(snx,sny,snz,1:3)

		do k=1,snz
			do j=1,sny
				do i=1,snx
					write(12,"(6(1X,E15.7E3))") x(i,j,k+359),y(i,j,k+359), &
						z(i,j,k+359), Tmp_U(i,j,k,1:3)
				end do
			end do
		end do
		! write(15) ((((Tmp_U(i,j,k,n),i=1,snx),j=1,sny),k=1,snz),n=1,3)
		close(12)
	end do
	!!!!!!!	 Reconstruct Flowfield
end subroutine



!!!!!!!!!!!!!!! 2-D POD modes 
subroutine CorvMCal_2D

	integer NS, snx,sny,snz,snp
	integer ts, m
	real*8,allocatable,dimension(:,:):: CorvM
	real*8,allocatable,dimension(:,:,:,:):: fluc_vel
	real*8,allocatable,dimension(:,:,:)::   samp_fluc
	! real*8,allocatable,dimension(:,:) :: Snapshots

	NS = 2

	sny = 120
	snz = 601
	snp = sny*snz*3

	allocate(fluc_vel(nx,ny,nz,5))
	allocate(samp_fluc(sny,snz,3))
	allocate(Snapshots(snp,NS))
	allocate(CorvM(NS,NS))
	CorvM = 0.0


	do m=1,NS
		ts	= 380000 + 500*(m-1)			
		call FlucCal(ts, fluc_vel)
		samp_fluc = fluc_vel(68, 1:120, 360:960, 2:4)
		Snapshots(:,m) = reshape(samp_fluc,(/snp/))
	end do

	open(11,file="Snapshots.uft",form='unformatted')
		write(11) Snapshots
	close(11)

	do j=1,NS
		do i=1,j
			CorvM(i,j) = dot_product(Snapshots(:,i),Snapshots(:,j))
			CorvM(j,i) = CorvM(i,j)
		end do
	end do		

	CorvM = CorvM/snp

	print*, "the Corvariance Matrix is : "
	print*, CorvM

	open(11,file="2DCorvM.txt")
	do j=1,NS
		do i=1,NS
			write(11,*) CorvM(i,j)
		end do
	end do

end subroutine


!!!!!!!!!!!!!!!  K-H instability wave-length 
subroutine KHWaveLen
	real*8 :: Ly,Lz,dy,dz
	integer :: jj,kk,ll(2)
	real*8,allocatable,dimension(:,:) :: Cor

	allocate(Cor(0:30,0:60))

	call ReadLamd2Field(400000)
	Ly = 3.0
	Lz = 6.0

	Cor = 0.0

	do kk=0,1
		do jj=0,1
			dz = Lz*kk/60
			dy = Ly*jj/30
			do k=600,900
			do j=1,150
			! ll=minloc( abs(z(nx/2,:,:)-(z(nx/2,j,k)+dz)) + &
			! 		   abs(y(nx/2,:,:)-(y(nx/2,j,k)+dy)) )					
			Cor(jj,kk) = Cor(jj,kk) + Lamd2(nx/2,j,k,1)*Lamd2(nx/2,j+jj,k+kk,1)
			end do
			end do
		end do
	end do

	open(11,file="Cor.txt")
	do kk=0,60
		do jj=0,30
			write(11,*) Cor(jj,kk)
		end do
	end do
	close(11)
end subroutine

end Module