module parameters
implicit none

	integer:: natom  !total atom number
	integer :: iatom,jatom
	integer, parameter :: dq=kind(1.0D0)
	integer :: nfram   !total t0,100ps
	integer :: ifram
	integer:: nwater      !total number of water
	integer :: iwater
	real(kind=dq) :: boxL   !the length of box , in anstrong
	character(len=50) :: temp
	integer,allocatable  :: NumMol(:),NumAtom(:)
	character(len=5),allocatable :: typeMol(:),typeAtom(:)  !atom type

	integer ii
	integer, parameter :: ndimention=3
	integer:: idimention
	real, allocatable :: posmd(:,:,:)   !huge array, need large memory

	!H-bond definition
	real(kind=dq):: R_OO          !in anstrong
	real(kind=dq):: theta_HOO       !in degree
	real(kind=dq):: theta
	real(kind=dq), parameter :: pi=3.1415926
	integer :: iHjump   !count jump times of every H atom
	real(kind=dq) :: r1(3),r2(3),roo,roh(2),thetaHOO(2)
	integer, allocatable :: nhbond(:,:)  !每一帧中每个H原子形成氢键的数目
	integer, allocatable :: hbond(:,:,:) !每一帧中与每个氢原子形成nhbond个氢键的氧原子的序号
	
	!search the stable state
	integer nf          !保存稳态发生跳跃的次数
	integer :: iHatom,nHatom  !save H atoms number   
	integer,allocatable :: nSS(:)
	integer,allocatable :: t0(:,:),t1(:,:),indexO(:,:)
 
	
	!calculate CRP in SSP
	integer :: nt  !时间截断值
	integer,allocatable :: nHjump(:)
	integer :: nswitch
	integer,allocatable :: tj(:)
	integer it,iswitch
	integer corr        !
	real(kind=dq),allocatable :: CRP(:)
	real:: tstep   !time step for saving
	character(len=40)::Package
	character(len=40)::gromacsfile
	character(len=40)::tinkerfile
	character(len=8):: inputfile='InputCar'
	
	namelist /input/ natom,nfram,nwater,R_OO,theta_HOO,nt,tstep,Package,gromacsfile,tinkerfile
	contains
	
	subroutine readinput()
	implicit none
		open(unit=12,file=inputfile,status='old')
		read(12,nml=input)
		close(12)
		nf=int(nfram/2)
		allocate(NumMol(natom),NumAtom(natom))
		allocate(typeMol(natom),typeAtom(natom))
		allocate(posmd(ndimention,natom,nfram))
		allocate(nhbond(2*nwater,nfram))      !大型数组
		allocate(hbond(10,2*nwater,nfram))    !设定一个氢原子最多能同时形成10个氢键，大型数组
		allocate(nHjump(2*nwater))
		allocate(nSS(2*nwater))
		allocate(t0(nf,2*nwater))
		allocate(t1(nf,2*nwater))
		allocate(indexO(nf,2*nwater))
		allocate(CRP(nt+1))

	end subroutine readinput
	
	subroutine readgro()
	use readxtc_read
	implicit none

		open (unit=11, file=trim(adjustl(gromacsfile)),status='unknown')
		read(11,*)temp
		read(11,*)temp  
		Do iatom=1, natom
			read (11,"(I5,2A5,I5)")NumMol(iatom),typeMol(iatom),typeAtom(iatom),NumAtom(iatom)
			!write (*,"(I5,2A5,I5)")NumMol(iatom),typeMol(iatom),typeAtom(iatom),NumAtom(iatom)
		End do
			read(11,'(f10.5)') boxL
			boxL = boxL*10.0
			rewind(11)

		!read position of every atom in all frame
		!allocate(posmd(ndimention,natom,nfram))
		posmd=0.0
		call read_xtcp(ndimention,natom,nfram,posmd)            
		posmd=posmd*10.0 
		close(11)
	end subroutine readgro
	
	subroutine readtinker()
	implicit none
		integer::atomtype
		!allocate(typeAtom(natom))
		open(unit=11,file=trim(adjustl(tinkerfile)),status='unknown')
		do ifram=1,nfram
			read(11,*) temp
			do iatom=1,natom
				read(11,"(11X,3F12.6,I6)") (posmd(idimention,iatom,ifram),idimention=1,ndimention),atomtype
				if(atomtype==1) then
					typeAtom(iatom)='   OW'
				elseif(atomtype==2) then
					typeAtom(iatom)='   HW'
				endif
			enddo
			read(11,"(F10.5)") boxL
		enddo
		close(11)
	end subroutine readtinker
	
	subroutine findHbond()
	implicit none
	theta=theta_HOO*pi/180.0   !change the unit to radian	
	nhbond=0
	hbond=0

	Do ifram=1,nfram
		iHatom=0
		DO iatom=1,natom
			If (typeAtom(iatom)=='   OW') then
				Do ii=1,2				!every water have two H atom
					iHatom=iHatom+1    !记录氢原子序号
					do jatom=1,natom
						if ((iatom/=jatom) .AND. typeAtom(jatom)=='   OW') then
							r1(:)=posmd(:,jatom,ifram)-posmd(:,iatom,ifram)    !O*Oa vector
							call pbc(r1,boxL)
							roo=sqrt(r1(1)**2+r1(2)**2+r1(3)**2)                !O*Oa scalar
							!write(*,"(i8,i8,i8,f9.3)")ifram,iatom,jatom,roo
							if(roo < R_OO) then								         		
								r2(:)=posmd(:,iatom+ii,ifram)-posmd(:,iatom,ifram)    !O*H vector
								call pbc(r2,boxL)
								roh(ii)=sqrt(r2(1)**2+r2(2)**2+r2(3)**2)               !O*H scalar
								thetaHOO(ii)=acos((r1(1)*r2(1)+r1(2)*r2(2)+r1(3)*r2(3))/(roo*roh(ii)))   ! theta H*O*Oa
								if(thetaHOO(ii) < theta) then
									nhbond(iHatom,ifram)=nhbond(iHatom,ifram)+1    !第ifram帧第iHatom个氢原子形成nhbond个氢键
									hbond(nhbond(iHatom,ifram),iHatom,ifram)=jatom     !与第ifram帧第iHatom个氢原子形成nhbond个氢键的氧原子序号									
								end if											
							end if						
						end if					
					end do  !jatom
					!write(301,"(i8,i8,i8)")ifram,iHatom,nhbond(iHatom,ifram)
				End do  !ii
			End if
		End do  !iatom
	End do   !ifram
	!deallocate(posmd)
	end subroutine findHbond
	
	subroutine findSS()
	implicit none
	t0=0
	t1=0
	indexO=0
	nHatom = 2*nwater
	Do iHatom=1,nHatom 
		nSS(iHatom)=0           !记录稳态数目
		do ifram=1,nfram
			If(nhbond(iHatom,ifram) == 1) then    !第ifram帧第iHatom的氢原子只形成一个氢键，稳态
				if(nSS(iHatom)==0) then  !表示找到第一个稳态
					nSS(iHatom)=nSS(iHatom)+1
					t0(nSS(iHatom),iHatom)=ifram           !稳态起始时刻对应的帧
					indexO(nSS(iHatom),iHatom)=hbond(1,iHatom,ifram)  !稳态对应的氧原子序号，Oa			
				else
					if(indexO(nSS(iHatom),iHatom) == hbond(1,iHatom,ifram)) then
						t1(nSS(iHatom),iHatom)=ifram           !稳态终止时刻对应的帧
					elseif(indexO(nSS(iHatom),iHatom) /= hbond(1,iHatom,ifram)) then
						nSS(iHatom)=nSS(iHatom)+1                             !找到了新的稳态
						t0(nSS(iHatom),iHatom)=ifram                    !新稳态起始时刻对应的帧
						indexO(nSS(iHatom),iHatom)=hbond(1,iHatom,ifram)   !新稳态对应的氧原子序号，Ob
					endif					
				end if						
			End if					
		end do !ifram			
	End do  !iHatom
	end subroutine findSS
	
	subroutine Cal_SSP_CF()
	implicit none
		open(unit=20,file='jump-frame.out',status='unknown')
		!open(unit=21,file='corr.out',status='unknown')
		write(20,*) 'jump correlation function for water-water:'
		!!由nSS（iHatom）确定nHjump(iHatom)
		nswitch=0   !统计switch次数
		Do iHatom=1,nHatom
			if (nSS(iHatom) == 0 ) then 
				nHjump(iHatom)=1		
			elseif (nSS(iHatom) == 1) then
				nHjump(iHatom)=2
			else
				nHjump(iHatom)=nSS(iHatom)-1
			end if
			nswitch=nswitch+nHjump(iHatom)        !记录发生switch的总数
		End do
		write(20,*)nswitch
		allocate(tj(nswitch))    !这个数组只能在这里allocate
		tj=0
		iswitch=0	
		Do iHatom=1,nHatom
				Do iHjump=1,nHjump(iHatom)
					iswitch=iswitch+1
					if(nSS(iHatom) == 0) then
						tj(iswitch)=nfram
					elseif(nSS(iHatom) == 1) then 
						tj(iswitch)=int(nfram/2)
					else 
						tj(iswitch)=t0(iHjump+1,iHatom)-t0(iHjump,iHatom)
					endif
				End do
		ENd do

		CRP=0.0
		Do it=0,nt
			corr=0
			Do iswitch=1,nswitch
				if(it <= tj(iswitch))then
					corr=corr+1
				end if
			End do
			!write(21,"(2(i8,2x))") corr,nswitch
			CRP(it)=real(corr)/real(nswitch)
			write(20,"(2(f10.5,2x))") real(it)*tstep,CRP(it)
		End do
		close(20)
	end subroutine Cal_SSP_CF
	
	subroutine cleanmemory()
	implicit none
	write(*,*)'OK! 6'
		deallocate(posmd)
		write(*,*)'OK! 7'
		!deallocate(nSS)
		!write(*,*)'OK! 12'
		deallocate(NumMol,NumAtom)
		write(*,*)'OK! 8'
		deallocate(typeMol,typeAtom)
		write(*,*)'OK! 9'
		deallocate(nhbond)
		write(*,*)'OK!10'
		deallocate(hbond)
		write(*,*)'OK! 11'

		deallocate(t0)
		write(*,*)'OK! 13'
		deallocate(t1)
		write(*,*)'OK! 14'
		deallocate(indexO)
		write(*,*)'OK! 15'
		deallocate(tj)
		write(*,*)'OK! 16'
		deallocate(nHjump)
		write(*,*)'OK! 17'
		!deallocate(CRP)
		!write(*,*)'OK! 18'
	end subroutine cleanmemory
	
	subroutine closefile()
	implicit none
		close(20)
	end subroutine closefile
	
	!!!pbc, move the vector to one box!!!
	subroutine pbc(vec,boxL)
	implicit none
	integer, parameter :: dq=kind(1.0D0)
	real(kind=dq) :: vec(3),boxL
	integer :: i
	
	do i=1,3
		vec(i)=MOD(vec(i),boxL)
		if (vec(i) < 0) then
			vec(i)=vec(i)+boxL
		else 
			vec(i)=vec(i)
		end if
	end do
	
	return
	end subroutine pbc
end module
