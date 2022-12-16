!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!change t0 sample
!change boxL to array													       	
!Last update 2018-3-15															                
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module parameters
implicit none

	integer:: natom  !total atom number
	integer :: iatom,jatom
	integer, parameter :: dq=kind(1.0D0)
	integer::inifram    !initial fram
	integer :: nfram   !total t0,100ps
	integer :: ifram
	integer:: nwater      !total number of water
	integer :: iwater
	real(kind=dq),allocatable :: boxL(:)   !the length of box , in anstrong
	character(len=50) :: temp
	integer,allocatable  :: NumMol(:),NumAtom(:)
	character(len=5),allocatable :: typeMol(:),typeAtom(:)  !atom type

	integer ii
	integer, parameter :: ndimention=3
	integer:: idimention
	real(kind=dq), allocatable :: posmd(:,:,:)   !huge array, need large memory

	!H-bond definition
	real(kind=dq):: R_OO          !in anstrong
	real(kind=dq):: theta_HOO       !in degree
	real(kind=dq):: theta
	real(kind=dq), parameter :: pi=3.1415926
	integer :: iHjump   !count jump times of every H atom
	real(kind=dq) :: r1(3),r2(3),roo,roh(2),thetaHOO(2)
	integer, allocatable :: nhbond(:,:)  !每一帧中每个H原子形成氢键的数目
	integer, allocatable :: hbond(:,:,:) !每一帧中与每个氢原子形成nhbond个氢键的氧原子的序号
	integer, allocatable :: nhbondO(:,:) !每一帧中每个O原子形成氢键的数目
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
	real(kind=dq):: corr        !
	real(kind=dq),allocatable :: CRP(:)
	real:: tstep   !time step for saving
	character(len=40)::Package
	character(len=40)::gromacsfile
	character(len=40)::tinkerfile
	character(len=8):: inputfile='InputCar'
	Logical::Lzq_SSP_CF
	integer :: ts     !计算几何结构参量的时候取的时间范围
	
	namelist /input/ natom,inifram,nfram,nwater,R_OO,theta_HOO,nt,tstep,Package,gromacsfile,tinkerfile,Lzq_SSP_CF,ts
	contains
	
	subroutine readinput()
	implicit none
		
		!!!!!!!inintial input parameters!!!!!!!!
		natom				=	10
		inifram     = 2
		nfram				=	100
		nwater			=	30
		R_OO				=	3.4
		theta_HOO		=	20.0
		nt 					=	5000
		tstep				=	0.01
		Package 		=	'tinker'
		gromacsfile	=	'water.gro'
		tinkerfile	=	'spce-nve.arc'
		Lzq_SSP_CF	= .TRUE.
		ts          = 199
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		open(unit=12,file=inputfile,status='old')
		read(12,nml=input)
		close(12)
		nf=int(nfram/2)
		allocate(boxL(nfram))
		allocate(NumMol(natom),NumAtom(natom))
		allocate(typeMol(natom),typeAtom(natom))
		allocate(posmd(ndimention,natom,nfram))
		allocate(nhbond(2*nwater,nfram))      !大型数组
		allocate(hbond(10,2*nwater,nfram))    !设定一个氢原子最多能同时形成10个氢键，大型数组
		allocate(nhbondO(2*nwater,nfram))
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

		open (unit=11, file=trim(adjustl(gromacsfile)),status='old')
		read(11,*)temp
		read(11,*)temp  
		Do iatom=1, natom
			read (11,"(I5,2A5,I5)")NumMol(iatom),typeMol(iatom),typeAtom(iatom),NumAtom(iatom)
		End do	
		rewind(11)
		close(11)
		
		!read position of every atom in all frame
		!read box length of all fram
		boxL=0.0
		posmd=0.0
		call read_xtcp(ndimention,natom,nfram,posmd,boxL)
		posmd=posmd*10.0    !change unit to anstrong
		boxL=boxL*10.0		
		
	end subroutine readgro
	
	subroutine readtinker()
	implicit none
		integer::atomtype
		!allocate(typeAtom(natom))
		open(unit=11,file=trim(adjustl(tinkerfile)),status='old')
		do ifram=1,inifram-1
			do iatom=1,natom+2
				read(11,*)temp
			end do
		end do
		do ifram=1,nfram
			read(11,*) temp
			read(11,"(1X,F12.6)") boxL(ifram)    !here boxL is an array
			do iatom=1,natom
				read(11,"(11X,3F12.6,I6)") (posmd(idimention,iatom,ifram),idimention=1,ndimention),atomtype
				if(atomtype==1) then
					typeAtom(iatom)='   OW'
				elseif(atomtype==2) then
					typeAtom(iatom)='   HW'
				endif				
			enddo
			!read(11,"(F10.5)") boxL
		enddo
		close(11)
	end subroutine readtinker
	
	subroutine findHbond()
	implicit none
	!theta=theta_HOO*pi/180.0   !change the unit to radian	
	nhbond=0
	hbond=0
	nhbondO=0
	Do ifram=1,nfram
		iHatom=0
		DO iatom=1,natom
			If (typeAtom(iatom)=='   OW') then
				Do ii=1,2				!every water have two H atom
					iHatom=iHatom+1    !记录氢原子序号
					do jatom=1,natom
						if ((jatom/=iatom) .AND. typeAtom(jatom)=='   OW') then
							r1(:)=posmd(:,jatom,ifram)-posmd(:,iatom,ifram)    !O*Oa vector
							call pbc(r1,boxL(ifram))
							roo=sqrt(sum(r1*r1))                !O*Oa scalar
							if(roo < R_OO) then								         		
								r2(:)=posmd(:,iatom+ii,ifram)-posmd(:,iatom,ifram)    !O*H vector
								call pbc(r2,boxL(ifram))
								roh(ii)=sqrt(sum(r2*r2))               !O*H scalar								
								thetaHOO(ii)=acos(dot_product(r1,r2)/(roo*roh(ii)))   ! theta H*O*Oa
								thetaHOO(ii)=57.29577951308232088d0*thetaHOO(ii)  !弧度转成角度
								if(thetaHOO(ii) < theta_HOO) then
									nhbond(iHatom,ifram)=nhbond(iHatom,ifram)+1    !第ifram帧第iHatom个氢原子形成nhbond个氢键
									hbond(nhbond(iHatom,ifram),iHatom,ifram)=jatom   !与第ifram帧第iHatom个氢原子形成nhbond个氢键的氧原子序号
									nhbondO(jatom,ifram)=nhbondO(jatom,ifram)+1     !注意数组下标和iHatom表示的意义不同
								end if											                !第ifram帧序号为jatom的氧原子形成的氢键数目
							end if						
						end if					
					end do  !jatom
				End do  !ii
			End if
		End do  !iatom	
	End do   !ifram
	!deallocate(posmd)
	end subroutine findHbond
	

	
	subroutine findSS()
	implicit none
	Logical::Lflag=.False. !用于判断MD是否经过了一次非稳态
	t0=0
	t1=0
	indexO=0
	nHatom = 2*nwater
	Do iHatom=1,nHatom 
		nSS(iHatom)=0           !记录稳态数目		
		do ifram=1,nfram
			Lflag= (Lflag .or. (nhbond(iHatom,ifram)/=1)) !确保经过一次非稳态以后，进入循环
			if(Lflag) then
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
			endif
		end do !ifram
	End do  !iHatom
	end subroutine findSS
	
	subroutine geometric_quantities()
	implicit none
		integer,parameter :: tsize=1000
		integer :: njump
		integer :: mid,Ostar,Hstar,Oa,Ob
		integer :: ti,tf,j0,angtol
		integer,allocatable :: icount(:),popu(:)
		real(kind=dq) :: phi,psai,angjump,theta1
		real(kind=dq),allocatable :: rOsOa(:),rOsOb(:),rOaOb(:)
		real(kind=dq),allocatable :: angtheta(:),angphi(:),angpsai(:),ang(:)
		real(kind=dq) :: eOsOa(3),eOsOb(3),eOaOb(3),ebisec(3),epn(3),eOsHs(3),ep(3),epj(3),en(3),OsHs(3)
		integer,allocatable :: nrOs(:),nrOa(:),nrOb(:)
		
		allocate(icount(-tsize:tsize))
		allocate(rOsOa(-tsize:tsize))
		allocate(rOsOb(-tsize:tsize))
		allocate(rOaOb(-tsize:tsize))
		allocate(angtheta(-tsize:tsize))
		allocate(angphi(-tsize:tsize))
		allocate(angpsai(-tsize:tsize))
		allocate(popu(500))        !think twice
		allocate(nrOs(-tsize:tsize))
		allocate(nrOa(-tsize:tsize))
		allocate(nrOb(-tsize:tsize))
		allocate(ang(-tsize:tsize))
		njump=0
		icount=0
		rOsOa=0.0
		rOsOb=0.0
		rOaOb=0.0
		angtheta=0.0
		angphi=0.0
		angpsai=0.0
		popu=0
		j0=0
		angjump=0.0
		nrOs=0
		nrOa=0
		nrOb=0
		
		open(unit=21,file='geometric_quantities.out',status='replace')
		open(unit=22,file='jump-angle.out',status='replace')
		open(unit=23,file='h-bonds.out',status='replace')
		Do iHatom=1,nHatom
			if (nSS(iHatom) == 0 .OR. nSS(iHatom) == 1) then 
				nHjump(iHatom)=0
			else
				nHjump(iHatom)=nSS(iHatom)-1
			end if
		End do
		
		Do iHatom=1,nHatom
			if(nHjump(iHatom) > 0) then
				Ostar=3*int((iHatom-1)/2)+1
				Hstar=ceiling(real(iHatom)/2.0)+iHatom   !向上取整
				Do iHjump=1,nHjump(iHatom)
					Oa=indexO(iHjump,iHatom)
					Ob=indexO(iHjump+1,iHatom)
					mid=int((t1(iHjump,iHatom)+t0(iHjump+1,iHatom))/2)
					ti=t0(iHjump,iHatom)			
					if(mid-ti > ts) ti=mid-ts
					tf=t1(iHjump+1,iHatom)
					if(tf-mid > ts) tf=mid+ts
				
					do it=ti,tf
						icount(it-mid)=icount(it-mid)+1
						r1(:)=posmd(:,Oa,it)-posmd(:,Ostar,it)
						call pbc(r1,boxL(it))
						rOsOa(it-mid)=rOsOa(it-mid)+sqrt(sum(r1*r1))   !O*Oa 矢量的长度
						eOsOa=r1/sqrt(sum(r1*r1))   !O*Oa单位矢量
					
						r1(:)=posmd(:,Ob,it)-posmd(:,Ostar,it)
						call pbc(r1,boxL(it))
						rOsOb(it-mid)=rOsOb(it-mid)+sqrt(sum(r1*r1))   !O*Ob 矢量的长度
						eOsOb=r1/sqrt(sum(r1*r1))   !O*Ob单位矢量
					
						r1(:)=posmd(:,Ob,it)-posmd(:,Oa,it)
						call pbc(r1,boxL(it))
						rOaOb(it-mid)=rOaOb(it-mid)+sqrt(sum(r1*r1))   !OaOb 矢量的长度
						eOaOb=r1/sqrt(sum(r1*r1))   !OaOb单位矢量
					
					  !!calculate the angle between the projection of the O*H* vector 
					  !&   on the OaO*Ob plane and the bisector of the OaO*Ob angle
						r1(:)=eOsOa(:)+eOsOb(:)
						ebisec=r1/sqrt(sum(r1*r1))  !OaO*Ob的角平分线矢量的单位向量
					
						epn(1)=eOsOa(2)*eOsOb(3)-eOsOa(3)*eOsOb(2)    !O*Oa叉乘O*Ob矢量，得到平面的法向量epn矢量
						epn(2)=eOsOa(3)*eOsOb(1)-eOsOa(1)*eOsOb(3)
						epn(3)=eOsOa(1)*eOsOb(2)-eOsOa(2)*eOsOb(1)
						epn=epn/sqrt(sum(epn*epn))                 !epn单位向量
					
						OsHs(:)=posmd(:,Hstar,it)-posmd(:,Ostar,it)
						call pbc(OsHs,boxL(it))
						eOsHs=OsHs/sqrt(sum(OsHs*OsHs))                    !O*H*单位向量
					
						ep(1)=ebisec(2)*epn(3)-ebisec(3)*epn(2)        !ebisec向量和epn向量叉乘得到向量ep
						ep(2)=ebisec(3)*epn(1)-ebisec(1)*epn(3)
						ep(3)=ebisec(1)*epn(2)-ebisec(2)*epn(1)
						ep=ep/sqrt(sum(ep*ep))                    !O*H*矢量与角平分面的夹角

						theta=57.29577951308232088d0*acos(dot_product(eOsHs,ep))
						if(dot_product(eOsOa,ep) > 0) then 
							theta=theta-90.0d0
						else
							theta=90.0d0-theta
						end if
						angtheta(it-mid)=angtheta(it-mid)+theta
						
						!!!!compara two definition for calculate theta
						r1=OsHs-dot_product(OsHs,epn)*epn
						epj=r1/sqrt(sum(r1*r1))     !O*H*的投影矢量
						theta1=57.29577951308232088d0*acos(dot_product(epj,ebisec))
						en(1)=ebisec(2)*epj(3)-ebisec(3)*epj(2)        !ebisec向量和epj向量叉乘得到向量en
						en(2)=ebisec(3)*epj(1)-ebisec(1)*epj(3)
						en(3)=ebisec(1)*epj(2)-ebisec(2)*epj(1)
						if (dot_product(en,epn) < 0) theta1=-theta1
						ang(it-mid)=ang(it-mid)+theta1
						
					
					  !calculate the out-of-plane angle phi between the O*H* bond and the OaO*Ob plane
						phi=90.0d0-57.29577951308232088d0*acos(abs(dot_product(eOsHs,epn)))
						angphi(it-mid)=angphi(it-mid)+phi
					
					  !calculate the OaO*Ob angle
						psai=57.29577951308232088d0*acos(dot_product(eOsOa,eOsOb))
						angpsai(it-mid)=angpsai(it-mid)+psai
						if (it==mid) then
							njump=njump+1             !记录发生jump的总次数
							angjump=angjump+psai      !用于求平均跳跃角度
							j0=int(psai)+1
							popu(j0)=popu(j0)+1       !记录跳跃角度分布
						end if
						
						nrOs(it-mid)=nrOs(it-mid)+nhbondO(Ostar,it)   !the number of h-bond accepted by O*
						nrOa(it-mid)=nrOa(it-mid)+nhbondO(Oa,it)      !the number of h-bond accepted by Oa
						nrOb(it-mid)=nrOb(it-mid)+nhbondO(Ob,it)      !the number of h-bond accepted by Ob
						
					end do
				End do
			end if
		End do
		
		write(21,*)'geometric quantities with time:'
		write(21,"(8(A12))")'t','RO*Oa','RO*Ob','ROaOb','theta','phi','psai','theta1'
		write(23,*)'The number of h-bonds accepted by O*, Oa, Ob,respectively:'
		Do it=-ts,ts
			write(21,"(8(F10.5,2X))")real(it)*tstep,rOsOa(it)/real(icount(it)),rOsOb(it)/real(icount(it)),&
															rOaOb(it)/real(icount(it)),angtheta(it)/real(icount(it)),&
															angphi(it)/real(icount(it)),angpsai(it)/real(icount(it)),ang(it)/real(icount(it))
			write(23,"(4(F10.5,2X))")real(it)*tstep,nrOs(it)/real(icount(it)),nrOa(it)/real(icount(it)),&
															nrOb(it)/real(icount(it))
		End do
		write(22,*)'average jump angle=',angjump/real(njump)
		write(22,*)'jump angle distribution:'
		angtol=0
		Do it=1,180
			angtol=angtol+popu(it)
		End do
		Do it=1,180
			write(22,"(2(F10.5,2X))")real(it),real(popu(it))/real(angtol)     
		End do
		
		deallocate(icount)
		deallocate(rOsOa)
		deallocate(rOsOb)
		deallocate(rOaOb)
		deallocate(angtheta)
		deallocate(angphi)
		deallocate(angpsai)
		deallocate(popu)
		deallocate(nrOs)
		deallocate(nrOa)
		deallocate(nrOb)
		
		close(21)
		close(22)
		close(23)
	end subroutine geometric_quantities
	
	subroutine Cal_SSP_CF()
	implicit none
		open(unit=20,file='jump-fram.out',status='replace')
		write(20,*) 'jump correlation function for water-water:'
		!!由nSS（iHatom）确定nHjump(iHatom)
		nswitch = 0   !统计switch次数
		nHjump = 0
		
		Do iHatom=1,nHatom
			if (nSS(iHatom) == 0 .OR. nSS(iHatom) == 1) then 
				nHjump(iHatom)=0
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
			if(nHjump(iHatom)>0) then
				Do iHjump=1,nHjump(iHatom)    
					iswitch=iswitch+1
					tj(iswitch)=t0(iHjump+1,iHatom)-t0(iHjump,iHatom)					
				End do
			endif
		ENd do
			
		CRP=0.0
		Do it=0,nt
			corr=0.0
			Do iswitch=1,nswitch
				if(it <= tj(iswitch))then
					corr=corr+real(tj(iswitch))-real(it)
				end if
			End do
			CRP(it)=corr/SUM(tj)
			write(20,"(2(f10.5,2x))") real(it)*tstep,CRP(it)
		End do
		close(20)
		!deallocate(nSS)
		!deallocate(CRP)
		deallocate(tj)
	end subroutine Cal_SSP_CF
	
	subroutine Cal_SSP_CF_ZQ()
	implicit none
	integer :: tend
		open(unit=20,file='jump-fram_zq.out',status='replace')
		write(20,*) 'jump correlation function for water-water:'
		!!由nSS（iHatom）确定nHjump(iHatom)
		nHjump = 0
			Do iHatom=1,nHatom
				if (nSS(iHatom)> 1) then 
					nHjump(iHatom)=nSS(iHatom)-1
				end if
			End do

		CRP=0.0	
		Do iHatom=1,nHatom
			if(nHjump(iHatom > 1))then
				Do iHjump=2,nHjump(iHatom)-1					
					tend=t0(iHjump+1,iHatom)-t0(iHjump,iHatom)  !注意t0的含义
					if(tend >= nt)tend=nt-1
					do it=0,tend
						do ii=it,tend
							CRP(ii-it)=CRP(ii-it)+1.0
						end do
					end do
				End do
			end if
		ENd do
		
		if (CRP(0)/=0)then	
			do it=1,nt
				CRP(it)=CRP(it)/CRP(0)
			end do
		end if
		CRP(0)=1
		do it=0,nt
				write(20,"(2(f10.5,2x))") real(it)*tstep,CRP(it)	
		end do
		close(20)
	
	end subroutine Cal_SSP_CF_ZQ
	
	subroutine Cal_SSP_CF_fram()
	implicit none
		integer :: tend,Ostar,Oa
		integer,allocatable ::iCRP(:)
		real :: dot
	
		allocate(iCRP(nt+1))
	
		open(unit=21,file='CF-fram.out',status='replace')
		write(21,*) 'O..O fram correlation function for water-water:'
		!!由nSS（iHatom）确定nHjump(iHatom)
		nHjump = 0
			Do iHatom=1,nHatom
				if (nSS(iHatom)> 1) then 
					nHjump(iHatom)=nSS(iHatom)-1
				end if
			End do
		iCRP=0
		CRP=0.0	
		Do iHatom=1,nHatom     !注意，这里是第几个氢原子，不是原子序号
			Ostar=3*int((iHatom-1)/2)+1   !这里要用原子序号
			if(nHjump(iHatom > 0))then
				Do iHjump=1,nHjump(iHatom)					
					Oa=indexO(iHjump,iHatom)				
					tend=t0(iHjump+1,iHatom)-t0(iHjump,iHatom)  !注意t0的含义
					if(tend >= nt) then
						tend=t0(iHjump,iHatom)+nt
					else
						tend=t0(iHjump+1,iHatom)
					end if									
					do it=t0(iHjump,iHatom),tend
						do ii=it,tend
							r1(:)=posmd(:,Oa,ii)-posmd(:,Ostar,ii)
							call pbc(r1,boxL(ii))
							roo=sqrt(sum(r1*r1))
							r1=r1/roo
							r2(:)=posmd(:,Oa,it)-posmd(:,Ostar,it)
							call pbc(r2,boxL(ii))
							roo=sqrt(sum(r2*r2))
							r2=r2/roo
							dot=dot_product(r1,r2)
							CRP(ii-it)=CRP(ii-it)+0.50d0*(3.0d0*dot**2-1.0d0)
							iCRP(ii-it)=iCRP(ii-it)+1
						end do
					end do
				End do
			end if
		ENd do
				
		do it=1,nt
			if (iCRP(it)/=0)then	
				CRP(it)=CRP(it)/iCRP(it)
			end if
		end do
		
		CRP(0)=1
		do it=0,nt
				write(21,"(2(f10.5,2x))") real(it)*tstep,CRP(it)	
		end do
		close(21)
		!deallocate(iCRP)
	end subroutine Cal_SSP_CF_fram
	
	subroutine cleanmemory()
	implicit none
		!deallocate(boxL)
		deallocate(posmd)
		!deallocate(nSS)
		deallocate(NumMol,NumAtom)
		deallocate(typeMol,typeAtom)
		deallocate(nhbond)
		deallocate(hbond)
		deallocate(nhbondO)
		deallocate(t0)
		deallocate(t1)
		deallocate(indexO)
		!deallocate(nHjump)
		
		!deallocate(CRP)
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
		if (vec(i) > boxL/2.0)then
			vec(i)=vec(i)-boxL
		elseif(vec(i) < -boxL/2.0) then
			vec(i)=vec(i)+boxL
		endif
	end do
	
	return
	end subroutine pbc
	
end module
