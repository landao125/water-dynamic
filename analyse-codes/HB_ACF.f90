 !!!pbc, move the vector to one box!!!
  subroutine pbc(vec,boxL)
    implicit none
    integer, parameter :: dp=kind(1.0D0)
    real(kind=dp) :: vec(3),boxL(3)
    integer :: i  
    do i=1,3
      !vec(i)=MOD(vec(i),boxL(i))
      if (vec(i) > boxL(i)/2.0)then
        vec(i)=vec(i)-boxL(i)
      elseif(vec(i) < -boxL(i)/2.0) then
        vec(i)=vec(i)+boxL(i)
      endif
    end do
  
    return
  end subroutine pbc

Program HB_exist_ACF
  use readxtc_read
  implicit none
  
  integer, parameter :: dp=kind(1.0D0)
  integer, parameter :: ndimention=3
  integer :: natom      !total atom number
  integer :: iatom,jatom
  integer :: iOatom,jOatom,nOatom
  integer :: npatom
  integer :: nfram      !total t0,100ps
  integer :: ifram
  integer :: nsite
  integer :: nts,its    !length of ACF 
  integer :: i_initial,n_initial   !different time rigin
  integer :: nh,ihb,ihb_0
  integer :: ndxatom    !atom numbers in ndx file 
  integer :: iOregin
  integer :: ierror,i,atomCount
  real(kind=dp) :: tstep
  real(kind=dp) :: rij(3),r1(3),r2(3)
  real(kind=dp) :: roo,roh1,roh2,thetaHOO1,thetaHOO2
  character(len=50)   :: temp 
  character(len=40)   :: ndxfile 
  integer,allocatable :: GalAtomNum(:),reginAtom(:,:),nhb(:,:),hb(:,:,:)
  integer,allocatable :: hAve1(:),h0htAve1(:),hAve2(:),h0htAve2(:)  
  real(kind=dp),   allocatable :: boxL(:)        !the length of box , in anstrong
  real(kind=dp),   allocatable :: posmd(:,:,:)   !huge array, need large memory
  real(kind=dp),   allocatable :: dis_ij(:)
  real(kind=dp),   allocatable :: ACF1(:),ACF2(:)
  
  namelist /input/ natom,npatom,nfram,nsite,&
                   ndxfile,ndxatom,n_initial,nts,tstep
  
  !!!!!!!inintial input parameters!!!!!!!!
    natom        =  51603
    npatom      = 359
    nfram        =  1000
    n_initial    =  200
    nts         =  500
    nsite       = 4
    ndxfile      =  'Gal.ndx'
    ndxatom     = 88                      !there is 88 atoms in Gal
    tstep       = 0.1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  open(unit=10,file="HB_exist_input.txt",status='old',iostat=ierror)
  if (ierror /= 0) then
    write(*,*) 'Can''t open file ""HB_exist_input.txt"", please check!'
    stop
  end if
  read(unit=10,nml=input)
  close(10)
  
  if(n_initial+nts > nfram) then
    write(*,*)"Please check your input time origin(n_initial) and  calculate step (nts)!!!"
    stop
  end if 

  allocate(posmd(ndimention,natom,nfram))
  allocate(boxL(nfram))
  allocate(GalAtomNum(ndxatom),reginAtom(natom,nfram))
  allocate(dis_ij(ndxatom))
  allocate(hAve1(nts),h0htAve1(nts),ACF1(0:nts))
  allocate(hAve2(nts),h0htAve2(nts),ACF2(0:nts))
  
  nOatom=(natom-npatom)/nsite
  allocate(nhb(nOatom,nfram),hb(10,nOatom,nfram))
  
    
    !read position of every atom in all frame
    !read box length of all fram
    boxL=0.0
    posmd=0.0
    call read_xtcp(ndimention,natom,nfram,posmd,boxL)
    posmd=posmd*10.0    !change unit to angstrom
    boxL=boxL*10.0
    
    open(unit=12,file=trim(adjustl(ndxfile)),status="old",iostat=ierror)
    if (ierror /= 0) then
      write(*,*) 'Can''t open ndxfile, please check!'
      stop
    end if
    read(12,*) 
    read(12,*) (GalAtomNum(i),i=1,ndxatom)
    close(12)
  
  
   ! any |O-Gal| < 3.5A, O atom belong to regin1  
   ! all |O-Gal| >10.0A and any |O-Gal| < 12.0A ,O atom belong to regin2  
    reginAtom=0
    dis_ij=0.0
    loop1:  do ifram=1,nfram
      loop2:  do iatom=npatom+1,natom,nsite
               atomCount=0
        loop3: do jatom=1,ndxatom
                rij=posmd(:,jatom,ifram)-posmd(:,iatom,ifram)
                call pbc(rij,boxL(ifram))
                dis_ij(jatom)=sqrt(sum(rij*rij))
                if (dis_ij(jatom)<3.5)then   
                  reginAtom(iatom,ifram)=1  
                  exit loop3
                elseif(dis_ij(jatom)>10.0)then  
                  atomCount=atomCount+1         
                  if(atomCount==ndxatom .and. ANY(dis_ij<12.0))reginAtom(iatom,ifram)=2
                end if
               end do loop3
        
              end do loop2      
    
            end do loop1
  
 
  
    nhb=0
    hb=0
    do ifram=1,nfram
      do iOatom=1,nOatom
        iatom = npatom+1+nsite*(iOatom-1)
        iOregin = reginAtom(iatom,ifram)
        if (reginAtom(iatom,ifram)/=0)then
          do jOatom=1,nOatom
            jatom = npatom+1+nsite*(jOatom-1)          
            if (jOatom/=iOatom .and. reginAtom(jatom,ifram)==iOregin)then
              rij=posmd(:,jatom,ifram)-posmd(:,iatom,ifram)
              call pbc(rij,boxL(ifram))
              roo=sqrt(sum(rij*rij))
              if(roo<0.005) write(*,*) "roo will be wrong."
              if(roo<3.5)then
                r1=posmd(:,iatom+1,ifram)-posmd(:,iatom,ifram)
                call pbc(r1,boxL(ifram))
                roh1=sqrt(sum(r1*r1))
                if(roh1<0.005) write(*,*) "roh1 will be wrong."
                thetaHOO1=57.29577951308232088d0*acos(dot_product(rij,r1)/(roo*roh1)) 
                r2=posmd(:,iatom+2,ifram)-posmd(:,iatom,ifram)
                call pbc(r2,boxL(ifram))
                roh2=sqrt(sum(r2*r2))
                if(roh2<0.005) write(*,*) "roh2 will be wrong."
                thetaHOO2=57.29577951308232088d0*acos(dot_product(rij,r2)/(roo*roh2)) 
                if (thetaHOO1< 30.0 .or. thetaHOO2< 30.0 )then
                  nhb(iOatom,ifram)=nhb(iOatom,ifram)+1
                  hb(nhb(iOatom,ifram),iOatom,ifram)=jatom
                end if
              end if
            end if
          end do
        end if
      end do
    end do
  
   
   deallocate(posmd)
   deallocate(boxL)
   deallocate(GalAtomNum)
   deallocate(dis_ij)
   
   write(*,*)"11 ok"
  
  hAve1=0
  h0htAve1=0
  hAve2=0
  h0htAve2=0
  do its=1,nts   !resolve hAve(its) and h0htAve(its)
    do i_initial=1,n_initial
      do iOatom=1,nOatom
        iatom = npatom+1+nsite*(iOatom-1)
        if(reginAtom(iatom,i_initial)==1) then
          hAve1(its)=hAve1(its)+nhb(iOatom,i_initial+its)
          do ihb=1,nhb(iOatom,i_initial+its)
            do ihb_0=1,nhb(iOatom,i_initial)
              if(hb(ihb,iOatom,i_initial+its)==hb(ihb_0,iOatom,i_initial)) then
                h0htAve1(its)=h0htAve1(its) + 1
              endif
            enddo
          enddo              
        elseif(reginAtom(iatom,i_initial)==2) then
          hAve2(its)=hAve2(its)+nhb(iOatom,i_initial+its)
          do ihb=1,nhb(iOatom,i_initial+its)
            do ihb_0=1,nhb(iOatom,i_initial)
              if(hb(ihb,iOatom,i_initial+its)==hb(ihb_0,iOatom,i_initial)) then
                h0htAve2(its)=h0htAve2(its) + 1
              endif
            enddo
          enddo              
        endif    
        !do jatom=iatom+1,oAtom
        !  hAve(its)=hAve(its)+hb(jatom,iatom,i_initial+its)
        !  h0htAve(its)=h0htAve(its)+hb(jatom,iatom,i_initial+its)*hb(jatom,iatom,i_initial)
        !end do
      end do
    end do
  end do
  ACF1=real(h0htAve1)/real(hAve1)
  ACF2=real(h0htAve2)/real(hAve2)
  write(*,*)"12 ok"
  
  
  open (unit=101, file='HB_exist_ACF.out',status='unknown')
  ACF1(0)=1.0
  ACF2(0)=1.0
  do its=0,nts
    write(101,'(3(f12.6,2x))') real(its)*tstep,ACF1(its),ACF2(its)
  end do
  write(*,*)"13 ok"


  deallocate(reginAtom)
  deallocate(nhb,hb)
  deallocate(hAve1,h0htAve1,ACF1)
  deallocate(hAve2,h0htAve2,ACF2)
  
  
  close(101)

End program HB_exist_ACF    



