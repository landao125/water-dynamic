! XDR Fortran Interface Example Program
! 2014 (c) James W. Barnett <jbarnet4@tulane.edu>
! https://github.com/wesbarnett/
module readxtc_read
!	use, intrinsic :: iso_c_binding, only: C_NULL_CHAR, C_PTR, c_f_pointer
!    use xtc
	implicit none
	contains
	
	subroutine read_xtc_natom(nnatom)
	use, intrinsic :: iso_c_binding, only: C_NULL_CHAR, C_PTR, c_f_pointer
    use xtc
    implicit none
		integer ,intent(out):: nnatom                !read parameter of allatoms and MD snape
		character (len=1024) :: filename
		integer :: NATOMS, STAT
		type(xdrfile), pointer :: xd
		logical :: ex

    ! Set the file name for C.
    filename = "traj.xtc"//C_NULL_CHAR

    inquire(file=trim(filename),exist=ex)

    if (ex .eqv. .false.) then
        write(0,*)
        write(0,'(a)') " Error: "//trim(filename)//" does not exist."
        write(0,*)
        stop
    end if

    STAT = read_xtc_natoms(filename,NATOMS)
	nnatom = NATOMS
	STAT = xdrfile_close(xd)
	end subroutine read_xtc_natom
	
subroutine read_xtcp(ndimention,nnatom,nfram,posmd,boxL)
	use, intrinsic :: iso_c_binding, only: C_NULL_CHAR, C_PTR, c_f_pointer
    use xtc
    implicit none
		integer, parameter :: dq=kind(1.0D0)
		integer ,intent(in):: ndimention,nnatom,nfram                 !read parameter of allatoms and MD snape
		real(kind=dq),intent(inout) :: posmd(ndimention,nnatom,nfram)
		real(kind=dq),intent(inout) :: boxL(nfram)
		character (len=1024) :: filename
		real, allocatable :: pos(:,:)
		
		integer :: NATOMS, STEP, STAT,ifram
		real :: box(3,3), prec, time, box_trans(3,3)
		type(C_PTR) :: xd_c
		type(xdrfile), pointer :: xd
		logical :: ex

    ! Set the file name for C.
    filename = "AFGP8_traj_10-30ns.xtc"//C_NULL_CHAR

    inquire(file=trim(filename),exist=ex)

    if (ex .eqv. .false.) then
        write(0,*)
        write(0,'(a)') " Error: "//trim(filename)//" does not exist."
        write(0,*)
        stop
    end if

    STAT = read_xtc_natoms(filename,NATOMS)
	
	if (NATOMS /= nnatom) then
		write(*,*) "THE number of atoms is wrong!!Please check!"
		stop
	end if
	
    allocate(pos(3,NATOMS))
	
    ! Open the file for reading. Convert C pointer to Fortran pointer.
    xd_c = xdrfile_open(filename,"r")
    call c_f_pointer(xd_c,xd)

    STAT = read_xtc(xd,NATOMS,STEP,time,box_trans,pos,prec)
	ifram = 0
    do while ( STAT == 0 .AND. ifram < nfram)
		ifram = ifram + 1
        ! C is row-major, whereas Fortran is column major. Hence the following.
        box = transpose(box_trans)
		boxL(ifram)=box(1,1)
        ! Just an example to show what was read in
        !write(*,'(a,f12.6,a,i0)') " Time (ps): ", time, "  Step: ", STEP
        !write(*,'(a,f12.6,a,i0)') " Precision: ", prec, "  No. Atoms: ", NATOMS
        !write(*,'(3f9.3)') pos
		posmd(:,:,ifram) = pos(:,:)

        ! This is the same order as found in the GRO format fyi
        !write(*,'(11f9.5)') box(1,1), box(2,2), box(3,3), &
        !                    box(1,2), box(1,3), & 
        !                    box(2,1), box(2,3), &
        !                    box(3,1), box(3,2) 

        STAT = read_xtc(xd,NATOMS,STEP,time,box_trans,pos,prec)

    end do

    STAT = xdrfile_close(xd)
    deallocate(pos)
	
	end subroutine read_xtcp
end module readxtc_read
