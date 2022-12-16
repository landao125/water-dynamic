Program jump_analy
use parameters
implicit none	
	call readinput()
	if(package=='gromacs') then
		call readgro()
	elseif(package=='tinker') then
		call readtinker()
	else
		write(*,*) "The package is wrong!!Please check file type(gromacs or tinker)."
		stop
	endif
	write(*,*)'OK! 111'
	call findHbond()
	write(*,*)'OK! 222'
	call findSS()
	write(*,*)'OK! 333'
	call Cal_SSP_CF()
	write(*,*)'OK! 444'
	call cleanmemory()
	write(*,*)'OK! 555'
	!call closefile()
End program jump_analy


