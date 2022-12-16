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
	call findHbond()
	call findSS()
	call geometric_quantities()
	if(Lzq_SSP_CF)then
		call Cal_SSP_CF_ZQ()
	else
		call Cal_SSP_CF()
	endif
	call Cal_SSP_CF_fram
	call cleanmemory()
	!call closefile()
	write(*,*) "Calculate is finished!"
stop
End program jump_analy


