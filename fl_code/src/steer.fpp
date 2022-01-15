      subroutine steer(choice,num)
c Choice = 1 for 1996 and 2 for 1997
	implicit none

#include "local.inc"
	integer null,choice,num,num2,num3
	character*1 a,b,c
	character*9 fname, fname2

	if (num.le.9) then
        write (a,'(i1)') num
        fname = 'st961_'//a
	  fname2 = 'st971_'//a
	else if (num.gt.9.and.num.le.18) then 
	  num2 = num-9
        write (b,'(i1)') num2
        fname = 'st962_'//b
        fname2 = 'st972_'//b
	else if (num.gt.18.and.num.le.26) then
	  num3 = num-18
        write (c,'(i1)') num3
        fname = 'st963_'//c
        fname2 = 'st973_'//c
	end if	

	if (choice.eq.1) then
      open(unit=20,file=fname,status='old')
	read(20,*,err=90) null
	read(20,*,err=90) lumiscale
	read(20,*,err=90) lumires
	read(20,*,err=90) enelow
	read(20,*,err=90) elglow
	read(20,*,err=90) empztothigh
	read(20,*,err=90) empztotlow
	read(20,*,err=90) empzcalcard
	read(20,*,err=90) zvertex
	read(20,*,err=90) bgdnormal
	read(20,*,err=90) boxcutcard1
	read(20,*,err=90) tag44mveto
	read(20,*,err=90) hadscale
	read(20,*,err=90) elenscale
	read(20,*,err=90) offlineuse
	read(20,*,err=90) bkspuse
	read(20,*,err=90) ybinhigh
	read(20,*,err=90) ybinlow
	read(20,*,err=90) srtdpressc
	read(20,*,err=90) diffractiv


90	continue
	close(20)
	end if

      if (choice.eq.2) then
	open(unit=20,file=fname2,status='old')
	read(20,*,err=80) null
	read(20,*,err=80) lumiscale
	read(20,*,err=80) lumires
	read(20,*,err=80) enelow
	read(20,*,err=80) elglow
	read(20,*,err=80) empztothigh
	read(20,*,err=80) empztotlow
	read(20,*,err=80) empzcalcard
	read(20,*,err=80) zvertex
	read(20,*,err=80) bgdnormal
	read(20,*,err=80) boxcutcard1
	read(20,*,err=80) tag44mveto
	read(20,*,err=80) hadscale
	read(20,*,err=80) elenscale
	read(20,*,err=80) offlineuse
	read(20,*,err=80) bkspuse
	read(20,*,err=80) ybinhigh
	read(20,*,err=80) ybinlow
	read(20,*,err=80) srtdpressc
	read(20,*,err=80) diffractiv

80	continue
	close(20)
	end if
	
	write(*,*) '------------------------------------------'
	write(*,*) '    Steering parameters for F2ISR         '
	write(*,*) ' '
	write(*,*) ' Lumi energy scale :',lumiscale
	write(*,*) ' Lumi res :',lumires
	write(*,*) ' Electron en lower cut:',enelow
	write(*,*) ' Lumi photon energy lower cut:',elglow
	write(*,*) ' empztot upper cut :',empztothigh
	write(*,*) ' empztot lower cut :',empztotlow
	write(*,*) ' empz cal cut :',empzcalcard
	write(*,*) ' z vertex cut :',zvertex
	write(*,*) ' bgd normalisation e-pz:',bgdnormal
	write(*,*) ' Do boxcut +- 0.5:',boxcutcard1
	write(*,*) ' Do 44m Tagger veto:',tag44mveto
	write(*,*) ' Had energy scale :',hadscale
	write(*,*) ' electron energy scale :',elenscale
	write(*,*) ' offline trigger:',offlineuse
	write(*,*) ' backsplash :',bkspuse
	write(*,*) ' y bin high:',ybinhigh
	write(*,*) ' y bin low:',ybinlow
	write(*,*) ' srtd pres scale:',srtdpressc
	write(*,*) ' diffractive frac:',diffractiv
	write(*,*) ' '
	write(*,*) '-----------------------------------------'

 	write(999,*) '------------------------------------------'
	write(999,*) '    Steering parameters for F2ISR         '
	write(999,*) ' '
	write(999,*) ' Lumi energy scale :',lumiscale
	write(999,*) ' Lumi res :',lumires
	write(999,*) ' Electron en lower cut:',enelow
	write(999,*) ' Lumi photon energy lower cut:',elglow
	write(999,*) ' empztot upper cut :',empztothigh
	write(999,*) ' empztot lower cut :',empztotlow
	write(999,*) ' empz cal cut :',empzcalcard
	write(999,*) ' z vertex cut :',zvertex
	write(999,*) ' bgd normalisation e-pz:',bgdnormal
	write(999,*) ' Do boxcut +- 0.5:',boxcutcard1
	write(999,*) ' Do 44m Tagger veto:',tag44mveto
	write(999,*) ' Had energy scale :',hadscale
	write(999,*) ' electron energy scale :',elenscale
	write(999,*) ' offline trigger:',offlineuse
	write(999,*) ' backsplash :',bkspuse
	write(999,*) ' y bin high:',ybinhigh
	write(999,*) ' y bin low:',ybinlow
	write(999,*) ' srtd pres scale:',srtdpressc
	write(999,*) ' diffractive frac:',diffractiv
	write(999,*) ' '
	write(999,*) '-----------------------------------------'

	return
	end
