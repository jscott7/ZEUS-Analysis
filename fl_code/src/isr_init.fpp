c     =============================
      subroutine isr_init(choice,num)
c     =============================
c
c --- Start of processing initialisations
c
      implicit none

#include "local.inc"

c ---------------------------------------------------------
c Define the bin for the Fl analysis
c ---------------------------------------------------------
c      real ybins(2)
c 0.11 - 0.21
c      data ybins /0.1, 0.23/

c      real q2bins(2)
c      data q2bins /1.0, 30.0/

c      real xmin, xmax, ymean

      logical first

      data first /.true./

      save first

      integer num, num2, num3

      character*1 a
      character*1 b
      character*1 c

      character*8 fname
      integer choice, istat, i, j, k, l, m, n
      real hmem
      common/pawc/hmem(2000000)

      if (first.eqv..true.) then
         call hlimit(2000000)
         first = .false.
      end if
c
c --- Open output hbook file
c
      write(*,*) 'Loop:', num

      if (num.le.9) then
         write (a,'(i1)') num
         fname = 'isr1_'//a
      else if (num.gt.9.and.num.le.18) then 
         num2 = num-9
         write (b,'(i1)') num2
         fname = 'isr2_'//b
      else if (num.gt.18.and.num.le.26) then
	  num3 = num-18
          write (c,'(i1)') num3
          fname = 'isr3_'//c
      end if

      if (choice.eq.1) then
         call hropen(50,'histos','israll96.hbook','N',1024,istat)
      elseif (choice.eq.2) then
         call hropen(50,'histos','israll97.hbook','N',1024,istat)
      elseif (choice.eq.3) then
         call hropen(50,'histos',fname,'N',1024,istat)
      endif

      if (istat.ne.0) then 
	 write (*,*) 'INITIALISE: Failed to open output file'
         write (999,*) 'INITIALISE: Failed to open output file'
         stop
      endif

c -----------------------------------
c Histogram bookings....
c -----------------------------------
c Data
c -----------------------------------

      call hbook1(100,'DATA E-Pz (Cal all cuts)',100,0.,100.,0.)
      call hbook1(101,'DATA E-Pz (Tot all cuts)',100,20.,120.,0.) 
      call hbook1(103,'DATA E-Pz (Zufos had)',50,0.,40.,0.)  
      call hbook1(104,'DATA E-Pz (Cal)',100,0.,100.,0.)
      call hbook1(105,'DATA E-Pz (total)',100,20.,120.,0.)
      call hbook1(108,'DATA Uncorrected El_energy',30,3.,27.,0.)
      call hbook1(109,'DATA Uncorrected Q2_electron',40,-1.,3.,0.)
      call hbook1(111,'DATA y el',40,-2.,0.,0.)
      call hbook1(112,'DATA y el (corrected for ISR)',40,-2.,0.,0.)
      call hbook1(113,'DATA y sig',15,-0.9,-0.1,0)
      call hbook1(114,'DATA y jb(Zufos)',60,-3.,0.,0)
      call hbook1(116,'DATA x_el',30,-5.,0.,0.)
      call hbook2(117,'DATA Q2vX',30,-5.,-0.,30,-1.,3.,0.)
      call hbook1(118,'DATA z-vtx',100,-120.,120.,0.) 
      call hbook1(120,'DATA electron energy',30,3.,27.,0.)
      call hbook1(121,'DATA Q2_electron',50,-1.2,4.,0.)
      call hbook1(122,'DATA E lumig',30,3.,25.,0.) 
      call hbook1(129,'Gamma Had',60,0.,180.,0.)
      call hbook1(131,'DATA eta_max',56,-7.,7.,0.)
      call hbook1(132,'DATA deltaisr',100,-5.,5.,0.)
      call hbook2(150,'x v y',100,-80.,80.,100,-80.,80.,0.)  
      call hbook1(151,'xpos',100,-80.,80.,0.)
      call hbook1(152,'ypos',100,-80.,80.,0.)
      call hbook1(153,'best theta',100,0.,3.2,0.)
      call hbook1(160,'ysig x z',60,-3.,-0.,0.)
      call hbook1(161,'y (sig*z)',15,0.,0.5,0.)
      call hbook1(162,'y (sig)',15,0.05,1.,0.)
      call hbook1(163,'y (sig) unweighted',15,0.05,1.,0.)
      call hbook1(171,'DATA y el',40,-2.,0.,0.)
      call hbook1(172,'DATA y el (corrected for ISR)',40,-2.,0.,0.)
      call hbook1(173,'DATA y sig',60,-3.,0.,0)
      call hbook1(174,'DATA y jb(Zufos)',60,-3.,0.,0)
      call hbook1(176,'Had E-Pz',60,0.,20.,0.)

c -----------------------------------
c Background
c -----------------------------------
      call hbook1(200,'Bak E-Pz (Cal all cuts)',100,0.,100.,0.)
      call hbook1(201,'Bak E-Pz (Tot all cuts)',100,20.,120.,0.) 
      call hbook1(203,'Bak E-Pz (Zufos had)',50,0.,40.,0.)	           
      call hbook1(204,'Bak E-Pz (Cal)',100,0.,100.,0.)
      call hbook1(205,'Bak E-Pz (total)',100,20.,120.,0.)
      call hbook1(208,'Bak Uncorrected El_energy',30,3.,27.,0.)
      call hbook1(209,'Bak Uncorrected Q2_electron',40,-1.,3.,0.)
      call hbook1(212,'Bak y_el',40,-2.,0.,0.)
      call hbook1(213,'Bak y_sig',15,-0.9,-0.1,0)
      call hbook1(214,'Bak y_jb(Zufos)',60,-3.,0.,0)
      call hbook1(216,'Bak x_el',30,-5.,0.,0.)
      call hbook2(217,'Bak Q2vX',30,-4.,-0.,30,0.,3.,0.)
      call hbook1(218,'Bak z-vtx',100,-120.,120.,0.)     
      call hbook1(220,'Bgd electron_energy(wt)',30,3.,27.,0.)
      call hbook1(221,'Bgd Q2_electron(wt)',50,-1.2,4.,0.)
      call hbook1(222,'Bgd E lumig',30,3.,25.,0.)
      call hbook1(223,'Bak ENE_LG',30,3.,25.,0.)
      call hbook1(224,'Bak LUMIGTMP',30,3.,25.,0.)
      call hbook1(229,'Bak Gamma Had',60,0.,180.,0.)
      call hbook1(231,'Bak eta_max',56,-7.,7.,0.)
      call hbook1(232,'Bak delta_isr',100,-5.,5.,0.)
      call hbook1(253,'Bak best theta',100,0.,3.2,0.)
      call hbook1(260,'ysig x z',60,-3.,-0.,0.)
      call hbook1(261,'y (sig*z)',15,0.,0.5,0.)
      call hbook1(262,'y (sig)',15,0.05,1.,0.)
      call hbook1(263,'y (sig) unweighted',15,0.05,1.,0.)
      call hbook1(272,'Bak y_el',40,-2.,0.,0.)
      call hbook1(273,'Bak y_sig',60,-3.,0.,0)
      call hbook1(274,'Bak y_jb(Zufos)',60,-3.,0.,0)
      call hbook1(276,'Bak Had E-Pz',60,0.,20.,0.)

c -----------------------------------
c MonteCarlo
c -----------------------------------
      call hbook1(300,'MC E-Pz (Cal all cuts)',100,0.,100.,0.)
      call hbook1(301,'MC E-Pz (Tot all cuts)',100,20.,120.,0.) 
      call hbook1(303,'MC E-Pz (Zufos had)',50,0.,40.,0.)	           
      call hbook1(304,'MC E-Pz (Cal)',100,0.,100.,0.)
      call hbook1(305,'MC E-Pz (total)',100,20.,120.,0.)
      call hbook1(308,'MC Uncorrected El_energy',30,3.,27.,0.)
      call hbook1(309,'MC Uncorrected Q2_electron',40,-1.,3.,0.)
      call hbook1(310,'MC y_sig fl',15,-0.9,-0.1,0)
      call hbook1(311,'MC y_sig fl',15,-0.9,-0.1,0)
      call hbook1(312,'MC y_el',40,-2.,0.,0.)
      call hbook1(313,'MC y_sig',15,-0.9,-0.1,0)
      call hbook1(314,'MC y_jb(Zufos)',60,-3.,0.,0)
      call hbook1(315,'MC y_sig fl',15,-0.9,-0.1,0)
      call hbook1(316,'MC x_el',30,-5.,0.,0.)
      call hbook2(317,'MC Q2vX',30,-4.,-0.,30,0.,3.,0.)
      call hbook1(318,'MC z-vtx',100,-120.,120.,0.)       
      call hbook1(320,'MC electron_energy(wt)',30,3.,27.,0.)
      call hbook1(321,'MC Q2_electron(wt)',50,-1.2,4.,0.)
      call hbook1(322,'MC E lumig(wt)',30,3.,25.,0.) 
      call hbook1(323,'MC y_sig(wt)',60,-3.,0.,0.) 
      call hbook1(329,'MC Gamma Had',60,0.,180.,0.)
      call hbook1(331,'MC eta_max',56,-7.,7.,0.)
      call hbook1(332,'MC delta_isr',100,-5.,5.,0.)
      call hbook1(353,'MC best theta',100,0.,3.2,0.)   
      call hbook1(360,'ysig x z',20,-2.,0.,0.)
      call hbook1(361,'y (sig*z)',15,0.,0.5,0.)
      call hbook1(362,'y (sig)',15,0.05,1.,0.)
      call hbook1(363,'y (sig) unweighted',15,0.05,1.,0.)
      call hbook1(364,'y (sig) R=1.4',15,0.05,1.,0.)
      call hbook1(372,'MC y_el',40,-2.,0.,0.)
      call hbook1(373,'MC y_sig',60,-3.,0.,0)
      call hbook1(374,'MC y_jb(Zufos)',60,-3.,0.,0)
      call hbook1(376,'MC Had E-Pz',60,0.,20.,0.)
      call hbook1(377,'MC electron_energy',30,3.,27.,0.)
      call hbook1(378,'MC E lumig',30,3.,25.,0.) 

c --------------------------------------
c Montecarlo True
c --------------------------------------
      call hbook1(401,'MC True electron_en',30,3.,27.,0.)
      call hbook1(402,'MC True Elumig',30,3.,25.,0.)
      call hbook1(403,'MC True Elumie',30,3.,25.,0.)
      call hbook1(404,'MC True Q2',40,-1.,3.,0.)
      call hbook1(405,'MC True X',30,-5.,0.,0.)
      call hbook1(406,'MC True Y',20,-2.,0.,0.)
      call hbook1(407,'MC True Y2',20,-2.,0.,0.)
      call hbook1(408,'MC RW True X',30,-5.,0.,0.)
      call hbook1(409,'MC RW True Y',20,-2.,0.,0.)
      call hbook1(410,'MC RW True Zvtx',100,-120.,120.,0.)

c y resolution plots
      call hbook1(430,'temp',40,-1.,1.,0.)
      call hbook1(431,'temp',40,-1.,1.,0.)
      call hbook1(432,'temp',40,-1.,1.,0.)
      call hbook1(433,'temp',40,-1.,1.,0.)
      call hbook1(434,'temp',40,-1.,1.,0.)
      call hbook1(435,'temp',40,-1.,1.,0.)
      call hbook1(436,'temp',40,-1.,1.,0.)
      call hbook1(437,'temp',40,-1.,1.,0.)
      call hbook1(438,'temp',40,-1.,1.,0.)
      call hbook1(439,'temp',40,-1.,1.,0.)

c Lumi resolution plots
      call hbook1(461,'temp',60,-5.,5.,0.)
      call hbook1(462,'temp',60,-5.,5.,0.)

c KP plots
      call hbook1(500,' ',30,20.,30.,0.)
      call hbook1(501,' ',30,20.,30.,0.)
      call hbook1(502,' ',30,20.,30.,0.)
      call hbook1(503,' ',30,20.,30.,0.)
      call hbook1(504,' ',30,20.,30.,0.)
      call hbook1(505,' ',30,20.,30.,0.)
      call hbook1(506,' ',30,20.,30.,0.)
      call hbook1(507,' ',30,20.,30.,0.)
      call hbook1(508,' ',30,20.,30.,0.)
      call hbook1(509,' ',30,20.,30.,0.)
      call hbook1(510,' ',30,20.,30.,0.)
      call hbook1(511,' ',30,20.,30.,0.)

c -----------------------------------
c Initialise lumis for weighting
c -----------------------------------
      lumda = 10.524
      lumda97 = 25.39882

c --------------------------------
c Non Diffractive > 0.1
      lumcnd1 = 0.111
      lumcnd197 = 0.212
c Diffractive > 0.1
      lumcd1 = 0.052
      lumcd197 = 0.0735
c --------------------------------
c Non Diffractive > 0.2
      lumcnd2 = 0.826
      lumcnd297 = 0.431
c Diffractive > 0.2
      lumcd2 = 0.65
      lumcd297 = 0.2099
c --------------------------------
c Non Diffractive > 0.5
      lumcnd3 = 4.7894
      lumcnd397 = 3.99
c Diffractive > 0.5
      lumcd3 = 0.741
      lumcd397 = 2.883
c --------------------------------
c Non Diffractive > 0.75 FL ONLY
      lumcnd3197 = 23.046
      lumcd3197 = 0.0
c --------------------------------
c Non Diffractive > 2.0
      lumcnd4 = 12.038
      lumcnd497 = 6.295
c Diffractive > 2.0
      lumcd4 = 0.0
      lumcd497 = 0.0
c --------------------------------
c If no Diff files then set to 0. 
c      fracdif = 0.15
c	fracdif97 = 0.15

      nmcfiles = 17
      nnodiffiles = 13

      nmcfiles97 = 35
      nnodiffiles97 = 29

c Without Q^2 > 0.75 MC
c      nmcfiles97 = 16
c	nnodiffiles97 = 10

      q2mc1 = 0.1
      q2mc2 = 0.2
      q2mc3 = 0.5

      q2mc31 = 0.75
      q2mc4 = 2.0
      
c -----------------------------------
c Initialise Bin
c Only 1 bin for Fl !!!!
c -----------------------------------
c      bindat = 0.
c      binbgd = 0.
c      binmc = 0.
c	bakraw = 0.
c	datraw = 0.
c	mcraw = 0.

c      q2min = q2bins(1)
c      q2max = q2bins(2)
c      ymin = ybins(1)
c      ymax = ybins(2)
c      xmin = q2min / (90200. * ymin)
c      xmax = q2max / (90200. * ymax)

c      write(*,*) '********************************'
c      write(*,*) '*** Q^2 and y bin boundaries ***'
c      write(*,*) '********************************'
c      write(999,*) '********************************'
c      write(999,*) '*** Q^2 and y bin boundaries ***'
c      write(999,*) '********************************'
 

c	ymean = (10.**((log10(ymin) + log10(ymax))/2.))
c      q2mean = (10.**((log10(q2min) + log10(q2max))/2.))
c	xmean = q2mean/(ymean*90200.)  

c	logxmean = log10(xmean)

c      write(*,*) xmean, q2mean    
c      write(*,*) q2min, q2max, ymin, ymax
c      write(999,*) xmean, q2mean    
c      write(999,*) q2min, q2max, ymin, ymax

c ----------------------------------
c PDF Init 
c ----------------------------------
      call PDFINIT

      return 
      end
