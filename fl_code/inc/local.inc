c -------------------------------------------------
c Local definitions for F2/ISR analysis
c -------------------------------------------------
c Kinematical variables
c -------------------------------------------------
      real empzcal, empztot, y_el, y_sig, y_elcorr, corrected_en,
     &     best_th, Q2_el, x_el, bestx, besty, bestz, electron_en, 
     &     electron_th, logyel, logxel, logq2elc, logysig, logq2el,
     &     logmcq2, logmcx, logmcy, elumie, elumig, q2corr, 
     &     zempzhad,  zy_jb, ccy_jb, zy_jbcorr, 
     &     logyjbz, best_th2, gamma,gamma2,z,logmcy2

      common /kinematics/ empzcal,empztot,y_el,y_sig,y_elcorr,gamma2, 
     & corrected_en,best_th,best_th2,Q2_el,x_el,bestx,besty,bestz,
     & electron_en,logyel,logxel,logq2elc,logysig,logmcq2,logmcx,
     & logmcy,elumie,elumig,q2corr,zempzhad,zy_jb,logq2el,
     & zy_jbcorr,logyjbz,electron_th,gamma,z,logmcy2

c -----------------------------------------------
c For reweighting
c -----------------------------------------------
      real lumcnd1, lumcd1, lumcnd2, lumcd2, lumcnd3, lumcd3,
     &     lumcnd4, lumcd4, lumda,fracdif,fracdif97,
     &     q2mc1, q2mc2, q2mc3, q2mc4,q2mc31,
     &     mcwtq2, f2weight,wt_vtx,bpcor
      real lumda97, lumcnd197, lumcnd297, lumcd297, lumcnd397, 
     &     lumcd397,lumcd197,lumcnd497,lumcd497,lumcnd3197,lumcd3197
      integer nnodiffiles,nmcfiles,mcinput,nmcfiles97,nnodiffiles97

      common /weighting/ lumcnd1, lumcd1, lumcnd2, lumcd2, lumcnd3,
     &     lumcd3, lumcnd4, lumcd4, lumda, fracdif, q2mc1,q2mc31,
     &     q2mc2, q2mc3, q2mc4, mcwtq2, f2weight, nnodiffiles,
     &     nmcfiles, mcinput,wt_vtx,lumda97,lumcd197,lumcnd197,
     &     lumcd297,lumcnd297,lumcd397,lumcnd397,lumcd497,
     &     lumcnd497,nmcfiles97,nnodiffiles97,fracdif97,lumcnd3197,
     &     lumcd3197,bpcor

c ------------------------------------------------
c Bins
c ------------------------------------------------
      real ymin, ymax, q2min, q2max, xmean, q2mean, logxmean,
     &bindat, binbgd, binmc, ysigz,mcraw,bakraw,datraw

      common /binstuff/ ymin, ymax, q2min, q2max,mcraw,bakraw,datraw,
     &xmean,q2mean,logxmean,bindat,binbgd,binmc,ysigz

c ------------------------------------------------
c Others
c ------------------------------------------------
      real acceptance, f2isr, fallm, f2mrsa, errors, en44m, q2tru, xtru,
     &deltaisr,logyelcorr,datanorm,backnorm,baknorm,zempzhad2,ysig_x_z
    
      common /others/ acceptance, f2isr, fallm, errors,
     &     f2mrsa, en44m, q2tru, xtru, deltaisr, logyelcorr, datanorm, 
     &     backnorm, baknorm,zempzhad2,ysig_x_z
                
c -------------------------------------------------
c Steering cards
c -------------------------------------------------

      real enelow,elglow,empztothigh,empztotlow,
     &empzcalcard,zvertex,bgdnormal,lumires,lumiscale,
     &hadscale,elenscale,ybinhigh,
     &ybinlow,srtdpressc,diffractiv

	  integer boxcutcard1, tag44mveto,offlineuse,bkspuse

      common /steeringcards/ enelow,elglow,empztothigh,empztotlow,
     &empzcalcard,zvertex,bgdnormal,lumires,lumiscale,
     &hadscale,elenscale,offlineuse,bkspuse,ybinhigh,
     &ybinlow,srtdpressc,diffractiv,boxcutcard1,tag44mveto         





