c -----------------------------------------------
      subroutine isr_histos97(histofill)
c -----------------------------------------------

      implicit none

#include "common97.inc"
#include "local.inc"
#include "constant.inc"
                      
      integer histofill
      real datawt, bakwt, mcwt

      real energy
      real zisr,zisrtrue,rfl,epsilon,flcor,flwt
      real temp

      real flcor2,flcor3,flwt2,flwt3

      if (acceptance.ne.0) then
         datawt = 1. / acceptance
         bakwt = baknorm / acceptance
      else 
         datawt = 0.
         bakwt = 0.
      end if

      zisr = (27.52-(elumig))/27.52
      zisrtrue=(27.52-hcl_eneg)/27.52
      rfl = 1.4
      epsilon = 2*(1-((y_tru)/zisrtrue))/(1+((1-(y_tru)/zisrtrue))**2)
      flcor = (1+epsilon*rfl)/(1+rfl)
      rfl = 0.2
      flcor2 = (1+epsilon*rfl)/(1+rfl)
      rfl = 10
      flcor3 = (1+epsilon*rfl)/(1+rfl)

      mcwt = mcwtq2 * f2weight * wt_vtx * bpcor

      flwt= mcwt*flcor
      flwt2 = mcwt*flcor2
      flwt3 = mcwt*flcor3

      energy = 0
      temp = elumig - ene_lg

c ---------------------------------        
c  Data
c ---------------------------------
      if (histofill.eq.1) then
         call hf1(100,empzcal,datawt)
         call hf1(101,empztot,datawt)	  
         call hf1(108,CAL_EE,datawt)  
         call hf1(109,logq2el,datawt)
         call hf1(111,logyel,datawt)
         call hf1(112,logyelcorr,datawt)
         call hf1(113,logysig,datawt)
         call hf1(114,logyjbz,datawt)
         call hf1(116,logxel,datawt)
         call hf2(117,logxel,logq2el,1.)
         call hf1(118,VCT_ZVC,datawt)
         call hf1(120,electron_en,datawt)
         call hf1(121,logq2elc,datawt)
         call hf1(122,elumig,datawt)
         call hf1(129,gamma,datawt)
         call hf1(130,gamma2,datawt)
         call hf1(131,ETAMAX,datawt)
         call hf1(132,deltaisr,datawt)
         call hf2(150,bestx,besty,1.)
         call hf1(151,bestx,datawt)
         call hf1(152,besty,datawt)
         call hf1(153,best_th,datawt)
         call hf1(160,ysig_x_z,datawt)
         call hf1(161,ysigz,datawt)
         call hf1(162,y_sig,datawt)
         call hf1(163,y_sig,1.)
         call hf1(171,logyel,1.)
         call hf1(172,logyelcorr,1.)
         call hf1(173,logysig,1.)
         call hf1(174,logyjbz,1.)
         call hf1(176,zempzhad,datawt)

c ------------------------------------
c Background
c ------------------------------------
      else if (histofill.eq.2) then
         call hf1(200,empzcal,bakwt)
         call hf1(201,empztot,bakwt)   
         call hf1(208,CAL_EE,bakwt)
         call hf1(209,logq2el,bakwt)
         call hf1(212,logyel,bakwt)
         call hf1(213,logysig,bakwt)
         call hf1(214,logyjbz,bakwt)
         call hf1(216,logxel,bakwt)
         call hf2(217,logxel,logq2el,1.)
         call hf1(218,VCT_ZVC,bakwt)
         call hf1(220,electron_en,bakwt)
         call hf1(221,logq2elc,bakwt)
         call hf1(222,elumig,bakwt)
         call hf1(223,ENE_LG,bakwt)
         call hf1(224,TEMPLUMG,bakwt)
         call hf1(229,gamma2,bakwt)
         call hf1(230,gamma2,bakwt)
         call hf1(231,ETAMAX,bakwt)
         call hf1(232,deltaisr,bakwt)
         call hf1(153,best_th,bakwt)
         call hf1(260,ysig_x_z,bakwt)
         call hf1(261,ysigz,bakwt)
         call hf1(262,y_sig,bakwt)
         call hf1(263,y_sig,1.)
         call hf1(272,logyel,1.)
         call hf1(273,logysig,1.)
         call hf1(274,logyjbz,1.)
         call hf1(276,zempzhad,bakwt)

c ------------------------------------   
c MonteCarlo
c ------------------------------------        
      else if (histofill.eq.3) then
         call hf1(300,empzcal,mcwt)
         call hf1(301,empztot,mcwt)	  
         call hf1(308,CAL_EE,mcwt)
         call hf1(309,logq2el,mcwt)
         call hf1(310,logysig,flwt2)   
         call hf1(311,logysig,flwt3)   
         call hf1(312,logyel,mcwt)
         call hf1(313,logysig,mcwt)
         call hf1(314,logyjbz,mcwt)
         call hf1(315,logysig,flwt)   
         call hf1(316,logxel,mcwt)
         call hf2(317,logxel,logq2el,1.)
         call hf1(318,VCT_ZVC,mcwt)
         call hf1(320,electron_en,mcwt)
         call hf1(321,logq2elc,mcwt)
         call hf1(322,elumig,mcwt) 
         call hf1(323,logysig,mcwt)
         call hf1(329,gamma2,mcwt)
         call hf1(330,gamma2,mcwt)
         call hf1(331,ETAMAX,mcwt)
         call hf1(332,deltaisr,mcwt)
         call hf1(353,best_th,mcwt)
         call hf1(360,ysig_x_z,mcwt)
         call hf1(361,ysigz,mcwt)
         call hf1(362,y_sig,mcwt)
         call hf1(363,y_sig,1.)
         call hf1(364,y_sig,flwt)
         call hf1(372,logyel,1.)
         call hf1(373,logysig,1.)
         call hf1(374,logyjbz,1.)
         call hf1(376,zempzhad,mcwt)
         call hf1(377,electron_en,1.)	
         call hf1(378,elumig,1.)	    

c -------------------------------------
c MC True
c -------------------------------------
         call hf1(401,MCELE,mcwt)
         call hf1(402,HCL_ENEG,mcwt)
         call hf1(403,HCL_ENEE,mcwt)
         call hf1(404,logmcq2,mcwt)
         call hf1(405,logmcx,mcwt)
         call hf1(406,logmcy,mcwt)
         call hf1(408,logmcx,1.)
         call hf1(409,logmcy,1.)
         call hf1(410,MC_ZV,1.)
         call hf1(462,temp,1.)
      end if

      return
      end
