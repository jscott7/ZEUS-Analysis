C
C     =========================================================== 
C

      SUBROUTINE f2allm(x,q2,f2)             


      REAL M02,M32,LAM2,M22
      DIMENSION PAR(35)
      DATA ALFA /112.2 /
      DATA XMP2 /0.8802/
      REAL X,Q2,F2

C      x=10.0**x
      W2=q2*(1./x -1.)+xmp2
      W=sqrt(w2)

C  POMERON
      S11   =   0.28067
      S12   =   0.22291
      S13   =   2.1979
      A11   =  -0.0808 
      A12   =  -0.44812 
      A13   =   1.1709
      B11   =   0.60243
      B12   =   1.3754
      B13   =   1.8439
      M12   =  49.457
 
C  REGGEON
      S21   =   0.80107
      S22   =   0.97307
      S23   =   3.4942  
      A21   =   0.58400
      A22   =   0.37888
      A23   =   2.6063
      B21   =   0.10711
      B22   =   1.9386
      B23   =   0.49338
      M22   =   0.15052
C
      M02   =   0.31985
      LAM2  =   0.065270
      Q02   =   0.46017 +LAM2
C     
      S=0.     
      IF(Q2.GT.0.)S=LOG(LOG((Q2+Q02)/LAM2)/LOG(Q02/LAM2))
                      
      Z=1.  
                 
                                                   
      IF(Q2.GT.0.)Z=(W2-XMP2)/(Q2+W2-XMP2)      
      IF(S.EQ.0.) THEN                        
                                  
C                                       
C   POMERON                               
C                                 
      XP=1./(1.+(W2-XMP2)/(Q2+M12))     
      AP=A11                           
      BP=B11**2                        
      SP=S11                              
      F2P=SP*XP**AP*Z**BP                  
C                                           
C   REGGEON                                  
C                                          
      XR=1./(1.+(W2-XMP2)/(Q2+M22))      
      AR=A21                             
      BR=B21**2                    
      SR=S21                       
      F2R=SR*XR**AR*Z**BR            
C                                    
      ELSE                           
C                                    
C   POMERON                             
C                               
      XP=1./(1.+(W2-XMP2)/(Q2+M12))   
      AP=A11+(A11-A12)*(1./(1.+S**A13)-1.)   
      BP=B11**2+B12**2*S**B13              
      SP=S11+(S11-S12)*(1./(1.+S**S13)-1.)   
      F2P=SP*XP**AP*Z**BP                  
C                                          
C   REGGEON                          
C                                       
      XR=1./(1.+(W2-XMP2)/(Q2+M22))    
      AR=A21+A22*S**A23                
      BR=B21**2+B22**2*S**B23         
      SR=S21+S22*S**S23               
      F2R=SR*XR**AR*Z**BR             
C                                     
      ENDIF
      CIN=ALFA/(Q2+M02)*(1.+4.*XMP2*Q2/(Q2+W2-XMP2)**2)/Z  
      SIGal=CIN*(F2P+F2R)                                              
      f2=sigal/alfa*(q2**2*(1.-x))/(q2+4.*xmp2*x**2)

      RETURN                                    
      END                               
C


