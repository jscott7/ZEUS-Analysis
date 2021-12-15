      real function fitfllin(x)
      common/pawpar/par(2)

      fitfllin=par(2)*((1+((2*(1-1.18*x))/
     &     (1+(1-1.18*x)**2))*par(1))/(1+par(1)))

      end
