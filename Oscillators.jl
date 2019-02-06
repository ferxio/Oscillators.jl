#==== OSCILADORES =====#
   
module Oscillators

   export Sprott1
   export Sprott2
   export Sprott3
   export Sprott4
   export Sprott5
   export Sprott6
   export Sprott7
   export Sprott8
   export Sprott9
   export Sprott10
   export Sprott11
   export Sprott12
   export Sprott13
   export Sprott14
   export Sprott15
   export Sprott16
   export Sprott17
   export Sprott18
   export Sprott19
      
   export Rossler
   export Lorenz
   export Logistic
   
   export AllSprott
   export AllOscillators

   
   ### Osciladores de Sprott
   function Sprott1(U)
      x, y, z = U[1], U[2], U[3]
      return [y,  -x.+y.*z,  1.0.-y.^2]
   end
   
   function Sprott2(U)
      x, y, z = U[1], U[2], U[3]
      return [y.*z,  x.-y,  1.0.-x.*y]
   end
   
   function Sprott3(U)
      x, y, z = U[1], U[2], U[3]
      return [y.*z,  x.-y,  1.0.-x.^2]
   end
   
   function Sprott4(U)
      x, y, z = U[1], U[2], U[3]
      return [-y,  x.+z,  x.*z.+3.0.*y.^2]
   end
   
   function Sprott5(U)
      x, y, z = U[1], U[2], U[3]
      return ([y.*z,  x.^2 .- y,  1.0.-4.0.*x])
   end

   function Sprott6(U)
      x, y, z = U[1], U[2], U[3]
      return ([y.+z,  -x.+0.5.*y,  x.^2.0 .- z])
   end
   
   function Sprott7(U)
      x, y, z = U[1], U[2], U[3]
      return ([0.4.*x.+z,  x.*z.-y,  -x.+y])
   end
   
   function Sprott8(U)
      x, y, z = U[1], U[2], U[3]
      return ([-y.+z.^2,  x.+0.5.*y,  x.-z])
   end
   
   function Sprott9(U)
      x, y, z = U[1], U[2], U[3]
      return ([0.2.*y,  x.+z,  x.+y.^2.0.-z])
   end
   
   function Sprott10(U)
      x, y, z = U[1], U[2], U[3]
      return ([2.0.*z,  -2.0.*y.+z,  -x.+y.+y.^2])
   end
   
   function Sprott11(U)
      x, y, z = U[1], U[2], U[3]
      return ([x.*y.-z,  x.-y,  x.+0.3.*z])
   end
   
   function Sprott12(U)
      x, y, z = U[1], U[2], U[3]
      return ([y.+3.9.*z,  0.9.*x.^2.0.-y,  1.0.-x])
   end
   
   function Sprott13(U)
      x, y, z = U[1], U[2], U[3]
      return ([-z,  -x.^2.0.-y,  1.7.+1.7.*x.+y])
   end
   
   function Sprott14(U)
      x, y, z = U[1], U[2], U[3]
      return ([-2.0.*y,  x.+z.^2,  1.0.+y.-2.0.*z])
   end
   
   function Sprott15(U)
      x, y, z = U[1], U[2], U[3]
      return ([y,  x.-z,  x.+z.*x.+2.7.*y])
   end
   
   function Sprott16(U)
      x, y, z = U[1], U[2], U[3]
      return ([2.7.*y.+x, -x.+y.^2,  x.+y])
   end
   
   function Sprott17(U)
      x, y, z = U[1], U[2], U[3]
      return ([-z,  -y,  3.1.*x.+y.^2+0.5.*z])
   end
   
   function Sprott18(U)
      x, y, z = U[1], U[2], U[3]
      return ([0.9.-y,  0.4.+z,  x.*y.-z])
   end
   
   function Sprott19(U)
      x, y, z = U[1], U[2], U[3]
      return ([-x.-4.0.*y,  x.+z.^2,  1.0.+z])
   end
   
   # Osciladores tipicos
   function Rossler(U,a=0.1,b=0.1,c=14.)
      return [-U[2].-U[3],  U[1].-a.*U[2],  b.+(U[1].-c).*U[3]]
   end

   function Lorenz(U,s=10.0,b=8/3,r=24.74)
      return [s.*(U[2].-U[1]),  -U[1].*U[3].+r.*U[1].-U[2],  U[1].*U[2].-b.*U[3]]
   end
   
   ###===== MAPEOS =====#
   function Logistic(x,r=3.99)
      #return Any[r.*x[1].*(1-x[1])]
      return -r.*x.*(x.-1.0)
   end

   AllSprott = [Sprott1 Sprott2 Sprott3 Sprott4 Sprott5 Sprott6 Sprott7 Sprott8 Sprott9 Sprott10 Sprott11 Sprott12 Sprott13 Sprott14 Sprott15 Sprott16 Sprott17 Sprott18 Sprott19]
   
   AllOscillators = [AllSprott Rossler Lorenz]
   
end # module
