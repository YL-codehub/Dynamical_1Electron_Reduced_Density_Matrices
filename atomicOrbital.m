classdef atomicOrbital
% Object representing an atomic orbital through its
% characteristics : contraction coefficients in line .d, associated
% primitives list .primitives (usually gaussian functions, see primitive class file),
% and atomic number .atomType.
%
% YL.
    
    properties
        d
        primitives
        atomType
    end
    
    methods
         function thisAtomicOrbital = atomicOrbital(origin,atom,angularMomentum,contractionCoefs,powerExponents)
            thisAtomicOrbital.atomType = atom;
            if nargin == 5
                thisAtomicOrbital.d = contractionCoefs;
                for i=1:length(contractionCoefs)
                    thisAtomicOrbital.primitives{i} = primitive(origin,angularMomentum,powerExponents(i));
                end
            end
        
         end
        
         function value = eval( theAtomicOrbital, r )
             
            value = 0;
             for i=1:length(theAtomicOrbital.d)
                  value = value + theAtomicOrbital.d(i)*theAtomicOrbital.primitives{i}.eval(r);
             end
               
         end
         
         function value = evalFourierTransform( theAtomicOrbital, p )
             
            value = 0;
             for i=1:length(theAtomicOrbital.d)
                  value = value + theAtomicOrbital.d(i)*theAtomicOrbital.primitives{i}.evalFourierTransform(p);
             end
               
         end
    end
end

