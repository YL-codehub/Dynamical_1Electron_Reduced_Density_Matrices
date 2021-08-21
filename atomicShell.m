classdef atomicShell
% object including 'atomicOrbital' (see atomicOrbital class) objects and
% passage matrices (computed here) as attributes.
%
% YL.

    properties
        atomicOrbitals
        CartesianToSpherical
    end
    
    methods
         function thisAtomicShell = atomicShell(origin,atom,type,contractionCoefs,powerExponents)
            
             if type == 's' % already in cartesian coordinates
                 angularMomentum = [0 0 0];
                 thisAtomicShell.atomicOrbitals{1} = atomicOrbital(origin,atom,angularMomentum,contractionCoefs,powerExponents);
                 thisAtomicShell.CartesianToSpherical = 1;
             elseif type =='sp' % already in cartesian coordinates
                 angularMomentum = [0 0 0];
                 thisAtomicShell.atomicOrbitals{1} = atomicOrbital(origin,atom,angularMomentum,contractionCoefs(1,:),powerExponents);
                 angularMomentum = [1 0 0];
                 thisAtomicShell.atomicOrbitals{2} = atomicOrbital(origin,atom,angularMomentum,contractionCoefs(2,:),powerExponents);
                 angularMomentum = [0 1 0];
                 thisAtomicShell.atomicOrbitals{3} = atomicOrbital(origin,atom,angularMomentum,contractionCoefs(2,:),powerExponents);
                 angularMomentum = [0 0 1];
                 thisAtomicShell.atomicOrbitals{4} = atomicOrbital(origin,atom,angularMomentum,contractionCoefs(2,:),powerExponents);
                 thisAtomicShell.CartesianToSpherical = eye(4);
             elseif type =='p' % already in cartesian coordinates
                 angularMomentum = [1 0 0];
                 thisAtomicShell.atomicOrbitals{1} = atomicOrbital(origin,atom,angularMomentum,contractionCoefs,powerExponents);
                 angularMomentum = [0 1 0];
                 thisAtomicShell.atomicOrbitals{2} = atomicOrbital(origin,atom,angularMomentum,contractionCoefs,powerExponents);
                 angularMomentum = [0 0 1];
                 thisAtomicShell.atomicOrbitals{3} = atomicOrbital(origin,atom,angularMomentum,contractionCoefs,powerExponents);
                 thisAtomicShell.CartesianToSpherical = eye(3);
             elseif type =='d'
                 angularMomentum = [2 0 0];
                 thisAtomicShell.atomicOrbitals{1} = atomicOrbital(origin,atom,angularMomentum,contractionCoefs,powerExponents);
                 angularMomentum = [1 1 0];
                 thisAtomicShell.atomicOrbitals{2} = atomicOrbital(origin,atom,angularMomentum,contractionCoefs,powerExponents);
                 angularMomentum = [1 0 1];
                 thisAtomicShell.atomicOrbitals{3} = atomicOrbital(origin,atom,angularMomentum,contractionCoefs,powerExponents);
                 angularMomentum = [0 2 0];
                 thisAtomicShell.atomicOrbitals{4} = atomicOrbital(origin,atom,angularMomentum,contractionCoefs,powerExponents);
                 angularMomentum = [0 1 1];
                 thisAtomicShell.atomicOrbitals{5} = atomicOrbital(origin,atom,angularMomentum,contractionCoefs,powerExponents);
                 angularMomentum = [0 0 2];
                 thisAtomicShell.atomicOrbitals{6} = atomicOrbital(origin,atom,angularMomentum,contractionCoefs,powerExponents);
                 thisAtomicShell.CartesianToSpherical = [...
                     -1/2           0   0   -1/2        0   1;
                      0             0   1    0          0   0;
                      0             0   0    0          1   0;
                      sqrt(3)/2     0   0   -sqrt(3)/2  0   0;
                      0             1   0    0          0   0;...
                      ];
             else
                 warning('not yet implemented for f')
             end
       
        end
    end
end

