function value = fourierAutoCorrelationAtomicOrbitals(AO1,AO2,e,q)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
d1 = AO1.d;
d2 = AO2.d;

primitives1 = AO1.primitives;
primitives2 = AO2.primitives;

value = 0;
    for i=1:length(d1)
        for j=1:length(d2)
            value = value + d1(i)*d2(j)*fourierAutoCorrelationPrimitives(primitives1{i},primitives2{j},e,q);
        end
    end
end

