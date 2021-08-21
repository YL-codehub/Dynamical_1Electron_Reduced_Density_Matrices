    function value = fourierProductAtomicOrbitals(AO1,AO2,Q)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
d1 = AO1.d;
d2 = AO2.d;

primitives1 = AO1.primitives;
primitives2 = AO2.primitives;

value = 0;
    for i=1:length(d1)
        for j=1:length(d2)
            value = value + d1(i)*d2(j)*fourierProductPrimitives(primitives1{i},primitives2{j},Q);
        end
    end
end
