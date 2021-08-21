function k = findNumPrim(atomicOrbitals)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

k=0;
for i=1:length(atomicOrbitals)
    for j=1:length(atomicOrbitals{i}.primitives)
         k=k+1;
    end
end
end

