function C = orthogonalisation(atomicOrbitals)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

OverlapMatrix = zeros(length(atomicOrbitals));

for i = 1:length(atomicOrbitals)
    for j = i:length(atomicOrbitals)
        OverlapMatrix(i,j) = fourierProductAtomicOrbitals(atomicOrbitals{i},atomicOrbitals{j},[0 0 0]);
        OverlapMatrix(j,i) = OverlapMatrix(i,j);
    end
end
C = OverlapMatrix^(-1/2);
end

