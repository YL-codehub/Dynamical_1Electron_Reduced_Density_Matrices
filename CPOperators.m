function CPOp = CPOperators(atomicOrbitals,C,e,q,R)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for i=1:size(e,1)
    e(i,:) = (R.'*(e(i,:).')).'; % .' = transpose
end

CPOp = zeros(length(atomicOrbitals),length(atomicOrbitals),size(e,1),length(q));

for i = 1:length(atomicOrbitals)
    for j = i:length(atomicOrbitals)
        for k=1:size(e,1)
            CPOp(i,j,k,:) = fourierAutoCorrelationAtomicOrbitals(atomicOrbitals{i},atomicOrbitals{j},e(k,:),q);
            CPOp(j,i,k,:) = conj(CPOp(i,j,k,:));
        end
    end
end

for i=1:size(e,1)
    for j=1:length(q)
    CPOp(:,:,i,j)=C.'*CPOp(:,:,i,j)*C;
    end
end


