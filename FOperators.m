function SFop = FOperators(atomicOrbitals,C,Q,R,T)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Qr = zeros(size(Q));

for i=1:size(Q,1)
    Qr(i,:)=(R.'*(Q(i,:).')).';
end

SFop = zeros(length(atomicOrbitals),length(atomicOrbitals),size(Q,1));

for i = 1:length(atomicOrbitals)
    for j = i:length(atomicOrbitals)
        SFop(i,j,:) = fourierProductAtomicOrbitals(atomicOrbitals{i},atomicOrbitals{j},Qr);
        SFop(j,i,:) = conj(SFop(i,j,:));
    end
end

for i=1:size(Q,1)
    SFop(:,:,i)=C.'*SFop(:,:,i)*C;
    exp(1i*dot(T,Q(i,:)));
    Q(i,:);
    i;
    SFop(:,:,i)=SFop(:,:,i)*exp(1i*dot(T,Q(i,:)));
end

