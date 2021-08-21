function value = fourierConvolution(primitive1,primitive2,e,q)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
A1 = primitive1.A;
A2 = primitive2.A;

a1 = primitive1.a;
a2 = primitive2.a;

alpha1 = primitive1.alpha;
alpha2 = primitive2.alpha;

N1 = primitive1.N;
N2 = primitive2.N;

beta = alpha1*alpha2/(alpha1+alpha2);

end

