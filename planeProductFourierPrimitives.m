function value = planeProductFourierPrimitives(primitive1,primitive2,e,q)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
A1 = primitive1.A;
A2 = primitive2.A;
R = A1-A2;

a1 = primitive1.a; L1=sum(a1);
a2 = primitive2.a; L2=sum(a2);
b = a1+a2;

alpha1 = primitive1.alpha;
alpha2 = primitive2.alpha;
beta = (alpha1+alpha2)/(4*alpha1*alpha2);

N1 = primitive1.N;
N2 = primitive2.N;

I = zeros(size(q,1),2);

if isequal(e,[1 0 0])
    for lambda=0:b(2)
        if mod(lambda,2)==0
        I(:,1) = I(:,1) + nchoosek(b(2),lambda)*(-1i/(2*beta)*R(2))^(b(2)-lambda)...
            *sqrt(pi/beta)*doubleFactorial(lambda-1)/(2*beta)^((lambda)/2);
        end
    end
    for lambda=0:b(3)
        if mod(lambda,2)==0
        I(:,2) = I(:,2) + nchoosek(b(3),lambda)*(-1i/(2*beta)*R(3))^(b(3)-lambda)...
            *sqrt(pi/beta)*doubleFactorial(lambda-1)/(2*beta)^((lambda)/2);
        end
    end
    value = prod(I,2);
    value = value.*q.^(b(1)).*exp(-beta*(q+1i/(2*beta)*R(1)).^2-1i*pi/2*(L1-L2)-norm(R)^2*alpha1*alpha2/(alpha1+alpha2));
    
    
elseif isequal(e,[0 1 0])
    for lambda=0:b(1)
        if mod(lambda,2)==0
        I(:,1) = I(:,1) + nchoosek(b(1),lambda)*(-1i/(2*beta)*R(1))^(b(1)-lambda)...
            *sqrt(pi/beta)*doubleFactorial(lambda-1)/(2*beta)^((lambda)/2);
        end
    end
    for lambda=0:b(3)
        if mod(lambda,2)==0
        I(:,2) = I(:,2) + nchoosek(b(3),lambda)*(-1i/(2*beta)*R(3))^(b(3)-lambda)...
            *sqrt(pi/beta)*doubleFactorial(lambda-1)/(2*beta)^((lambda)/2);
        end
    end
     value = prod(I,2);
     value = value.*q.^(b(2)).*exp(-beta*(q+1i/(2*beta)*R(2)).^2-1i*pi/2*(L1-L2)-norm(R)^2*alpha1*alpha2/(alpha1+alpha2));
     
     
elseif isequal(e,[0 0 1])
    for lambda=0:b(1)
        if mod(lambda,2)==0
        I(:,1) = I(:,1) + nchoosek(b(1),lambda)*(-1i/(2*beta)*R(1))^(b(1)-lambda)...
            *sqrt(pi/beta)*doubleFactorial(lambda-1)/(2*beta)^((lambda)/2);
        end
    end
    for lambda=0:b(2)
        if mod(lambda,2)==0
        I(:,2) = I(:,2) + nchoosek(b(2),lambda)*(-1i/(2*beta)*R(2))^(b(2)-lambda)...
            *sqrt(pi/beta)*doubleFactorial(lambda-1)/(2*beta)^((lambda)/2);
        end
    end
     value = prod(I,2);
     value = value.*q.^(b(3)).*exp(-beta*(q+1i/(2*beta)*R(3)).^2-1i*pi/2*(L1-L2)-norm(R)^2*alpha1*alpha2/(alpha1+alpha2));
end



value = value * N1*N2/(2^(3+L1+L2)*alpha1^(3/2+L1)*alpha2^(3/2+L2));
end

