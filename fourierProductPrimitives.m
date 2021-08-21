function value = fourierProductPrimitives(primitive1,primitive2,Q)
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

I = zeros(size(Q,1),3);
for i=1:3
    for lambda=0:a1(i)
        for mu=0:a2(i)
            if mod(lambda+mu,2)==0
            I(:,i) = I(:,i) + nchoosek(a1(i),lambda)*nchoosek(a2(i),mu)*...
                ((alpha1*(A1(i)-A2(i))+1i*Q(:,i)/2)./(alpha1+alpha2)).^(a2(i)-mu).*...
                ((alpha2*(A2(i)-A1(i))+1i*Q(:,i)/2)./(alpha1+alpha2)).^(a1(i)-lambda)*...
                sqrt(pi/(alpha1+alpha2))*doubleFactorial(lambda+mu-1)/(2*(alpha1+alpha2))^((lambda+mu)/2);
            end
        end
    end
    I(:,i) = I(:,i).*exp(...
    1/(alpha1+alpha2)*...
    (...
    - alpha1*alpha2*(A1(i)-A2(i))^2 ...
    + 1i*(alpha1*A1(i)+alpha2*A2(i))*Q(:,i)...
    - Q(:,i).^2/4 ...
    )...
    );
end
value = N1*N2*I(:,1).*I(:,2).*I(:,3);
end
