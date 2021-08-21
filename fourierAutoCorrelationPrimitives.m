



function value = fourierAutoCorrelationPrimitives(primitive1,primitive2,e,q)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
A1 = primitive1.A;
A2 = primitive2.A;
R = A1-A2;

a1 = primitive1.a; 
a2 = primitive2.a; 

alpha1 = primitive1.alpha;
alpha2 = primitive2.alpha;
beta = alpha1*alpha2/(alpha1+alpha2);
gammax = @(lambda,mu) nchoosek(a1(1),lambda)*nchoosek(a2(1),mu)*(-alpha1/(alpha2+alpha1))^(a2(1)-mu)*(alpha2/(alpha2+alpha1))^(a1(1)-lambda);
gammay = @(lambda,mu) nchoosek(a1(2),lambda)*nchoosek(a2(2),mu)*(-alpha1/(alpha2+alpha1))^(a2(2)-mu)*(alpha2/(alpha2+alpha1))^(a1(2)-lambda);
gammaz = @(lambda,mu) nchoosek(a1(3),lambda)*nchoosek(a2(3),mu)*(-alpha1/(alpha2+alpha1))^(a2(3)-mu)*(alpha2/(alpha2+alpha1))^(a1(3)-lambda);

deltax = @(lambda,mu) a2(1) - lambda + a1(1) - mu;
deltay = @(lambda,mu) a2(2) - lambda + a1(2) - mu;
deltaz = @(lambda,mu) a2(3) - lambda + a1(3) - mu;

LAMBDA = @(lambda,mu) sqrt(pi/(alpha1+alpha2))*doubleFactorial(lambda+mu-1)/(2*(alpha1+alpha2))^((lambda+mu)/2);
N1 = primitive1.N;
N2 = primitive2.N;

value=0;
for lambdax=0:a1(1)
    for mux=0:a2(1)
        if mod(lambdax+mux,2)==0
            for lambday=0:a1(2)
                 for muy=0:a2(2)
                     if mod(lambday+muy,2)==0
                        for lambdaz=0:a1(3)
                            for muz=0:a2(3)
                                if mod(lambdaz+muz,2)==0
                                    I=0;
                                    for kappax=0:deltax(lambdax,mux)
                                        for kappay=0:deltay(lambday,muy)
                                            for kappaz=0:deltaz(lambdaz,muz)
                                                if mod(kappax+kappay+kappaz,2)==0
                                                    I = I +  nchoosek(deltax(lambdax,mux),kappax)*nchoosek(deltay(lambday,muy),kappay)*nchoosek(deltaz(lambdaz,muz),kappaz)...
                                                        *(e(1)*(dot(e,R)-1i*q/(2*beta))-R(1)).^(deltax(lambdax,mux)-kappax)...
                                                        .*(e(2)*(dot(e,R)-1i*q/(2*beta))-R(2)).^(deltay(lambday,muy)-kappay)...
                                                        .*(e(3)*(dot(e,R)-1i*q/(2*beta))-R(3)).^(deltaz(lambdaz,muz)-kappaz)...
                                                        *e(1)^kappax*e(2)^kappay*e(3)^kappaz...
                                                        *sqrt(pi/beta)*doubleFactorial(kappax+kappay+kappaz-1)/(2*beta)^((kappax+kappay+kappaz)/2);
                                                end
                                            end
                                        end
                                    end        
                                    H = 1/(2*pi)*exp(-beta*(norm(R)^2-(dot(R,e)-1i*q/(2*beta)).^2)).*I;
                                    value = value + gammax(lambdax,mux)*gammay(lambday,muy)*gammaz(lambdaz,muz) ...
                                        *LAMBDA(lambdax,mux)*LAMBDA(lambday,muy)*LAMBDA(lambdaz,muz)*H;    
                                end
                            end
                        end
                     end
                 end
            end
        end
    end
end

value = N1*N2*value;

end

