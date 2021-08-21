classdef primitive
% object representing a gaussian primitive function, made of the following
% attributes :
% - .A : coordinates of the gaussian function
% - .a : exponents of (X-A_X), X = x,y or z
% - .alpha : exponent of the gaussian function.
% - .N : normalization coefficient 
%
% YL.

    
    properties
        A
        a
        alpha
    end
    
    properties (Dependent)
        N
    end
    
    methods
        
        function thisPrimitive = primitive(origin,angularMomentum,powerExponent)
            
        if nargin == 3
            thisPrimitive.A = origin;
            thisPrimitive.a = angularMomentum;
            thisPrimitive.alpha = powerExponent;
        end
        
        end
        
        function N = get.N( thePrimitive )
             a1 = thePrimitive.a(1);
             a2 = thePrimitive.a(2);
             a3 = thePrimitive.a(3);
             d = a1 + a2 + a3;
             factor = sqrt(((2^d)*factorial(a1)*factorial(a2)*factorial(a3))/((factorial(2*a1)*factorial(2*a2)*factorial(2*a3))));
             N = factor*((4*thePrimitive.alpha)^(d/2))*((2*thePrimitive.alpha/pi)^(3/4));
        end
        
        function value = eval( thePrimitive, r )
            A1 = thePrimitive.A(1);
            A2 = thePrimitive.A(2);
            A3 = thePrimitive.A(3);
            a1 = thePrimitive.a(1);
            a2 = thePrimitive.a(2);
            a3 = thePrimitive.a(3);
            value = thePrimitive.N ...
                .*(r(:,1)-A1).^a1 ...
                .*(r(:,2)-A2).^a2 ...
                .*(r(:,3)-A3).^a3 ...
                .*exp(-thePrimitive.alpha*sum((r-thePrimitive.A).^2,2));
        end
        
        function value = evalFourierTransform( thePrimitive, p )
            
            I=zeros(size(p,1),3);
            for i=1:3
                for lambda=0:thePrimitive.a(i)
                    if mod(lambda,2)==0
                    I(:,i) = I(:,i) ...
                        + nchoosek(thePrimitive.a(i),lambda)*(-1i*p(:,i)/(2*thePrimitive.alpha)).^(thePrimitive.a(i)-lambda)...
                        .*sqrt(pi/thePrimitive.alpha)*doubleFactorial(lambda-1)/(2*thePrimitive.alpha)^((lambda)/2);
                    end
                end
                I(:,i) = exp(-p(:,i).^2/(4*thePrimitive.alpha)-1i*p(:,i)*thePrimitive.A(i)).*I(:,i);
            end
             
            value = 1/sqrt(2*pi)^3*thePrimitive.N*prod(I,2);
        end
    end
end

