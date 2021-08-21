function [SF,CP] = addNoise(SF,CP,percent)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% percent : percentage/100

if nargin < 3
  disp('Noise amplitude set to 3% by default.');
  percent = 0.03;
end

ind = logical([0;abs(SF.value(2:end,1))>=10^-(10)]);
%ind = logical([0;abs(SF.value(2:end,1))>=10^(-20)]);
sigmaSF = zeros(size(SF.value(:,1)));
sigmaSF(ind) = abs(SF.value(ind,1))*percent;
sigmaSF(~ind) = 1;

alpha = 1./(percent^2*CP.value(1,:));
sigmaCP = sqrt(CP.value)./sqrt(alpha);

SF.value(ind,1) = SF.value(ind,1) + sigmaSF(ind).*randn(size(SF.value(ind,1)));

CP.value = CP.value + sigmaCP.*randn(size(CP.value));

SF.sigma = sigmaSF;
CP.sigma = sigmaCP;


end

