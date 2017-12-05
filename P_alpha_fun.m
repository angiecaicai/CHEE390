function [ P_alpha ] = P_alpha_fun(z,K,alphaP,N)
%   z is a vector containing z values from all speicies
%   K is a vector containing K values from all sepeicies
%   a is alpha
%   n is number of species
P_alpha =0;
for i=1:N
    P_alpha = P_alpha + z(i)*(K(i)-1)/(1+alphaP*(K(i)-1));
end
end

