function [C] = Call_SINC(S0,K,r,q,tau,theta,kappa,sigma,rho,V,N,Xh0)
% utilize the SINC algorithm to calculate European call options under the Heston model
X_min = Xh0;
X_max = -X_min;
n = 1:N;
k = log(K./S0)-(r-q).*tau;
kn = (2.*n-1)./(X_max - X_min);
f1 = charac(tau,2.*pi.*kn,theta,kappa,sigma,rho,V);
f2 = charac(tau,2.*pi.*kn - 1i,theta,kappa,sigma,rho,V);
P = 0.5.*(K.*exp(-r.*tau) - S0.*exp(-q.*tau)) + ...
    2./pi.*sum(1./(2.*n - 1).*(sin(2.*pi.*k.*kn).*(K.*exp(-r.*tau).*real(f1) - ...
    S0.*exp(-q.*tau).*real(f2)) - ...
    cos(2.*pi.*k.*kn).*(K.*exp(-r.*tau).*imag(f1) - S0.*exp(-q.*tau).*imag(f2))));
C = P + S0.*exp(-q.*tau) - K.*exp(-r.*tau);

end






