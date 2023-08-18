function [f] = charac(tau,phi,theta,kappa,sigma,rho,V)
% Characteristic function of the logarithmic return of Heston model to maturity
d = sqrt((kappa - rho.*sigma.*1i.*phi).^2 + sigma.^2.*(1i.*phi + phi.^2));
g2 = (kappa - rho.*sigma.*1i.*phi - d)./(kappa - rho.*sigma.*1i.*phi+ d);
f = exp(kappa.*theta./(sigma.^2).*((kappa - rho.*sigma.*1i.*phi - d).*tau - ...
    2.*log((1 - g2.*exp(-d.*tau))./(1 - g2)))).*...
    exp(V./(sigma.^2).*(kappa - rho.*sigma.*1i.*phi - d).*...
    (1 - exp(-d.*tau))./(1 - g2.*exp(-d.*tau)));

end