function [f] = charac_T(tau,phi,theta,kappa,sigma,rho,V)
%Sensitivity of the characteristic function of the logarithmic return of Heston model to maturity
f = -((exp(((theta).*(kappa).*(tau).*((kappa)-1i.*(sigma).*(rho).*(phi)-sqrt((sigma).^2.*(phi).*...
    (1i+(phi))+((kappa)-1i.*(sigma).*(rho).*(phi)).^2))-(V.*(sigma).^2.*(phi).*(1i+(phi)))./...
    ((kappa)-1i.*(sigma).*(rho).*(phi)+sqrt((sigma).^2.*(phi).*(1i+(phi))+((kappa)-1i.*(sigma).*...
    (rho).*(phi)).^2).*coth(1./2.*(tau).*sqrt((sigma).^2.*(phi).*(1i+(phi))+((kappa)-1i.*(sigma).*...
    (rho).*(phi)).^2)))-2.*(theta).*(kappa).*log((((kappa)-1i.*(sigma).*(rho).*(phi)+sqrt((sigma).^2.*...
    (phi).*(1i+(phi))+((kappa)-1i.*(sigma).*(rho).*(phi)).^2)).*(1+(exp(-(tau).*sqrt((sigma).^2.*...
    (phi).*(1i+(phi))+((kappa)-1i.*(sigma).*(rho).*(phi)).^2)).*(-(kappa)+1i.*(sigma).*(rho).*...
    (phi)+sqrt((sigma).^2.*(phi).*(1i+(phi))+((kappa)-1i.*(sigma).*(rho).*(phi)).^2)))./((kappa)-1i.*...
    (sigma).*(rho).*(phi)+sqrt((sigma).^2.*(phi).*(1i+(phi))+((kappa)-1i.*(sigma).*(rho).*...
    (phi)).^2))))./(2.*sqrt((sigma).^2.*(phi).*(1i+(phi))+((kappa)-1i.*(sigma).*(rho).*...
    (phi)).^2))))./(sigma).^2).*(-(kappa)+1i.*(sigma).*(rho).*(phi)+sqrt((sigma).^2.*(phi).*...
    (1i+(phi))+((kappa)-1i.*(sigma).*(rho).*(phi)).^2)).*(-(theta).*(kappa).*(sigma).^2.*(phi).*...
    (1i+(phi))+exp(2.*(tau).*sqrt((kappa).^2-2.*1i.*(kappa).*(sigma).*(rho).*(phi)+(sigma).^2.*...
    (phi).*(1i+(phi)+1i.^2.*(rho).^2.*(phi)))).*(theta).*(kappa).*(2.*(kappa).^2-4.*1i.*(kappa).*...
    (sigma).*(rho).*(phi)+(sigma).^2.*(phi).*(1i+(phi)+2.*1i.^2.*(rho).^2.*(phi))+2.*(kappa).*...
    sqrt((kappa).^2-2.*1i.*(kappa).*(sigma).*(rho).*(phi)+(sigma).^2.*(phi).*(1i+(phi)+1i.^2.*...
    (rho).^2.*(phi)))-2.*1i.*(sigma).*(rho).*(phi).*sqrt((kappa).^2-2.*1i.*(kappa).*(sigma).*...
    (rho).*(phi)+(sigma).^2.*(phi).*(1i+(phi)+1i.^2.*(rho).^2.*(phi))))+2.*exp((tau).*...
    sqrt((kappa).^2-2.*1i.*(kappa).*(sigma).*(rho).*(phi)+(sigma).^2.*(phi).*(1i+(phi)+1i.^2.*...
    (rho).^2.*(phi)))).*((kappa)-1i.*(sigma).*(rho).*(phi)+sqrt((kappa).^2-2.*1i.*(kappa).*...
    (sigma).*(rho).*(phi)+(sigma).^2.*(phi).*(1i+(phi)+1i.^2.*(rho).^2.*(phi)))).*((theta).*...
    (kappa).*(-(kappa)+1i.*(sigma).*(rho).*(phi))+V.*((kappa).^2-2.*1i.*(kappa).*(sigma).*(rho).*...
    (phi)+(sigma).^2.*(phi).*(1i+(phi)+1i.^2.*(rho).^2.*(phi))))))./((sigma).^2.*(-(kappa)+1i.*...
    (sigma).*(rho).*(phi)+exp((tau).*sqrt((sigma).^2.*(phi).*(1i+(phi))+((kappa)-1i.*(sigma).*...
    (rho).*(phi)).^2)).*((kappa)-1i.*(sigma).*(rho).*(phi))+sqrt((sigma).^2.*(phi).*(1i+(phi))+...
    ((kappa)-1i.*(sigma).*(rho).*(phi)).^2)+exp((tau).*sqrt((sigma).^2.*(phi).*(1i+(phi))+((kappa)-...
    1i.*(sigma).*(rho).*(phi)).^2)).*sqrt((sigma).^2.*(phi).*(1i+(phi))+((kappa)-1i.*(sigma).*(rho).*...
    (phi)).^2)).^2));

end













