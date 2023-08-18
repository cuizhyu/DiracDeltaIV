function [f] = European_call(S0,K,r,tau,sigma_all,q)
% the European call option prices for the BS model
f_d1 = (log(S0./K) + (r-q+sigma_all.^2./2).*tau)./(sigma_all.*sqrt(tau));
f_d2 = f_d1 - sigma_all.*sqrt(tau);
f = S0.*exp(-q.*tau).*normcdf(f_d1,0,1) - K.*exp(-r.*tau).*normcdf(f_d2,0,1);

end