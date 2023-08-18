function [f] = Bachelier_call(F0,K,tau,sigmaN)
% the European call option prices for the Bachelier model
d = (F0 - K)./(sigmaN.*sqrt(tau));
f = (F0 - K).*normcdf(d) + sigmaN.*sqrt(tau).*normpdf(d);

end