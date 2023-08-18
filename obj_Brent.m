function [f] = obj_Brent(sigmaN,F0,K,tau,C_real)
% Optimization objective function for Bachelier volatility under Brent algorithm
d = (F0 - K)./(sigmaN.*sqrt(tau));
C_pre = (F0 - K).*normcdf(d) + sigmaN.*sqrt(tau).*normpdf(d);
f = C_pre - C_real;

end

