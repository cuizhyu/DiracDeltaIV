function [sigmaN] = Jackel_Bachelier_IV(F0,K,tau,C_ture)
% the rational approximation formula for volatility under the Bachelier model proposed by Jackel(2017)
eta_C = -0.001882039271;
eta_0 = -abs(C_ture - max(F0 - K,0))./abs(K - F0);

g = 1./(eta_0 - 1/2);
epsilon_1 = (0.032114372355 - g.^2.*(0.016969777977 - g.^2.*(...
    2.6207332461e-3 - 9.6066952861e-5.*g.^2)))./(1 - ...
    g.^2.*(0.6635646938 - g.^2.*(0.14528712196 - 0.010472855461.*g.^2)));
x_1_a = g.*(1./sqrt(2.*pi) + epsilon_1.*g.^2);

h = sqrt(-log(-eta_0));
x_1_b = (9.4883409779 - h.*(9.6320903635 - h.*(0.58556997323 + 2.1464093351.*h)))./...
    (1 - h.*(0.65174820867 + h.*(1.5120247828 + 6.6437847132e-5.*h)));

x_1 = (eta_0 < eta_C).*x_1_a + (1 - (eta_0 < eta_C)).*x_1_b;

q = (normcdf(x_1) + normpdf(x_1)./x_1 - eta_0)./normpdf(x_1);
x = x_1 + (3.*q.*x_1.^2.*(2 - q.*x_1.*(2 + x_1.^2)))./...
    (6 + q.*x_1.*(-12 + x_1.*(6.*q + x_1.*(-6 + q.*x_1.*(3 + x_1.^2)))));

sigmaN = abs(F0 - K)./abs(x.*sqrt(tau));


end