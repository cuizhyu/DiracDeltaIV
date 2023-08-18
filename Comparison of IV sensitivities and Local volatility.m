%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Tighter bounds for implied volatility based on the Dirac delta family method
% Author: Zhenyu Cui,Yanchu Liu and Yuhang Yao
% Date: August 13, 2023
% Note: This document corresponds to secton 3.3-3.4 of the paper, with figure(11)-(13)
% MATLAB version: Matlab2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;
format long 
% initial parameters
S0 = 100;K = 100;r = 0.03;q = 0;
tau = 1/12;%1/12,,1/2

% heston model parameters
kappa = 1.5;theta = 0.2;sigma = 0.25;rho = -0.8;V = 0.2;

% SINC algorithm parameters
N0 = 2^9;% Number of items in the expansion
% The upper boundary of the truncated integra
% where it is the lower boundary of the logarithmic rate of return
% Tighter in the case of short maturity, wider in the case of long maturity
Xh0 = -10*(tau >= 0.5) + -2*(tau < 0.5);

% generate discrete points
moneyness_all = -0.2:0.4/20:0.2;
K = S0./exp(moneyness_all);

% Pre allocate memory to related variables
C_SINC = zeros(1,length(K));
IV_SINC = C_SINC;C_skew_finite_difference = C_SINC;IV_skew_finite_difference = C_SINC;
C_skew_SINC = C_SINC;C_convexity_SINC = C_SINC;C_convexity_finite_difference = C_SINC;
IV_convexity_finite_difference = C_SINC;C_time_SINC = C_SINC;
IV_chenxu = C_SINC;IV_skew_Chenxu = C_SINC;IV_convexity_Chenxu = C_SINC; 

for i = 1:length(K)
X = log(S0./K(i));
% utilize SINC algorithm to calculate option prices
C_SINC(i) = Call_SINC(S0,K(i),r,q,tau,theta,kappa,sigma,rho,V,N0,Xh0);
IV_SINC(i) = blsimpv(S0,K(i),r,tau,C_SINC(i),'Limit',1,'Yield',q,'Class', {'Call'},'method','jackel2016');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% untilize finite difference(central difference) to compute IV sensitivity
dX = 0.0001;
K1 = S0./exp(X + dX);
C_SINC1 = Call_SINC(S0,K1,r,q,tau,theta,kappa,sigma,rho,V,N0,Xh0);
IV_SINC1 = blsimpv(S0,K1,r,tau,C_SINC1,'Limit',1,'Yield',q,'Class', {'Call'},'method','jackel2016');
K2 = S0./exp(X - dX);
C_SINC2 = Call_SINC(S0,K2,r,q,tau,theta,kappa,sigma,rho,V,N0,Xh0);
IV_SINC2 = blsimpv(S0,K2,r,tau,C_SINC2,'Limit',1,'Yield',q,'Class', {'Call'},'method','jackel2016');

C_skew_finite_difference(i) = (C_SINC1 - C_SINC2)/(2.*dX);
IV_skew_finite_difference(i) = (IV_SINC1 - IV_SINC2)/(2.*dX);
C_convexity_finite_difference(i) = (C_SINC1 + C_SINC2 - 2.*C_SINC(i))/(dX.^2);
IV_convexity_finite_difference(i) = (IV_SINC1 + IV_SINC2 - 2.*IV_SINC(i))/(dX.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% utilize SINC algorithm to calculate option price sensitivities
X_min = Xh0;
X_max = -X_min;
n = 1:N0;
kn = (2.*n-1)./(X_max - X_min);
f1 = charac(tau,2.*pi.*kn,theta,kappa,sigma,rho,V);
f2 = charac(tau,2.*pi.*kn - 1i,theta,kappa,sigma,rho,V);
C_skew_SINC(i) = 0.5.*S0.*exp(-(X + r.*tau)) + 2./pi.*sum(1./(2.*n - 1).*...
    (-2.*pi.*kn.*cos(2.*pi.*(X + r.*tau).*kn).*(real(S0.*exp(-(X + r.*tau)).*f1 - S0.*f2)) + ...
    sin(2.*pi.*(X + r.*tau).*kn).*real(S0.*exp(-(X + r.*tau)).*f1) +.... 
    2.*pi.*kn.*sin(2.*pi.*(X + r.*tau).*kn).*(imag(S0.*exp(-(X + r.*tau)).*f1 - S0.*f2)) +...
    cos(2.*pi.*(X + r.*tau).*kn).*imag(S0.*exp(-(X + r.*tau)).*f1)));
C_convexity_SINC(i) = -0.5.*S0.*exp(-(X + r.*tau)) + 2./pi.*sum(1./(2.*n - 1).*...
    ((2.*pi.*kn).^2.*sin(2.*pi.*(X + r.*tau).*kn).*(real(S0.*exp(-(X + r.*tau)).*f1 - S0.*f2)) +...
    4.*pi.*kn.*cos(2.*pi.*(X + r.*tau).*kn).*real(S0.*exp(-(X + r.*tau)).*f1) -...
    sin(2.*pi.*(X + r.*tau).*kn).*real(S0.*exp(-(X + r.*tau)).*f1) +...
    (2.*pi.*kn).^2.*cos(2.*pi.*(X + r.*tau).*kn).*(imag(S0.*exp(-(X + r.*tau)).*f1 - S0.*f2)) -...
    4.*pi.*kn.*sin(2.*pi.*(X + r.*tau).*kn).*imag(S0.*exp(-(X + r.*tau)).*f1) -...
    cos(2.*pi.*(X + r.*tau).*kn).*imag(S0.*exp(-(X + r.*tau)).*f1)));
f11 = charac_T(tau,2.*pi.*kn,theta,kappa,sigma,rho,V);
f22 = charac_T(tau,2.*pi.*kn - 1i,theta,kappa,sigma,rho,V);
C_time_SINC(i) = 0.5.*r.*S0.*exp(-(X + r.*tau)) + 2./pi.*sum(1./(2.*n - 1).*...
    (-2.*pi.*kn.*r.*cos(2.*pi.*(X + r.*tau).*kn).*(real(S0.*exp(-(X + r.*tau)).*f1 - S0.*f2)) + ...
    r.*sin(2.*pi.*(X + r.*tau).*kn).*real(S0.*exp(-(X + r.*tau)).*f1) - ...
    sin(2.*pi.*(X + r.*tau).*kn).*(real(S0.*exp(-(X + r.*tau)).*f11 - S0.*f22))+.... 
    2.*pi.*kn.*r.*sin(2.*pi.*(X + r.*tau).*kn).*(imag(S0.*exp(-(X + r.*tau)).*f1 - S0.*f2)) +...
    r.*cos(2.*pi.*(X + r.*tau).*kn).*imag(S0.*exp(-(X + r.*tau)).*f1) - ...
    cos(2.*pi.*(X + r.*tau).*kn).*(imag(S0.*exp(-(X + r.*tau)).*f11 - S0.*f22))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ait-Sahalia method
sigma_00 = sqrt(V);
sigma_01 = rho.*sigma./(4.*sqrt(V));
sigma_02 = -(5.*rho.^2 - 2).*sigma.^2./(48.*V.^(3/2));
sigma_20 = (sigma.*(24.*rho.*(q-r) + sigma.*(rho.^2 - 4)) + V.*(12.*sigma.*rho - 24.*kappa) + ...
    24.*kappa.*theta)./(96.*sqrt(V));
sigma_21 = -sigma.*(16.*(2 - 5.*rho.^2).*(r - q).*sigma + rho.*(40.*theta.*kappa + ...
    3.*(3.*rho.^2 - 4).*sigma.^2 + V.*(4.*rho.*sigma - 8.*kappa)))./(384.*V.^(3/2));
sigma_22 = sigma.^2.*(sigma.*(1440.*rho.*(8.*rho.^2 - 5).*(q - r) + sigma.*(1883.*rho.^4 - ...
    2716.*rho.^2 + 608)) + 240.*theta.*kappa.*(23.*rho.^2 - 8) - ...
    120.*V.*(7.*rho.^2 - 2).*(2.*kappa - sigma.*rho))./(46080.*V.^(5/2));

varepsilon = sqrt(tau);
k = log(S0./K(i));
% IV second-order expansion
IV_chenxu(i) = sigma_00 - sigma_01.*k + sigma_02.*k.^2 + sigma_20.*varepsilon.^2 - ...
    sigma_21.*varepsilon.^2.*k + sigma_22.*varepsilon.^2.*k.^2;
% IV sensitivities based on second-order expansion
IV_skew_Chenxu(i) = - sigma_01 + sigma_02.*2.*k - sigma_21.*varepsilon.^2 + ...
    + sigma_22.*varepsilon.^2.*2.*k;
IV_convexity_Chenxu(i) = sigma_02.*2 + sigma_22.*varepsilon.^2.*2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 500; % the number of abscissas for calculating epsilon
N2 = 1999; % the number of abscissas for calculating boundary
sigma_down = 1e-3;sigma_up = 1; % naive boundary
epsilon0 = 1e-4; % naive limiting parameter
k1 = 1e-3; % adjusted coefficient related to naive boundary

% compute the IV bound
[sigma_min_dirac,sigma_max_dirac,epsilon] = Dirac_delta_bound(S0,K,r,q,tau,C_SINC,...
    N,N2,sigma_down,sigma_up,epsilon0,k1);

% compute the IV sensitivities based on the formula 19 and 20
k2 = 1e-3; % adjusted coefficient related to tighter boundary
epsilon_dirac = epsilon.*k2; % improved limiting parameter
N = N2;
X = log(S0./K);

% generate discrete points
i = (1:N+1)';
sigma_all = sigma_min_dirac + (i-1).*(sigma_max_dirac - sigma_min_dirac)./N;

% calculate variables related to IV sensitivities
% the vega of the BS model and its derivatives with respect to X(=log(S./K) and T
vega = (exp(-r.*tau-1./2.*(sqrt(tau).*(sigma_all)-(tau.*(-q+r+(sigma_all).^2./2)+...
    log(S0./K))./(sqrt(tau).*(sigma_all))).^2).*K.*(tau.*(-q+r+(sigma_all).^2./2)+...
    log(S0./K)))./(sqrt(2.*(pi)).*sqrt(tau).*(sigma_all).^2)-(exp(-q.*tau-(tau.*...
    (-q+r+(sigma_all).^2./2)+log(S0./K)).^2./(2.*tau.*(sigma_all).^2)).*S0.*...
    (-(sqrt(tau)./sqrt(2))+(tau.*(-q+r+(sigma_all).^2./2)+log(S0./K))./(sqrt(2).*...
    sqrt(tau).*(sigma_all).^2)))./sqrt((pi));
vega_X = -(1./(2.*sqrt(2.*(pi)).*sqrt(tau).*(sigma_all).^2)).*exp(1./8.*(-4.*((q+r).*tau+...
    X)-(4.*(-q.*tau+r.*tau+X).^2)./(tau.*(sigma_all).^2)-tau.*(sigma_all).^2)).*S0.*(2.*...
    X+tau.*(-2.*q+2.*r+(sigma_all).^2));
vega_X_X = (1./(4.*sqrt(2.*(pi)).*tau.^(3./2).*(sigma_all).^4)).*exp(1./8.*(-4.*((q+r).*tau+X)-...
    (4.*(-q.*tau+r.*tau+X).^2)./(tau.*(sigma_all).^2)-tau.*(sigma_all).^2)).*S0.*(4.*(-q.*tau+...
    r.*tau+X).^2+4.*tau.*(-1-q.*tau+r.*tau+X).*(sigma_all).^2+tau.^2.*(sigma_all).^4);
vega_T = (exp(1./8.*(-4.*((q+r).*tau+X)-(4.*(-q.*tau+r.*tau+X).^2)./(tau.*sigma_all.^2)-...
    tau.*sigma_all.^2)).*S0.*(4.*X.^2+tau.*(-4.*(q-r).^2.*tau-4.*(-1+(q+r).*tau).*sigma_all.^2-...
    tau.*sigma_all.^4)))./(8.*sqrt(2.*pi).*tau.^(3./2).*sigma_all.^2);

% the derivatives of BS model option price with respect to X and T
C_X = 1./2.*exp(-r.*tau-X).*S0.*erfc((2.*q.*tau-2.*(r.*tau+X)+tau.*(sigma_all).^2)./...
    (2.*sqrt(2).*sqrt(tau).*(sigma_all)));
C_X_X = (exp(1./8.*(-4.*((q+r).*tau+X)-(4.*(-q.*tau+r.*tau+X).^2)./(tau.*(sigma_all).^2)-...
    tau.*(sigma_all).^2)).*S0)./(sqrt(2.*(pi)).*sqrt(tau).*(sigma_all))-1./2.*exp(-r.*tau-X).*...
    S0.*erfc((2.*q.*tau-2.*(r.*tau+X)+tau.*(sigma_all).^2)./(2.*sqrt(2).*sqrt(tau).*(sigma_all)));
C_T = 1./4.*S0.*((exp(1./8.*(-4.*((q+r).*tau+X)-(4.*(-q.*tau+r.*tau+X).^2)./(tau.*sigma_all.^2)-...
    tau.*sigma_all.^2)).*sqrt(2./pi).*sigma_all)./sqrt(tau)+2.*exp(-r.*tau-X).*r.*erfc((2.*q.*tau-...
    2.*(r.*tau+X)+tau.*sigma_all.^2)./(2.*sqrt(2).*sqrt(tau).*sigma_all))-2.*exp(-q.*tau).*q.*...
    erfc(-((2.*X+tau.*(-2.*q+2.*r+sigma_all.^2))./(2.*sqrt(2).*sqrt(tau).*sigma_all))));

% intermediate variable
varsigma = vega.*(-(European_call(S0,K,r,tau,sigma_all,q) - C_SINC)./(2.*epsilon_dirac)).*...
    (C_X - C_skew_SINC) + vega_X;

% compute the IV and its sensitivities with respect to X and T
f = sigma_all.*exp(-(European_call(S0,K,r,tau,sigma_all,q) - C_SINC).^2./...
    (4.*epsilon_dirac)).*(vega);
IV_Dirac_Delta = (1./(2.*sqrt(pi.*epsilon_dirac))).*(sigma_max_dirac - sigma_min_dirac)/...
    N.*(0.5.*(f(1,:) + f(end,:)) + sum(f(2:end-1,:)));

f1 = sigma_all.*exp(-(European_call(S0,K,r,tau,sigma_all,q) - C_SINC).^2./...
    (4.*epsilon_dirac)).*...
    (vega.*(-(European_call(S0,K,r,tau,sigma_all,q) - C_SINC)./(2.*epsilon_dirac)).*...
    (C_X - C_skew_SINC) + vega_X);
IV_skew_Dirac_Delta = (1./(2.*sqrt(pi.*epsilon_dirac))).*(sigma_max_dirac - sigma_min_dirac)/...
    N.*(0.5.*(f1(1,:) + f1(end,:)) + sum(f1(2:end-1,:)));

f2 = sigma_all.*exp(-(European_call(S0,K,r,tau,sigma_all,q) - C_SINC).^2./...
    (4.*epsilon_dirac)).*...
    ((-(European_call(S0,K,r,tau,sigma_all,q) - C_SINC)./(2.*epsilon_dirac)).*...
    ((C_X - C_skew_SINC).*(varsigma + vega_X) + vega.*(C_X_X - C_convexity_SINC)) + ...
    (-0.5./epsilon_dirac).*vega.*(C_X - C_skew_SINC).^2 + vega_X_X);
IV_convexity_Dirac_Delta = (1./(2.*sqrt(pi.*epsilon_dirac))).*(sigma_max_dirac - sigma_min_dirac)/...
    N.*(0.5.*(f2(1,:) + f2(end,:)) + sum(f2(2:end-1,:)));

f3 = sigma_all.*exp(-(European_call(S0,K,r,tau,sigma_all,q) - C_SINC).^2./...
    (4.*epsilon_dirac)).*...
    (vega.*(-(European_call(S0,K,r,tau,sigma_all,q) - C_SINC)./(2.*epsilon_dirac)).*...
    (C_T - C_time_SINC) + vega_T);
IV_time_Dirac_Delta = (1./(2.*sqrt(pi.*epsilon_dirac))).*(sigma_max_dirac - sigma_min_dirac)/...
    N.*(0.5.*(f3(1,:) + f3(end,:)) + sum(f3(2:end-1,:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the local volatility
% calculate the benchmark of local volatility based on Dupire formula 24
C_T = C_time_SINC;
C_K = C_skew_SINC.*(-1./S0.*exp(X));
C_K_K = (C_convexity_SINC + C_skew_SINC).*(1./S0.*exp(X)).^2;
Local_vol_Dupire = (C_T + q.*C_SINC + (r - q).*K.*C_K)./(0.5.*K.^2.*C_K_K);

% calculate the local volatility based on proposed formula 21
IV_K = IV_skew_Dirac_Delta.*(-1./K);
IV_K_K = (IV_convexity_Dirac_Delta + IV_skew_Dirac_Delta).*(1./K.^2);
y = log(K./(S0.*exp((r-q).*tau)));
Local_vol_Dirac_Delta = (IV_Dirac_Delta.^2 + 2.*IV_Dirac_Delta.*tau.*...
    (IV_time_Dirac_Delta + r.*K.*IV_K))./...
    ((1 - K.*y./IV_Dirac_Delta.*IV_K).^2 + K.*IV_Dirac_Delta.*tau.*...
    (IV_K - 1./4.*K.*IV_Dirac_Delta.*tau.*IV_K.^2 + K.*IV_K_K));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch tau
    case 1/12
figure(10)
    case 2
 figure(11)
end

subplot(1,2,1)
plot(X,abs(IV_skew_finite_difference - IV_skew_Dirac_Delta),'r*-','linewidth',1)
hold on
plot(X,abs(IV_skew_finite_difference - IV_skew_Chenxu),'b+-','linewidth',1)
legend('Dirac Delta method','Aït-Sahalia''s method');
xlabel('$X$','Interpreter','latex')
ylabel('Absolute error of $\partial \sigma /\partial X$','Interpreter','latex')
set(gca,'YScale','log');

subplot(1,2,2)
plot(X,abs(IV_convexity_finite_difference - IV_convexity_Dirac_Delta),'r*-','linewidth',1)
hold on
plot(X,abs(IV_convexity_finite_difference - IV_convexity_Chenxu),'b+-','linewidth',1)
legend('Dirac Delta method','Aït-Sahalia''s method');
xlabel('$X$','Interpreter','latex')
ylabel('Absolute error of ${\partial ^2}\sigma /\partial {X^2}$','Interpreter','latex')
set(gca,'YScale','log');

switch tau
    case 1/12
sgtitle('$T=1/12$','Interpreter','latex')
    case 2
sgtitle('$T=1/2$','Interpreter','latex')
end

figure(12)
switch tau
    case 1/12
subplot(1,2,1)
plot(X,Local_vol_Dupire,'bo-','linewidth',1)
hold on
plot(X,Local_vol_Dirac_Delta,'r*-','linewidth',1)
legend('Dupire formula','Dirac Delta method');
xlabel('$X$','Interpreter','latex')
ylabel('Local volatility')
title('$T=1/12$','Interpreter','latex')

    case 1/2
subplot(1,2,2)
plot(X,Local_vol_Dupire,'bo-','linewidth',1)
hold on
plot(X,Local_vol_Dirac_Delta,'r*-','linewidth',1)
legend('Dupire formula','Dirac Delta method');
xlabel('$X$','Interpreter','latex')
ylabel('Local volatility')
title('$T=1/2$','Interpreter','latex')
end










