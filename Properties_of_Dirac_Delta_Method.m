%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Tighter bounds for implied volatility based on the Dirac delta family method
% Author: Zhenyu Cui,Yanchu Liu and Yuhang Yao
% Date: August 13, 2023
% Note: This document corresponds to chapter 2 of the paper, with figure (1)-(3)
% MATLAB version: Matlab2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;
% initial parameter
S0 = 100;K = 100;r = 0.03;q = 0;tau = 1;sigma_real = 0.2;
% generate real option price
C_real = European_call(S0,K,r,tau,sigma_real,q);

N = 2000; % the number of abscissas
sigma_down = 1e-3;sigma_up = 1; % naive boundary
epsilon = 1e-4; % naive limiting parameter
% generate discrete points
sigma_all = (sigma_down:(sigma_up - sigma_down)/N:sigma_up);
C_all = European_call(S0,K,r,tau,sigma_all,q);

% compute the vega and its derivative
vega_all = (exp((-((tau.^2.*(4.*(q-r).^2+4.*(q+r).*(sigma_all).^2+(sigma_all).^4)+...
    4.*log(S0./K).^2)./(8.*tau.*(sigma_all).^2)))).*K.*(S0./K).^(1./2+(q-r)./(sigma_all).^2).*...
    sqrt(tau))./sqrt(2.*(pi));
vvega_all = -(1./(4.*sqrt(2.*(pi)).*sqrt(tau).*(sigma_all).^3)).*exp((-((tau.^2.*...
    (4.*(q-r).^2+4.*(q+r).*(sigma_all).^2+(sigma_all).^4)+4.*log(S0./K).^2)./(8.*tau.*...
    (sigma_all).^2)))).*K.*(S0./K).^(1./2+(q-r)./(sigma_all).^2).*(tau.*(2.*q-2.*r+(sigma_all).^2)-...
    2.*log(S0./K)).*(tau.*(-2.*q+2.*r+(sigma_all).^2)+2.*log(S0./K));
vvvega_all = (1./(16.*sqrt(2.*(pi)).*tau.^(3./2).*(sigma_all).^6)).*exp(-((tau.^2.*(4.*(q-r).^2+4.*...
    (q+r).*(sigma_all).^2+(sigma_all).^4)-4.*log(S0./K).^2)./(8.*tau.*(sigma_all).^2))).*K.*(S0./K).^(1./2+...
    q./(sigma_all).^2-(r.*tau+log(S0./K))./(tau.*(sigma_all).^2)).*(tau.^3.*(16.*(q-r).^4.*tau-48.*(q-...
    r).^2.*(sigma_all).^2-8.*(q-r).^2.*tau.*(sigma_all).^4-4.*(sigma_all).^6+tau.*(sigma_all).^8)-8.*(2.*(q-r).*...
    tau-log(S0./K)).*log(S0./K).*(tau.*(4.*(q-r).^2.*tau-6.*(sigma_all).^2-tau.*(sigma_all).^4)+2.*...
    log(S0./K).*(2.*(-q+r).*tau+log(S0./K))));

% compute the intergrand of Dirac delta family method and its derivative
f =  sigma_all./sqrt(2.*pi.*epsilon).*exp(-(C_all - C_real).^2./(4.*epsilon)).*vega_all;
f1 = 1./sqrt(2.*pi.*epsilon).*exp(-(C_all - C_real).^2./(4.*epsilon)).*vega_all + ...
    sigma_all./sqrt(2.*pi.*epsilon).*exp(-(C_all - C_real).^2./(4.*epsilon)).*vega_all.*...
    (-2.*(C_all - C_real)./(4.*epsilon).*vega_all) + ...
    sigma_all./sqrt(2.*pi.*epsilon).*exp(-(C_all - C_real).^2./(4.*epsilon)).*vvega_all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch K
    case 100
p = 280:520;% Show the right tail of the curve
figure(1)
subplot(2,1,1)
plot(sigma_all,f,'linewidth',1)
xlabel('\sigma')
ylabel('G(\sigma)')
subplot(2,1,2)
plot(sigma_all(p),f(p),'linewidth',1)
axis([-inf,inf,-inf,inf])
xlabel('\sigma')
ylabel('G(\sigma)')

figure(2)
subplot(2,1,1)
plot(sigma_all,f1,'linewidth',1)
xlabel('\sigma')
ylabel('G^\prime(\sigma)')
subplot(2,1,2)
plot(sigma_all(p),f1(p),'linewidth',1)
xlabel('\sigma')
ylabel('G^\prime(\sigma)')
axis([-inf,inf,-inf,inf])

    case 150
p = 1:200;% Show the left tail of the curve
figure(3)
subplot(3,2,1)
plot(sigma_all,vega_all,'r','linewidth',1)
xlabel('\sigma')
ylabel('f^\prime_{BS}(\sigma)')
subplot(3,2,2)
plot(sigma_all(p),vega_all(p),'b','linewidth',1)
xlabel('\sigma')
ylabel('f^\prime_{BS}(\sigma)')
axis([-inf,inf,-inf,inf])
subplot(3,2,3)
plot(sigma_all,vvega_all,'r','linewidth',1)
xlabel('\sigma')
ylabel('f^\prime^\prime_{BS}(\sigma)')
subplot(3,2,4)
plot(sigma_all(p),vvega_all(p),'b','linewidth',1)
xlabel('\sigma')
ylabel('f^\prime^\prime_{BS}(\sigma)')
axis([-inf,inf,-inf,inf])
subplot(3,2,5)
plot(sigma_all,vvvega_all,'r','linewidth',1)
xlabel('\sigma')
ylabel('f^\prime^\prime^\prime_{BS}(\sigma)')
subplot(3,2,6)
plot(sigma_all(p),vvvega_all(p),'b','linewidth',1)
xlabel('\sigma')
ylabel('f^\prime^\prime^\prime_{BS}(\sigma)')
axis([-inf,inf,-inf,inf])
end






