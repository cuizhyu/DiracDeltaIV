%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Tighter bounds for implied volatility based on the Dirac delta family method
% Author: Zhenyu Cui,Yanchu Liu and Yuhang Yao
% Date: August 13, 2023
% Note: This document corresponds to secton 3.1 of the paper, with figure (4)-(7)
% MATLAB version: Matlab2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;
type = 1;
switch type
    case 1 % figure 4-5
% initial parameter
S0 = 100;r = 0.03;q = 0;tau = 1;
N = 500; % the number of abscissas for calculating epsilon
N2 = 2000; % the number of abscissas for calculating boundary
sigma_down = 1e-3;sigma_up = 1; % naive boundary
epsilon0 = 1e-4; % naive limiting parameter
k1 = 1e-3; % adjusted coefficient related to boundary
K_all = (105:(800-105)/100:800);% generate discrete strike corresponding OTM options
sigma_real_all = 0.2.*ones(1,length(K_all));

    case 2 % figure 6
S0 = 100;r = 0;q = 0;tau = 16;
N = 500; 
N2 = 2000; 
sigma_down = 1e-3;sigma_up = 1; 
epsilon0 = 1e-4;
k1 = 1e-3;
sigma_real_all = 0.05:0.02:0.95;
k = 0.2;
K_all = exp(k).*S0.*exp((r-q).*tau).*ones(1,length(sigma_real_all));

    case 3 % figure 7
S0 = 1;r = 0;q = 0;tau = 9;
N = 500; 
N2 = 4000; 
sigma_down = 1e-3;sigma_up = 1; 
epsilon0 = 1e-10;
k1 = 1e-4;
K_all = (0.05:(5-0.01)/100:5);
sigma_real_all = 0.5.*ones(1,length(K_all));
end

% Pre allocate memory to related variables to improve calculation speed
sigma_min_dirac = zeros(1,length(K_all));sigma_max_dirac = zeros(1,length(K_all));
sigma_min_Tehranchi = zeros(1,length(K_all));sigma_max_Tehranchi = zeros(1,length(K_all));
sigma_min_Choi = zeros(1,length(K_all));sigma_max_Choi = zeros(1,length(K_all));
c_all = zeros(1,length(K_all));

for i = 1:length(K_all)
% Option prices with different stirke are calculated separately here for a more intuitive display. 
% In the next section, multiple options will be calculated in parallel
K = K_all(i);
sigma_real = sigma_real_all(i);
C_real = European_call(S0,K,r,tau,sigma_real,q);

% generate discrete points
sigma_all = (sigma_down:(sigma_up - sigma_down)/N:sigma_up)';
C_all = European_call(S0,K,r,tau,sigma_all,q);

% locate the positions to compute epsilon according to formula 11
temp0 = sum(C_all<C_real);
temp1 = [temp0;temp0+1];

% compute the vega and its derivative
vega1 = (exp((-((tau.^2.*(4.*(q-r).^2+4.*(q+r).*(sigma_all(temp1)).^2+(sigma_all(temp1)).^4)+...
    4.*log(S0./K).^2)./(8.*tau.*(sigma_all(temp1)).^2)))).*K.*(S0./K).^(1./2+(q-r)./(sigma_all(temp1)).^2).*...
    sqrt(tau))./sqrt(2.*(pi));
vvega1 = -(1./(4.*sqrt(2.*(pi)).*sqrt(tau).*(sigma_all(temp1)).^3)).*exp((-((tau.^2.*...
    (4.*(q-r).^2+4.*(q+r).*(sigma_all(temp1)).^2+(sigma_all(temp1)).^4)+4.*log(S0./K).^2)./(8.*tau.*...
    (sigma_all(temp1)).^2)))).*K.*(S0./K).^(1./2+(q-r)./(sigma_all(temp1)).^2).*(tau.*(2.*q-2.*r+(sigma_all(temp1)).^2)-...
    2.*log(S0./K)).*(tau.*(-2.*q+2.*r+(sigma_all(temp1)).^2)+2.*log(S0./K));
vvvega1 = (1./(16.*sqrt(2.*(pi)).*tau.^(3./2).*(sigma_all(temp1)).^6)).*exp(-((tau.^2.*(4.*(q-r).^2+4.*...
    (q+r).*(sigma_all(temp1)).^2+(sigma_all(temp1)).^4)-4.*log(S0./K).^2)./(8.*tau.*(sigma_all(temp1)).^2))).*K.*(S0./K).^(1./2+...
    q./(sigma_all(temp1)).^2-(r.*tau+log(S0./K))./(tau.*(sigma_all(temp1)).^2)).*(tau.^3.*(16.*(q-r).^4.*tau-48.*(q-...
    r).^2.*(sigma_all(temp1)).^2-8.*(q-r).^2.*tau.*(sigma_all(temp1)).^4-4.*(sigma_all(temp1)).^6+tau.*(sigma_all(temp1)).^8)-8.*(2.*(q-r).*...
    tau-log(S0./K)).*log(S0./K).*(tau.*(4.*(q-r).^2.*tau-6.*(sigma_all(temp1)).^2-tau.*(sigma_all(temp1)).^4)+2.*...
    log(S0./K).*(2.*(-q+r).*tau+log(S0./K))));

% compute the reasonable epsilon according to formula 13
epsilon = min(min(abs(sigma_all(temp1).*vega1.^3./(4.*vvega1 + ...
    2.*sigma_all(temp1).*vvvega1)).*k1),epsilon0);

% regenerate more discrete points to compute the boundary
sigma_all = (sigma_down:(sigma_up - sigma_down)/N2:sigma_up)';
C_all = European_call(S0,K,r,tau,sigma_all,q);

% recompute the vega and its derivative
vega_all = (exp((-((tau.^2.*(4.*(q-r).^2+4.*(q+r).*(sigma_all).^2+(sigma_all).^4)+...
    4.*log(S0./K).^2)./(8.*tau.*(sigma_all).^2)))).*K.*(S0./K).^(1./2+(q-r)./(sigma_all).^2).*...
    sqrt(tau))./sqrt(2.*(pi));
vvega_all = -(1./(4.*sqrt(2.*(pi)).*sqrt(tau).*(sigma_all).^3)).*exp((-((tau.^2.*...
    (4.*(q-r).^2+4.*(q+r).*(sigma_all).^2+(sigma_all).^4)+4.*log(S0./K).^2)./(8.*tau.*...
    (sigma_all).^2)))).*K.*(S0./K).^(1./2+(q-r)./(sigma_all).^2).*(tau.*(2.*q-2.*r+(sigma_all).^2)-...
    2.*log(S0./K)).*(tau.*(-2.*q+2.*r+(sigma_all).^2)+2.*log(S0./K));

% compute the intergrand
f1 = 1./sqrt(2.*pi.*epsilon).*exp(-(C_all - C_real).^2./(4.*epsilon)).*vega_all + ...
    sigma_all./sqrt(2.*pi.*epsilon).*exp(-(C_all - C_real).^2./(4.*epsilon)).*vega_all.*...
    (-2.*(C_all - C_real)./(4.*epsilon).*vega_all) + ...
    sigma_all./sqrt(2.*pi.*epsilon).*exp(-(C_all - C_real).^2./(4.*epsilon)).*vvega_all;

% compute Dirac delta bound according to formula 14
% the abscissas corresponding to the maximum and minimum values
temp4 = (f1==max(f1));
sigma_min_dirac(i)  = sigma_all(temp4);
temp4 = (f1==min(f1));
sigma_max_dirac(i)  = sigma_all(temp4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute Tehranchi bound
c = C_real./(exp(-q.*tau).*S0);
k = log(K./(S0.*exp((r-q).*tau)));
sigma_min_Tehranchi(i) = (k>=0).*(-2./sqrt(tau).*norminv((1-c)./2)) + ...
    (1 - (k>=0)).*(-2./sqrt(tau).*norminv((1-c)./(2.*exp(k))));
sigma_max_Tehranchi(i) = -2./sqrt(tau).*norminv((1-c)./(1+exp(k)));
c_all(i) = c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute Choi bound
C_Choi = exp((r - q).*tau).*C_real;
c = (C_Choi - max(S0.*exp((r-q).*tau) - K,0))./min(K,S0.*exp((r-q).*tau));
k = abs(log(S0.*exp((r-q).*tau)./K));
temp1 = min((1+c)./2,c + exp(k).*normcdf(-sqrt(2.*k)));
sigma_max_Choi(i) = (norminv(temp1) - norminv((temp1 - c)./exp(k)));
d1 = -k./sigma_max_Choi(i) + sigma_max_Choi(i)./2;
d2 = -k./sigma_max_Choi(i) - sigma_max_Choi(i)./2;
temp2 =  1 - exp(k).*normcdf(d2)./normcdf(d1);
temp3 = norminv(c./temp2);
sigma_min_Choi(i) = (temp3 + sqrt(temp3.^2 + 2.*k))./sqrt(tau);
sigma_max_Choi(i) = (norminv(temp1) - norminv((temp1 - c)./exp(k)))./sqrt(tau);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch type
    case 1 % figure 4-5
figure(4)
subplot(1,3,1)
plot(K_all,ones(1,length(K_all)).*sigma_real,'k','linewidth',1)
hold on
plot(K_all,sigma_min_dirac,'r-.','linewidth',1)
hold on
plot(K_all,sigma_max_dirac,'r--','linewidth',1)
xlabel('$K$','Interpreter','latex');
ylabel('$\sigma $','Interpreter','latex');
legend('benchmark','${a^{Dirac}}$','${b^{Dirac}}$','Interpreter','latex');
title('Dirac Delta Bound')
axis([-inf inf 0.18 0.22])

subplot(1,3,2)
plot(K_all,ones(1,length(K_all)).*sigma_real,'k','linewidth',1)
hold on
plot(K_all,sigma_min_Tehranchi,'b-.','linewidth',1)
hold on
plot(K_all,sigma_max_Tehranchi,'b--','linewidth',1)
xlabel('$K$','Interpreter','latex');
ylabel('$\sigma $','Interpreter','latex');
legend('benchmark','${a^{Teh}}$','${b^{Teh}}$','Interpreter','latex');
title('Tehranchi Bound')
axis([-inf inf -inf inf])

subplot(1,3,3)
plot(K_all,ones(1,length(K_all)).*sigma_real,'k','linewidth',1)
hold on
plot(K_all,sigma_min_Choi,'m-.','linewidth',1)
hold on
plot(K_all,sigma_max_Choi,'m--','linewidth',1)
axis([-inf inf -inf inf])
xlabel('$K$','Interpreter','latex');
ylabel('$\sigma $','Interpreter','latex');
legend('benchmark','${a^{Choi}}$','${b^{Choi}}$','Interpreter','latex');
title('Choi Bound')
axis([-inf inf -inf inf])

figure(5)
plot(K_all,ones(1,length(K_all)).*sigma_real,'k','linewidth',1)
hold on
plot(K_all,sigma_min_dirac,'r-.','linewidth',1)
hold on
plot(K_all,sigma_min_Choi,'m-.','linewidth',1)
hold on
axis([-inf inf -inf 0.208])
xlabel('$K$','Interpreter','latex');
ylabel('$\sigma $','Interpreter','latex');
legend('benchmark','${a^{Dirac}}$','${a^{Choi}}$','Interpreter','latex');

    case 2 % figure 6
figure(6)
subplot(2,2,1)
plot(c_all,sqrt(tau).*sigma_real_all,'k','linewidth',1)
hold on
plot(c_all,sqrt(tau).*sigma_min_dirac,'r-.','linewidth',1)
hold on
plot(c_all,sqrt(tau).*sigma_max_dirac,'r--','linewidth',1)
legend('benchmark','${a^{Dirac}}\sqrt T$','${b^{Dirac}}\sqrt T$','Interpreter','latex');
xlabel('$c$','Interpreter','latex');
ylabel('$\sigma \sqrt T $','Interpreter','latex');
axis([-inf inf 0 4])
title('Dirac Delta Bound')

subplot(2,2,2)
plot(c_all,sqrt(tau).*sigma_real_all,'k','linewidth',1)
hold on
plot(c_all,sqrt(tau).*sigma_min_Tehranchi,'b-.','linewidth',1)
hold on
plot(c_all,sqrt(tau).*sigma_max_Tehranchi,'b--','linewidth',1)
legend('benchmark','${a^{Teh}}\sqrt T$','${b^{Teh}}\sqrt T$','Interpreter','latex');
xlabel('$c$','Interpreter','latex');
ylabel('$\sigma \sqrt T $','Interpreter','latex');
axis([-inf inf 0 4])
title('Tehranchi Bound')

subplot(2,2,3)
plot(c_all,sqrt(tau).*(sigma_min_dirac - sigma_real_all),'r','linewidth',1)
hold on
plot(c_all,sqrt(tau).*(sigma_min_Tehranchi - sigma_real_all),'b-.','linewidth',1)
legend('${a^{Dirac}}\sqrt T$','${a^{Teh}}\sqrt T$','Interpreter','latex');
xlabel('$c$','Interpreter','latex');
ylabel('$\Delta \sigma \sqrt T $','Interpreter','latex');
axis([-inf inf -0.2 0.01])
title('Lower Bound')

subplot(2,2,4)
plot(c_all,sqrt(tau).*(sigma_max_dirac - sigma_real_all),'r','linewidth',1)
hold on
plot(c_all,sqrt(tau).*(sigma_max_Tehranchi - sigma_real_all),'b--','linewidth',1)
legend('${b^{Dirac}}\sqrt T$','${b^{Teh}}\sqrt T$','Interpreter','latex');
xlabel('$c$','Interpreter','latex');
ylabel('$\Delta \sigma \sqrt T $','Interpreter','latex');
axis([-inf inf -0.005 0.1])
title('Upper Bound')


    case 3 % figure 7
subplot(2,2,1)
plot(K_all,sqrt(tau).*sigma_real.*ones(length(K_all),1),'k','linewidth',1)
hold on
plot(K_all,sqrt(tau).*sigma_min_dirac,'r-.','linewidth',1)
hold on
plot(K_all,sqrt(tau).*sigma_max_dirac,'r--','linewidth',1)
legend('benchmark','${a^{Dirac}}\sqrt T$','${b^{Dirac}}\sqrt T$','Interpreter','latex');
xlabel('$K$','Interpreter','latex');
ylabel('$\sigma \sqrt T $','Interpreter','latex');
axis([-inf inf -inf inf])
title('Dirac Delta Bound')

subplot(2,2,2)
plot(K_all,sqrt(tau).*sigma_real.*ones(length(K_all),1),'k','linewidth',1)
hold on
plot(K_all,sqrt(tau).*sigma_min_Choi,'b-.','linewidth',1)
hold on
plot(K_all,sqrt(tau).*sigma_max_Choi,'b--','linewidth',1)
legend('benchmark','${a^{Choi}}\sqrt T$','${b^{Choi}}\sqrt T$','Interpreter','latex');
xlabel('$K$','Interpreter','latex');
ylabel('$\sigma \sqrt T $','Interpreter','latex');
axis([-inf inf -inf inf])
title('Choi Bound')

subplot(2,2,3)
plot(K_all,sqrt(tau).*(sigma_min_dirac - sigma_real),'r','linewidth',1)
hold on
plot(K_all,sqrt(tau).*(sigma_min_Choi - sigma_real),'b-.','linewidth',1)
legend('${a^{Dirac}}\sqrt T$','${a^{Choi}}\sqrt T$','Interpreter','latex');
xlabel('$K$','Interpreter','latex');
ylabel('$\Delta \sigma \sqrt T $','Interpreter','latex');
axis([-inf inf -0.02 0.01])
title('Lower Bound')

subplot(2,2,4)
plot(K_all,sqrt(tau).*(sigma_max_dirac - sigma_real),'r','linewidth',1)
hold on
plot(K_all,sqrt(tau).*(sigma_max_Choi - sigma_real),'b--','linewidth',1)
legend('${b^{Dirac}}\sqrt T$','${b^{Choi}}\sqrt T$','Interpreter','latex');
xlabel('$K$','Interpreter','latex');
ylabel('$\Delta \sigma \sqrt T $','Interpreter','latex');
axis([-inf inf -0.005 0.1])
title('Upper Bound')
end




