%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Tighter bounds for implied volatility based on the Dirac delta family method
% Author: Zhenyu Cui,Yanchu Liu and Yuhang Yao
% Date: August 13, 2023
% Note: This document corresponds to secton 3.2 of the paper, with figure(8)-(10) and table 2
% MATLAB version: Matlab2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;
format long 
% initial parameter
S0 = 100;r = 0.03;q = 0;

type = 3;%1,2,3
switch type
    case 1  % figure(8) and table 2
% generate the three dimension dataset
num_K = 39;
num_tau = 39;
num_sigma = 39;
K0 = 105:(800-105)/num_K:800;
tau0 = 0.01:(2-0.01)/num_tau:2;
sigma0 = 0.01:(0.99-0.01)/num_sigma:0.99;
[x,y,z] = ndgrid(K0,tau0,sigma0);
K = reshape(x,1,(num_K+1)*(num_tau+1)*(num_sigma+1));
tau = reshape(y,1,(num_K+1)*(num_tau+1)*(num_sigma+1));
sigma = reshape(z,1,(num_K+1)*(num_tau+1)*(num_sigma+1));

% filter option prices that are too small to match market prices
C_real = European_call(S0,K,r,tau,sigma,q);
temp0 = find(C_real<1e-20);
C_real(temp0) = [];
tau(temp0) = [];
K(temp0) = [];
sigma(temp0) = [];

    case 2 % figure(9)
sigma = 0.2;
num_K = 19;
num_tau = 19;
tau0 = 0.1:(2-0.1)/num_tau:2;
K0 = 105:(180-105)/num_K:180;
[x,y] = ndgrid(K0,tau0);
K = reshape(x,1,(num_K+1)*(num_tau+1));
tau = reshape(y,1,(num_K+1)*(num_tau+1));

C_real = European_call(S0,K,r,tau,sigma,q);
temp0 = find(C_real<1e-20);
C_real(temp0) = [];
tau(temp0) = [];
K(temp0) = [];

    case 3 % figure(10)
sigma = 0.8;
num_K = 19;
num_tau = 19;
tau0 = 0.1:(2-0.1)/num_tau:2;
K0 = 105:(800-105)/num_K:800;
[x,y] = ndgrid(K0,tau0);
K = reshape(x,1,(num_K+1)*(num_tau+1));
tau = reshape(y,1,(num_K+1)*(num_tau+1));

C_real = European_call(S0,K,r,tau,sigma,q);
temp0 = find(C_real<1e-20);
C_real(temp0) = [];
tau(temp0) = [];
K(temp0) = [];   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
N = 500; % the number of abscissas for calculating epsilon
N2 = 1999; % the number of abscissas for calculating boundary
sigma_down = 1e-3;sigma_up = 1; % naive boundary
epsilon0 = 1e-4; % naive limiting parameter
k1 = 1e-3; % adjusted coefficient related to boundary

% compute tighter bound
[sigma_min_dirac,sigma_max_dirac,epsilon] = Dirac_delta_bound(S0,K,r,q,tau,C_real,...
    N,N2,sigma_down,sigma_up,epsilon0,k1);

% midpoint of the upper and lower boundaries is chosen as the initial point of Newton-Raphson method
N_iterative = 101;
IV_DBNR = zeros(N_iterative,length(C_real));% Pre allocate memory to related variables
IV_DBNR(1,:) = (sigma_min_dirac + sigma_max_dirac)./2;

for i = 2:N_iterative
    sigma_Dirac_bound = IV_DBNR(i-1,:);
    C_Dirac_bound = European_call(S0,K,r,tau,sigma_Dirac_bound,q);
    vega_Dirac_bound = (exp((-((tau.^2.*(4.*(q-r).^2+4.*(q+r).*(sigma_Dirac_bound).^2+(sigma_Dirac_bound).^4)+...
        4.*log(S0./K).^2)./(8.*tau.*(sigma_Dirac_bound).^2)))).*K.*(S0./K).^(1./2+(q-r)./(sigma_Dirac_bound).^2).*...
        sqrt(tau))./sqrt(2.*(pi));
    % first-order form of Newton-Raphson method
    IV_DBNR(i,:) = sigma_Dirac_bound - (C_Dirac_bound - C_real)./vega_Dirac_bound;
end

% absolute error and its statistics
t2_a = abs(IV_DBNR - sigma);
t2_a_statistics = [mean(t2_a(end,:)),std(t2_a(end,:)),max(t2_a(end,:)),min(t2_a(end,:))];
% average error of each iteration step
t2_b = mean(abs(IV_DBNR - sigma),2);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choi method
tic
C_real_2 = exp((r - q).*tau).*C_real;
c = (C_real_2 - max(S0.*exp((r-q).*tau) - K,0))./min(K,S0.*exp((r-q).*tau));
k = abs(log(S0.*exp((r-q).*tau)./K));
temp1 = min((1+c)./2,c + exp(k).*normcdf(-sqrt(2.*k)));
sigma_max_Choi = (norminv(temp1) - norminv((temp1 - c)./exp(k)));
d1 = -k./sigma_max_Choi + sigma_max_Choi./2;
d2 = -k./sigma_max_Choi - sigma_max_Choi./2;
temp2 =  1 - exp(k).*normcdf(d2)./normcdf(d1);
temp3 = norminv(c./temp2);
sigma_min_Choi = (temp3 + sqrt(temp3.^2 + 2.*k))./sqrt(tau);
temp0 = find(isnan(sigma_min_Choi)==1);
sigma_min_Choi(temp0) = 0;
sigma_max_Choi = (norminv(temp1) - norminv((temp1 - c)./exp(k)))./sqrt(tau);

N_iterative = 101;
IV_Choi = zeros(N_iterative,length(K));
IV_Choi(1,:) = sigma_min_Choi;
for i = 2:N_iterative
    sigma_normal = IV_Choi(i-1,:).*sqrt(tau);
    d1 = -k./sigma_normal + sigma_normal./2;
    d2 = -k./sigma_normal - sigma_normal./2;
    CV = (normcdf(d1,0,1) - exp(k).*normcdf(d2,0,1))./normpdf(d1,0,1);
    IV_Choi(i,:) = 1./sqrt(tau).*...
        (sigma_normal + (d1.^2/2 - log(CV) + log(c.*sqrt(2.*pi))).*CV);
end
toc

% absolute error and its statistics
t3_a = abs(IV_Choi - sigma);
t3_a_statistics = [mean(t3_a(end,:)),std(t3_a(end,:)),max(t3_a(end,:)),min(t3_a(end,:))];
% average error of each iteration step
t3_b = mean(abs(IV_Choi - sigma),2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inversion formula 2
type2 = 1;%1,2
switch type2
    case 1 % naive bound and epsilon
k2 = 1e-4;
epsilon_inverse = epsilon.*k2;
N = N2;
    case 2 % improved bound and epsilon
epsilon_inverse = 1e-4;
N = 3999;
sigma_min_dirac = 1e-3;sigma_max_dirac = 1;
end

tic
% generate discrete points
i = (1:N+1)';
sigma_all = sigma_min_dirac + (i-1).*(sigma_max_dirac - sigma_min_dirac)./N;

% compute IV according to formula 2
vega = (exp(-r.*tau-1./2.*(sqrt(tau).*(sigma)-(tau.*(-q+r+(sigma).^2./2)+...
    log(S0./K))./(sqrt(tau).*(sigma))).^2).*K.*(tau.*(-q+r+(sigma).^2./2)+...
    log(S0./K)))./(sqrt(2.*(pi)).*sqrt(tau).*(sigma).^2)-(exp(-q.*tau-(tau.*...
    (-q+r+(sigma).^2./2)+log(S0./K)).^2./(2.*tau.*(sigma).^2)).*S0.*...
    (-(sqrt(tau)./sqrt(2))+(tau.*(-q+r+(sigma).^2./2)+log(S0./K))./(sqrt(2).*...
    sqrt(tau).*(sigma).^2)))./sqrt((pi));
f = sigma_all.*exp(-(European_call(S0,K,r,tau,sigma_all,q) - C_real).^2./...
    (4.*epsilon_inverse)).*...
    (vega);
IV_inverse = (1./(2.*sqrt(pi.*epsilon_inverse))).*(sigma_max_dirac - sigma_min_dirac)/N.*...
    (0.5.*(f(1,:) + f(end,:)) + sum(f(2:end,:)));
toc

% absolute error and its statistics
t4_a = abs(IV_inverse - sigma);
t4_a_statistics = [mean(t4_a),std(t4_a),max(t4_a),min(t4_a)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab algorithm
tic
IV_brent = blsimpv(S0,K,r,tau,C_real,'Method','search'); % Dekker-Brent method
toc

tic
IV_jackel = blsimpv(S0,K,r,tau,C_real,'Method','jackel2016');% Jackel’s method
toc

% absolute error and its statistics
t5_a = abs(IV_brent - sigma);
t5_a_statistics = [mean(t5_a),std(t5_a),max(t5_a),min(t5_a)];

t6_a = abs(IV_jackel - sigma);
t6_a_statistics = [mean(t6_a),std(t6_a),max(t6_a),min(t6_a)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch type
    case 1 
figure(8)
num_iterative = 0:N_iterative-1;
subplot(1,2,1)
plot(num_iterative(1:11),t2_b(1:11),'r-*','markersize',4,'linewidth',1)
hold on
plot(num_iterative(1:11),t3_b(1:11),'b-o','markersize',4,'linewidth',1)
axis([0,10,1e-20,1])
legend('DBNR method','Choi method')
xlabel('Number of Iterations')
ylabel('Average AEIV')
set(gca,'YScale','log');

subplot(1,2,2)
plot(num_iterative(end-9:end),t2_b(end-9:end),'r-*','markersize',4,'linewidth',1)
hold on
plot(num_iterative(end-9:end),t3_b(end-9:end),'b-o','markersize',4,'linewidth',1)
axis([num_iterative(end-9),num_iterative(end),3e-16,5.0e-16])
legend('DBNR method','Choi method')
xlabel('Number of Iterations')
ylabel('Average AEIV')

    otherwise % figure(9)-(10)
% Increased AEIV of Dekker-Brent method relative to DBNR method(of which error
% of each output IV in the final iteration)
accuracy_up_brent = t5_a - t2_a(end,:);
figure(9)
subplot(1,2,1)
K_all = num2cell(round(K0,0));
tau_all = num2cell(round(tau0,2));
accuracy_up_brent_reshape = reshape(accuracy_up_brent,num_K+1,num_tau+1);
h = heatmap(tau_all,K_all,accuracy_up_brent_reshape,...
    'ColorLimits',[-max(accuracy_up_brent) max(accuracy_up_brent)]);
colormap jet
ylabel('\itK');
xlabel('\itT');
title(['Increased AEIV of Dekker-Brent method(\sigma = ',num2str(sigma),')'])

% Increased AEIV of Jackel’s method relative to DBNR method(of which error
% of each output IV in the final iteration)
accuracy_up_jackel = t6_a - t2_a(end,:);
temp1 = find(accuracy_up_jackel>1e-1);
accuracy_up_jackel((temp1)) = nan;
subplot(1,2,2)
accuracy_up_jackel_reshape = reshape(accuracy_up_jackel,num_K+1,num_tau+1);
h = heatmap(tau_all,K_all,accuracy_up_jackel_reshape,...
    'ColorLimits',[-max(accuracy_up_jackel) max(accuracy_up_jackel)]);
colormap jet
title(['Increased AEIV of Jackel method(\sigma = ',num2str(sigma),')'])
ylabel('\itK');
xlabel('\itT');
h.MissingDataLabel='0.1';
end

