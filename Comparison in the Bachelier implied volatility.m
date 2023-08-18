%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Tighter bounds for implied volatility based on the Dirac delta family method
% Author: Zhenyu Cui,Yanchu Liu and Yuhang Yao
% Date: August 13, 2023
% Note: This document corresponds to secton 3.5 of the paper, with table 5
% MATLAB version: Matlab2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;
format long 
% initial parameter for Bachelier model
F0 = 1;q = 0;

% generate the three dimension dataset
num_K = 39;
num_tau = 39;
num_sigma = 39;
K0 = 1.01:(10-1.01)/num_K:10;
tau0 = 0.01:(2-0.01)/num_tau:2;
sigma0 = 0.01:(0.99-0.01)/num_sigma:0.99;
[x,y,z] = ndgrid(K0,tau0,sigma0);
K = reshape(x,1,(num_K+1)*(num_tau+1)*(num_sigma+1));
tau = reshape(y,1,(num_K+1)*(num_tau+1)*(num_sigma+1));
sigma = reshape(z,1,(num_K+1)*(num_tau+1)*(num_sigma+1));

% filter option prices that are too small to match market prices
C_real = Bachelier_call(F0,K,tau,sigma);
temp0 = find(abs(C_real)<1e-20);
C_real(temp0) = [];
tau(temp0) = [];
K(temp0) = [];
sigma(temp0) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
N = 500; % the number of abscissas for calculating epsilon
N2 = 1999; % the number of abscissas for calculating boundary
sigma_down = 1e-3;sigma_up = 1; % naive boundary
epsilon0 = 1e-6; % naive limiting parameter
k1 = 1e-3; % adjusted coefficient related to boundary

% generate discrete points
sigma_all = (sigma_down:(sigma_up - sigma_down)/N2:sigma_up)';
C_all = Bachelier_call(F0,K,tau,sigma_all);

% locate the positions to compute epsilon according to formula 11
temp00 = 1:round(N2/N):N2;
temp0 = sum(C_all(temp00,:)<C_real);
temp1 = [temp0.*4;temp0.*4+1];

% calculate the derivative of option prices with respect to volatility under the Bachelier model
vega=(exp(-((F0-K).^2./(2.*tau.*sigma_all(temp1).^2))).*sqrt(tau))./sqrt(2.*(pi));
vvega=(exp(-((F0-K).^2./(2.*tau.*sigma_all(temp1).^2))).*(F0-K).^2)./(sqrt(2.*(pi)).*sqrt(tau).*sigma_all(temp1).^3);
vvvega=(exp(-((F0-K).^2./(2.*tau.*sigma_all(temp1).^2))).*((F0-K).^4-3.*(F0-K).^2.*tau.*sigma_all(temp1).^2))./...
    (sqrt(2.*(pi)).*tau.^(3./2).*sigma_all(temp1).^6);

% compute the reasonable epsilon according to formula 13
epsilon = min(min(abs(sigma_all(temp1).*vega.^3./(4.*vvega + ...
    2.*sigma_all(temp1).*vvvega)).*k1),epsilon0);

% recompute the vega and its derivative
vega_all=(exp(-((F0-K).^2./(2.*tau.*sigma_all.^2))).*sqrt(tau))./sqrt(2.*(pi));
vvega_all=(exp(-((F0-K).^2./(2.*tau.*sigma_all.^2))).*(F0-K).^2)./(sqrt(2.*(pi)).*sqrt(tau).*sigma_all.^3);

% compute the intergrand
f1 = 1./sqrt(2.*pi.*epsilon).*exp(-(C_all - C_real).^2./(4.*epsilon)).*vega_all + ...
    sigma_all./sqrt(2.*pi.*epsilon).*exp(-(C_all - C_real).^2./(4.*epsilon)).*vega_all.*...
    (-2.*(C_all - C_real)./(4.*epsilon).*vega_all) + ...
    sigma_all./sqrt(2.*pi.*epsilon).*exp(-(C_all - C_real).^2./(4.*epsilon)).*vvega_all;

% compute Dirac delta bound according to formula 14
% the abscissas corresponding to the maximum and minimum values
sigma_all_2 = repmat(sigma_all,1,length(K));
temp4 = (f1==max(f1));
sigma_min_dirac  = sigma_all_2(temp4)';
temp4 = (f1==min(f1));
sigma_max_dirac  = sigma_all_2(temp4)';

% midpoint of the upper and lower boundaries is chosen as the initial point of Newton-Raphson method
N_iterative = 101;
IV_DBNR = zeros(N_iterative,length(K));
IV_DBNR(1,:) = (sigma_min_dirac + sigma_max_dirac)./2;

for i = 2:N_iterative
    sigma_Dirac_bound = IV_DBNR(i-1,:);
    C_Dirac_bound = Bachelier_call(F0,K,tau,sigma_Dirac_bound);
    vega_Dirac_bound = (exp(-((F0-K).^2./(2.*tau.*sigma_Dirac_bound.^2))).*sqrt(tau))./sqrt(2.*(pi));
    IV_DBNR(i,:) = sigma_Dirac_bound - (C_Dirac_bound - C_real)./vega_Dirac_bound;
    
end
toc
% absolute error and its statistics
t1 = abs(IV_DBNR(end,:) - sigma);
t1_statistics = [mean(t1(1:end)),std(t1(1:end)),max(t1),min(t1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab algorithm
% Jackel’s method
tic
IV_Jackel = Jackel_Bachelier_IV(F0,K,tau,C_real);
toc
% absolute error and its statistics
t2 = abs(IV_Jackel - sigma);
t2_statistics = [mean(t2(1:end)),std(t2(1:end)),max(t2),min(t2)];

% Dekker-Brent method
IV_Brent = zeros(1,length(sigma));
tic
for i = 1:length(sigma)
    options = optimset('Display','off');
    sigmaN_0 = 1;
    IV_Brent(i) = fzero(@(sigmaN)obj_Brent(sigmaN,F0,K(i),tau(i),C_real(i)),1,options);
end
toc
% absolute error and its statistics
t3 = abs(IV_Brent - sigma);
t3_statistics = [mean(t3(1:end)),std(t3(1:end)),max(t3),min(t3)];





