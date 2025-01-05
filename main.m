clear; close all; clc;
% profile on 
% profile clear
%% parameters 

n = 1023; % length of sinusoid signal 
% n1 = (n+1)/2;
% n2 = n+1-n1;
s = 128; %  dim of channels s =128 r = 6
r = 10; % number of point sources
maxit = 100; % maximum iteratioins
p = 0.125; % sampling ratio 0.025
m = round(n*s*p);% number of observations
sep = 0 ; % frequency separation 
tol = 1e-7;
mode = 0; % run mode. 0, run time test mode; 1 phase transition test mode
%% generate multi-channel signals
[fs, cs, H, X_star, Xs,Omega] = getSignals_mulc(r, s, n, m, sep);
% 
[fiht_err,fiht_time] = FIHT_fh(Xs,X_star,Omega,p,n,r,s,maxit,tol,mode);
% 
[HTC_err,HTC_time] = HTC_Fast(Xs,X_star,Omega,p,n,r,s,maxit,tol,mode);

% [ScalGD_err,ScalGD_time] = ScaledGD(Xs,X_star,Omega,p,n,r,s,maxit,tol,mode);

fiht_time(end)

HTC_time(end)
% ScalGD_time(end)

% profile viewer

clrs = {[.5,0,.5], [1,.5,0], [1,0,0], [0,.5,0], [0,0,1]};
mks = {'o', 'x', 'p', 's', 'd'};
figure('Position', [0,0,800,600], 'DefaultAxesFontSize', 20);
lgd = {};
    %% plot HTC
    errors = HTC_err;
    t_subs = 1:10:length(errors);
    semilogy(t_subs, errors(t_subs), 'Color', clrs{1}, 'Marker', mks{1}, 'MarkerSize', 9);
    hold on; grid on;
    lgd{end+1} = sprintf('$\\mathrm{ScalHT}$');
    
    %% plot FIHT
    errors = fiht_err;
    t_subs = 1:10:length(errors);
    semilogy(t_subs, errors(t_subs), 'Color', clrs{2}, 'Marker', mks{2}, 'MarkerSize', 9);
    hold on; grid on;
    lgd{end+1} = sprintf('$\\mathrm{FIHT}$');
    
%     %% plot ScalGD
%     errors = ScalGD_err;
%     t_subs = 1:10:length(errors);
%     semilogy(t_subs, errors(t_subs), 'Color', clrs{3}, 'Marker', mks{3}, 'MarkerSize', 9);
%     hold on; grid on;
%     lgd{end+1} = sprintf('$\\mathrm{ScaledGD}$');
   
   
    xlabel('Iteration count');
    ylabel('Relative error');
    legend(lgd, 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 24);
    fig_name = sprintf('test_performance_n=%d_r=%d_s=%d_p=%g', n, r,s, p);