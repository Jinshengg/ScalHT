clear; close all; clc;
% profile on 
% profile clear
%% parameters 

n = 1023; % length of sinusoid signal 
% n1 = (n+1)/2;
% n2 = n+1-n1;
s = 512; %  dim of channels s =128 r = 6
r = 6; % number of point sources
maxit = 80; % maximum iteratioins
p = 0.12; % sampling ratio 
m = round(n*s*p);% number of observations
sep = 0 ; % frequency separation 
tol = 1e-7;
mode = 0; % run mode. 0, run time test mode; 1 phase transition test mode
eta = 0.4;
%% generate multi-channel signals
[fs, cs, H, X_star, Xs,Omega] = getSignals_mulc(r, s, n, m, sep);
% 
% 
[~,ScalHT_err,ScalHT_time] = ScalHT(Xs,X_star,Omega,p,n,r,s,maxit,tol,mode,eta);


% fiht_time(end)

% HTC_time(end)
% ScalGD_time(end)

% profile viewer

clrs = {[.5,0,.5], [1,.5,0], [1,0,0], [0,.5,0], [0,0,1]};
mks = {'o', 'x', 'p', 's', 'd'};
figure('Position', [0,0,800,600], 'DefaultAxesFontSize', 20);
lgd = {};
    %% plot ScalHT versus iteration
    errors = ScalHT_err;
    t_subs = 1:10:length(errors);
    semilogy(t_subs, errors(t_subs), 'Color', 'r', 'Marker', mks{1}, 'MarkerSize', 10,'LineWidth',1.5);
    hold on; grid on;
    lgd{end+1} = sprintf('$\\mathrm{ScalHT}$');
    
     
    xlabel('Iteration count');
    ylabel('Relative error');
    legend(lgd, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 24);
%     fig_name = sprintf('test_performance_n=%d_r=%d_s=%d_p=%g', n, r,s, p);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 16 12 ]);
    set(gca,'FontName','times new roman','FontSize',24,'Layer','top','linewidth',1.5,'Gridalpha',0.1,'MinorGridAlpha',0.1);
    
    figure('Position', [0,0,800,600], 'DefaultAxesFontSize', 20);
    %% plot ScalHT versus time
    errors = ScalHT_err;
    t_subs = 1:10:length(errors);
    semilogy(ScalHT_time(t_subs), errors(t_subs), 'Color', 'r', 'Marker', mks{1}, 'MarkerSize', 10,'LineWidth',1.5);
    hold on; grid on;
%     lgd{end+1} = sprintf('$\\mathrm{ScalHT}$');
    
     
    xlabel('Run time (secs)');
    ylabel('Relative error');
    legend(lgd, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 24);
%     fig_name = sprintf('test_performance_n=%d_r=%d_s=%d_p=%g', n, r,s, p);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 16 12 ]);
    set(gca,'FontName','times new roman','FontSize',24,'Layer','top','linewidth',1.5,'Gridalpha',0.1,'MinorGridAlpha',0.1);