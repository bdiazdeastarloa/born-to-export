% -------------------------------------------------------------------------
% Born to export 
%
% plot_trajs: plots trajectories (outcoes of simulation) of shocks and main
% variables (mostly aggregates) of the model.
%
% Written by Bernardo Diaz de Astarloa (2011-13).
% -------------------------------------------------------------------------


%% Prepare variables

global T
close all
% Load shocks, product appeals, etc.
load traj_par

TT = T;
tt = 1:1:TT;

% No sunk cost & no scrap value
load traj_spec0;
exporters_00 = log(totexdum(1:TT));         % log # of exporters
exports_00   = log(totexp(1:TT));           % log total exports
avg_exp_00   = totexp(1:TT)./totexdum(1:TT);% avg. exports per firm
agg_aa_00    = squeeze(agg_aa);
agg_bb_00    = squeeze(agg_bb);
hazz_00      = mean(hazard_f);
active_00    = active;
clear totexdum totexp agg_aa agg_bb hazard_f active

% Sunk cost & scrap value
load traj_spec2;
exporters_11 = log(totexdum(1:TT));         % # of exporters
exports_11   = log(totexp(1:TT));           % total exports
avg_exp_11   = totexp(1:TT)./totexdum(1:TT);% avg. exports per firm
agg_aa_11    = squeeze(agg_aa);
agg_bb_11    = squeeze(agg_bb);
hazz_11      = mean(hazard_f);
active_11    = active;
clear totexdum totexp agg_aa agg_bb hazard_f active

% Sunk cost & no scrap value
load traj_spec1;
exporters_10 = log(totexdum(1:TT));         % # of exporters
exports_10   = log(totexp(1:TT));           % total exports
avg_exp_10   = totexp(1:TT)./totexdum(1:TT);% avg. exports per firm
agg_aa_10    = squeeze(agg_aa);
agg_bb_10    = squeeze(agg_bb);
hazz_10      = mean(hazard_f);
active_10    = active;
clear totexdum totexp agg_aa agg_bb hazard_f active



%% Plot figures

% -------------------------------------------------------------------------
% Product appeal distribution (mu_f, mu_h)
% -------------------------------------------------------------------------
% figure('visible','off')
% scatterhist(mu_h,mu_f);
% xlabel('Home product appeal','Fontsize',20);
% ylabel('Foreign product appeal','Fontsize',20);
% title('Distribution of product appeal','Fontsize',20);
% set(gca,'fontsize',20);
% h = gcf;
% set(h,'PaperOrientation','landscape');
% set(h,'PaperUnits','normalized');
% set(h,'PaperPosition', [0 0 1 1]);
% print(gcf, '-dpdf', 'figures/mu_dist.pdf')

% -------------------------------------------------------------------------
% Sales
% -------------------------------------------------------------------------
% totsal = sale_h + export;
% pos_sal = totsal > 0;
% totsal2 = log(totsal(pos_sal));
% msal = mean(totsal2);
% s_sal = std(totsal2);
% totsal2 = (totsal2-msal)/s_sal;
% bins = -5:0.02:5;
% 
% figure(2)
% hist(totsal2,bins);  
% title('Histogram of normalized log sales')
% 
% totsal2 = sort(totsal2);
% [nr,~] = size(totsal2);
% counter = 1:1:nr;
% log_cum = log(counter'/nr);
% 
% figure(3)
% scatterhist(totsal2,log_cum);
% title('Normalized log sales cumulative distribution')


% -------------------------------------------------------------------------
% Total exports and foreign shock trajectories
% -------------------------------------------------------------------------
figure('visible','off');
[AX,x00,B] = plotyy(tt,exports_00,tt,xr_traj,'plot');
hold on;
[x10] = plot(tt,exports_10);
[x11] = plot(tt,exports_11);
title('Total exports and foreign shock trajectories','Fontsize',20)
set(get(AX(1),'Ylabel'),'String','Total exports (log)','Fontsize',20) 
set(get(AX(2),'Ylabel'),'String','Foreign shock (log)','Fontsize',20)
set(AX(2),'xticklab',[],'xtick',[])
set(AX(2),'Fontsize',20)
set(x00,'LineStyle','-','Color','b','Linewidth',1.5)
set(x10,'LineStyle','--','Color','r','Linewidth',1.5)
set(x11,'LineStyle','-.','Color','b','LineWidth',2)
set(B,'LineStyle','-','Marker','o','Linewidth',1.5)
xlabel('Period','Fontsize',20)
legend('K=0,psi=0','K>0,psi=0','K>0,psi>0','Location','NorthWest');
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/traj_exports.pdf')
hold off;

% -------------------------------------------------------------------------
% # of exporters and foreign shock trajectories
% -------------------------------------------------------------------------
figure('visible','off')
[AX,x00,B] = plotyy(tt,exporters_00,tt,xr_traj,'plot');
hold on;
[x10] = plot(tt,exporters_10);
[x11] = plot(tt,exporters_11);
title('Number of exporters and foreign shock trajectories','Fontsize',20)
set(get(AX(1),'Ylabel'),'String','Exporters (log)','Fontsize',20) 
set(get(AX(2),'Ylabel'),'String','Foreign shock (log)','Fontsize',20) 
set(AX(1),'YLim',[0 10])
set(AX(1),'YTick',0:2:10)
set(AX(2),'xticklab',[],'xtick',[])
set(AX(2),'Fontsize',20)
set(x00,'LineStyle','-','Color','b','Linewidth',1.5)
set(x10,'LineStyle','--','Color','r','Linewidth',1.5)
set(x11,'LineStyle','-.','Color','b','LineWidth',2)
set(B,'LineStyle','-','Marker','o','Linewidth',1.5)
xlabel('Period','Fontsize',20)
legend('K=0,psi=0','K>0,psi=0','K>0,psi>0','Location','NorthWest')
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/traj_exporters.pdf')
hold off


% -------------------------------------------------------------------------
% # of active firms
% -------------------------------------------------------------------------
figure%('visible','off')
[AX,x00,B] = plotyy(tt,active_10,tt,xr_traj,'plot');
hold on
[x11] = plot(tt,active_11);
set(get(AX(1),'Ylabel'),'String','Number of firms','Fontsize',20) 
set(get(AX(2),'Ylabel'),'String','Foreign shock','Fontsize',20)
set(AX(2),'xticklab',[],'xtick',[])
set(x00,'LineStyle','-','Color','b','Linewidth',1.5)
set(x11,'LineStyle','--','Color','r','Linewidth',1.5)
%set(x10,'LineStyle','-.','Color','b','LineWidth',1.5)
set(B,'LineStyle','-','Marker','o')
title('Number of active firms and foreign shock trajectories','Fontsize',20)
legend('K>0,phi=0','K>0,phi>0','Location','NorthWest')
xlabel('Period','Fontsize',20)
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/traj_activefirms.pdf')
hold off


% -------------------------------------------------------------------------
% Posterior theta distribution
% -------------------------------------------------------------------------
% figure('visible','off')
% plot(tt(1:TT),mean(mu_f_hat(:,1:TT)));
% title('Posterior theta distribution','Fontsize',20);
% xlabel('Period','Fontsize',20);
% set(gca,'fontsize',20);
% h = gcf;
% set(h,'PaperOrientation','landscape');
% set(h,'PaperUnits','normalized');
% set(h,'PaperPosition', [0 0 1 1]);
% print(gcf, '-dpdf', 'figures/post_theta.pdf')


% -------------------------------------------------------------------------
% # of sucessful and failed matches
% -------------------------------------------------------------------------
figure('visible','off')
plot(tt(1:TT),agg_aa_00,'Linewidth',1.5);
hold on;
plot(tt(1:TT),agg_bb_00,'LineStyle','--','Linewidth',1.5);
title('Evolution of succesful and failed matches (K=0,psi=0)','Fontsize',20);
legend('Successes','Failures','Location','NorthWest');
xlabel('Period','Fontsize',20);
ylabel('Successes and failures','Fontsize',20);
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/traj_ab-evol_00.pdf')
hold off;

figure('visible','off')
title('Evolution of succesful and failed matches (K>0,psi=0)','Fontsize',20);
plot(tt(1:TT),agg_aa_10,'Linewidth',1.5);
hold on;
plot(tt(1:TT),agg_bb_10,'LineStyle','--','Linewidth',1.5);
legend('Successes','Failures','Location','NorthWest');
xlabel('Period','Fontsize',20);
ylabel('Successes and failures','Fontsize',20);
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/traj_ab-evol_10.pdf')
hold off;

figure('visible','off')
plot(tt(1:TT),agg_aa_11,'Linewidth',1.5);
hold on;
plot(tt(1:TT),agg_bb_11,'LineStyle','--','Linewidth',1.5);
title('Evolution of succesful and failed matches (K>0,psi>0)','Fontsize',20);
legend('Successes','Failures','Location','NorthWest');
xlabel('Period','Fontsize',20);
ylabel('Successes and failures','Fontsize',20);
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/traj_ab-evol_11.pdf')
hold off;
