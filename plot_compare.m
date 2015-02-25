function plot_compare(mm)

% -------------------------------------------------------------------------
% Born to export 
%
% plot_compare: plots policy and value functions comparing different sunk
% costs/scrap values specifications.
%
% Written by Bernardo Diaz de Astarloa (2011).
% -------------------------------------------------------------------------

close all

% Load shocks
load erg_dist
n_size = mm.n_size;


%% 1. Search policy and value function for K=0 and K>0

% Load K=0 objects
load sol_bte_spec0
lambda_f_00 = lambda_f;
v_new_00    = v_new;
chi_c_00    = chi_c;
chi_e_00    = chi_e;
clear v_new lambda_f chi_c chi_e;

% Load K>0, psi=0 objects
load sol_bte_spec1
%lambda_f_10 = lambda_f;
%v_new_10    = v_new;
%chi_c_10    = chi_c;
chi_e_10    = chi_e;
clear v_new lambda_f chi_c chi_e;

% Load K>0, psi>0 objects
load sol_bte_spec2
lambda_f_11 = lambda_f;
v_new_11    = v_new;
chi_c_11    = chi_c;
chi_e_11    = chi_e;
clear v_new lambda_f chi_c chi_e;

% High productivity search and value
lambda_avg_00  = zeros(n_size,n_size);
v_avg_00       = zeros(n_size,n_size);
lambda_avg_11  = zeros(n_size,n_size);
v_avg_11       = zeros(n_size,n_size);
lambda_high_00 = zeros(n_size,n_size);
v_high_00      = zeros(n_size,n_size);
lambda_high_11 = zeros(n_size,n_size);
v_high_11      = zeros(n_size,n_size);
%lambda_low_00 = zeros(n_size,n_size);
v_low_00      = zeros(n_size,n_size);
%lambda_low_11 = zeros(n_size,n_size);
v_low_11      = zeros(n_size,n_size);

% Set desired productivity level and desired global type
% (in general global type is dim 1!)
prod_hi = 20; 
prod_lo = 5;
typ  = 1;

for i = 1:n_size
    for j = 1:n_size
        % Average over productivity and foreign shocks
        v_avg_00(i,j)      = erg_phi'*v_new_00{i,j,typ}*erg_xf;
        lambda_avg_00(i,j) = erg_phi'*lambda_f_00{i,j,typ}*erg_xf;
        v_avg_11(i,j)      = erg_phi'*v_new_11{i,j,typ}*erg_xf;
        lambda_avg_11(i,j) = erg_phi'*lambda_f_11{i,j,typ}*erg_xf;
        % Average over foreign shocks for high prod
        v_high_00(i,j)      = v_new_00{i,j,typ}(prod_hi,:)*erg_xf;
        lambda_high_00(i,j) = lambda_f_00{i,j,typ}(prod_hi,:)*erg_xf;
        v_high_11(i,j)      = v_new_11{i,j,typ}(prod_hi,:)*erg_xf;
        lambda_high_11(i,j) = lambda_f_11{i,j,typ}(prod_hi,:)*erg_xf;
        % Average over foreign shocks for low prod
        v_low_00(i,j)      = v_new_00{i,j,typ}(prod_lo,:)*erg_xf;
        %lambda_low_00(i,j) = lambda_f_00{i,j,typ}(prod_lo,:)*erg_xf;
        v_low_11(i,j)      = v_new_11{i,j,typ}(prod_lo,:)*erg_xf;
        %lambda_low_11(i,j) = lambda_f_11{i,j,typ}(prod_lo,:)*erg_xf;
    end
end
  
% Differences in search and value 
diff_l_avg  = lambda_avg_11 - lambda_avg_00;
diff_v_avg  = v_avg_11 - v_avg_00;
diff_l_high = lambda_high_11 - lambda_high_00;
diff_v_high = v_high_11 - v_high_00;
% diff_l_low  = lambda_low_11 - lambda_low_00;
diff_v_low  = v_low_11 - v_low_00;

% Diffs of diffs to illustrate small differences in value functions
dif_dif_av_hi = diff_v_avg - diff_v_high; 
dif_dif_lo_hi = diff_v_low - diff_v_high;

% Generate figures
% policies_diff(diff_l_avg,diff_v_avg,diff_l_high,diff_v_high);
policies_dif_dif(dif_dif_av_hi, dif_dif_lo_hi);

% Create figure
figure1 = figure('visible','off');
plot1 = subplot(1,1,1,'Parent',figure1,'YTick',[0 10 20 30 40 50],...
    'XTick',[0 10 20 30 40 50],...
    'FontSize',20);
zlim(plot1,[-0.05 0.07]);
view(plot1,[-37.5 30]);
grid(plot1,'on');
hold(plot1,'all');
mesh(diff_l_avg,'Parent',plot1);
title('Search intensity (average)','FontSize',20);
xlabel('Failures','FontSize',20);
ylabel('Successes','FontSize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/diff_l_avg.pdf');

% Create figure
figure1 = figure('visible','off');
plot1 = subplot(1,1,1,'Parent',figure1,'YTick',[0 10 20 30 40 50],...
    'XTick',[0 10 20 30 40 50],...
    'FontSize',20);
zlim(plot1,[-0.05 0.07]);
view(plot1,[-37.5 30]);
grid(plot1,'on');
hold(plot1,'all');
mesh(diff_l_high,'Parent',plot1);
title('Search intensity (high productivity)','FontSize',20);
xlabel('Failures','FontSize',20);
ylabel('Successes','FontSize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/diff_l_high.pdf');


%% 2. Entry, exit comparison and hysteresis band

% Entering and continuing indicators
% Fixed (a,b)=(1,10), for all prod and foreign shocks
a = 1;
b = 10;
chi_c_diff = chi_c_11{a,b,1} - chi_c_00{a,b,1};

% Entry policy: difference between specifications
chi_e_diff_1 = chi_e_11{1} - chi_e_00{1};
chi_e_diff_2 = 2*chi_e_11{1} - chi_e_10{1};

% Hysteresis band
hyst = chi_c_11{1,1,1} - chi_e_11{1};

% Continuing policy (average over foreign shock and prod)
% Difference between specifications
chi_c_diff_avg = zeros(n_size,n_size);
for i = 1:n_size
    for j = 1:n_size 
        chi_c_diff_avg(i,j) = erg_phi'*chi_c_11{i,j,typ}*erg_xf - erg_phi'*chi_c_00{i,j,typ}*erg_xf;
    end
end

figure('visible','off')
mesh(chi_c_diff);
title('Continuing policy (a,b)=(1,10)','FontSize',20);
xlabel('Foreign shock','FontSize',20);
ylabel('Productivity','FontSize',20);
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/cont_diff.pdf');

figure('visible','off')
mesh(chi_e_diff_1);
title('Entry policy','FontSize',20);
xlabel('Foreign shock','FontSize',20);
ylabel('Productivity','FontSize',20);
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/entry_diff_1.pdf');

figure('visible','off')
mesh(chi_e_diff_2);
title('Entry policy','FontSize',20);
xlabel('Foreign shock','FontSize',20);
ylabel('Productivity','FontSize',20);
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/entry_diff_2.pdf');

figure('visible','off')
mesh(hyst);
title('Hysteresis band','FontSize',20);
xlabel('Foreign shock','FontSize',20);
ylabel('Productivity','FontSize',20);
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/hysteresis_band.pdf');

end


%% Helper functions to graph subplots flexibly

function policies_diff(zdata1, zdata2, zdata3, zdata4)
%CREATEFIGURE(ZDATA1, ZDATA2, ZDATA3, ZDATA4)
%  ZDATA1:  surface zdata
%  ZDATA2:  surface zdata
%  ZDATA3:  surface zdata
%  ZDATA4:  surface zdata

%  Auto-generated by MATLAB on 04-Feb-2015 21:44:39

% Create figure
figure1 = figure('visible','off');

% Create subplot
subplot1 = subplot(2,2,1,'Parent',figure1,'YTick',[0 10 20 30 40 50],...
    'XTick',[0 10 20 30 40 50],...
    'FontSize',16);

zlim(subplot1,[-0.05 0.05]);
view(subplot1,[-37.5 30]);
grid(subplot1,'on');
hold(subplot1,'all');

mesh(zdata1,'Parent',subplot1);
title('Search intensity (average)','FontSize',16);
xlabel('Failures','FontSize',16);
ylabel('Successes','FontSize',16);

% Create subplot
subplot2 = subplot(2,2,2,'Parent',figure1,'YTick',[0 10 20 30 40 50],...
    'XTick',[0 10 20 30 40 50],...
    'FontSize',16);
%zlim(subplot2,[-20 1000]);
view(subplot2,[-37.5 30]);
grid(subplot2,'on');
hold(subplot2,'all');

mesh(zdata2,'Parent',subplot2);
title('Value of search (average)','FontSize',16);
xlabel('Failures','FontSize',16);
ylabel('Successes','FontSize',16);

% Create subplot
subplot3 = subplot(2,2,3,'Parent',figure1,'YTick',[0 10 20 30 40 50],...
    'XTick',[0 10 20 30 40 50],...
    'FontSize',16);
zlim(subplot3,[-0.05 0.05]);
view(subplot3,[-37.5 30]);
grid(subplot3,'on');
hold(subplot3,'all');

mesh(zdata3,'Parent',subplot3);
title('Search intensity (high productivity)','FontSize',16);
xlabel('Failures','FontSize',16);
ylabel('Successes','FontSize',16);

% Create subplot
subplot4 = subplot(2,2,4,'Parent',figure1,'YTick',[0 10 20 30 40 50],...
    'XTick',[0 10 20 30 40 50],...
    'FontSize',16);
%zlim(subplot4,[-20 1000]);
view(subplot4,[-37.5 30]);
grid(subplot4,'on');
hold(subplot4,'all');

mesh(zdata4,'Parent',subplot4);
title('Value of search (high productivity)','FontSize',16);
xlabel('Failures','FontSize',16);
ylabel('Successes','FontSize',16);

fig = gcf;

set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/policies_diff.pdf');
end

function policies_dif_dif(zdata1, zdata2)
%CREATEFIGURE(ZDATA1, ZDATA2, ZDATA3, ZDATA4)
%  ZDATA1:  surface zdata
%  ZDATA2:  surface zdata

% Create figure
figure1 = figure('visible','off');

% Create subplot
subplot1 = subplot(1,2,1,'Parent',figure1,'YTick',[0 10 20 30 40 50],...
    'XTick',[0 10 20 30 40 50],...
    'FontSize',15);

zlim(subplot1,[-20 140]);
view(subplot1,[-37.5 30]);
grid(subplot1,'on');
hold(subplot1,'all');

mesh(zdata1,'Parent',subplot1);
title('Diff. of diff.: deltav(avg)-deltav(high)','FontSize',16);
xlabel('Failures','FontSize',16);
ylabel('Successes','FontSize',16);

% Create subplot
subplot2 = subplot(1,2,2,'Parent',figure1,'YTick',[0 10 20 30 40 50],...
    'XTick',[0 10 20 30 40 50],...
    'FontSize',15);
zlim(subplot2,[-20 140]);
view(subplot2,[-37.5 30]);
grid(subplot2,'on');
hold(subplot2,'all');

mesh(zdata2,'Parent',subplot2);
title('Diff. of diff.: deltav(low)-deltav(high)','FontSize',16);
xlabel('Failures','FontSize',16);
ylabel('Successes','FontSize',16);

fig = gcf;

set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/policies_dif_dif.pdf');
end
