function []=plot_policies(spec,mm)

% -------------------------------------------------------------------------
% Born to export 
%
% plot_policies: plot search policy and value function.
%
% Written by Bernardo Diaz de Astarloa (2011).
% -------------------------------------------------------------------------


%% Prepare search policy and value functions

% Load shocks and solution
load erg_dist;

if spec==0
    load sol_bte_spec0;
elseif spec==1
    load sol_bte_spec1;
else
    load sol_bte_spec2;
end
    
n_size = mm.n_size;

% Set matrices for plots
lambda_avg = zeros(n_size,n_size); 
v_avg      = zeros(n_size,n_size);      
lambda_hi  = zeros(n_size,n_size); 
v_hi       = zeros(n_size,n_size); 
lambda_lo  = zeros(n_size,n_size); 
v_lo       = zeros(n_size,n_size); 
% Set desired prod and type (in general type dim is 1!)
prod_hi = 20;
prod_lo = 5;
typ  = 1;

for i = 1:n_size
    for j = 1:n_size
        % Average over prod and foreign shocks
        v_avg(i,j) = erg_phi'*v_new{i,j,typ}*erg_xf;                %#ok
        lambda_avg(i,j) = erg_phi'*lambda_f{i,j,typ}*erg_xf;        %#ok
        % Average over foreign shock for high productivity
        v_hi(i,j) = v_new{i,j,typ}(prod_hi,:)*erg_xf;
        lambda_hi(i,j) = lambda_f{i,j,typ}(prod_hi,:)*erg_xf;
        % Average over foreign shock for low productivity
        v_lo(i,j) = v_new{i,j,typ}(prod_lo,:)*erg_xf;
        lambda_lo(i,j) = lambda_f{i,j,typ}(prod_lo,:)*erg_xf;
    end
end


%% Search and value functions int. over shocks

figure('visible','off')
mesh(lambda_avg(:,:));
title('Search intensity (integrated over shocks)');
xlabel('Unsuccessful matches');
ylabel('Successful matches');
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/search_avg.pdf')

figure('visible','off')
mesh(v_avg(:,:));
title('Value function (integrated over shocks)');
xlabel('Unsuccessful matches');
ylabel('Successful matches');
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/valuef_avg.pdf')

figure('visible','off')
mesh(lambda_hi(:,:));
title('Search intensity (high productivity)');
xlabel('Unsuccessful matches');
ylabel('Successful matches');
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/search_hi.pdf')

figure('visible','off')
mesh(v_hi(:,:));
title('Value function (high productivity)');
xlabel('Unsuccessful matches');
ylabel('Successful matches');
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/valuef_hi.pdf')

figure('visible','off')
mesh(lambda_lo(:,:));
title('Search intensity (low productivity)');
xlabel('Unsuccessful matches');
ylabel('Successful matches');
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/search_lo.pdf')

figure('visible','off')
mesh(v_lo(:,:));
title('Value function (low productivity)');
xlabel('Unsuccessful matches');
ylabel('Successful matches');
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/valuef_lo.pdf')


%% Search policy for fixed successes and failures (varying shocks)

figure('visible','off')
mesh(lambda_f{20,20,typ});
title('Search intensity: (a,b)=(20,20)');
xlabel('Foreign shock');
ylabel('Productivity');
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/search_shocks_a20b20.pdf')

figure('visible','off')
mesh(lambda_f{5,20,typ});
title('Search intensity, (a,b)=(5,20)');
xlabel('Foreign shock');
ylabel('Productivity');
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/search_shocks_a5b20.pdf')

figure('visible','off')
mesh(lambda_f{1,1,typ});
title('Search intensity, (a,b)=(1,1)');
xlabel('Foreign shock');
ylabel('Productivity');
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/search_shocks_a1b1.pdf')

figure('visible','off')
mesh(v_new{20,1,typ});
title('Value function, (a,b)=(5,20)');
xlabel('Foreign shock');
ylabel('Productivity');
set(gca,'fontsize',20);
h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', 'figures/valuef_shocks_a5b20.pdf')

