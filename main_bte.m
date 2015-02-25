% -------------------------------------------------------------------------
% Born to export 
%
% main_bte: setup parameters, call solver, simulate trajectories, compute
% moments.
%
% This routine can be used in estimation or standalone to solve the model,
% simulate and generate figures for a specific parametrization.
%
% Written by James Tybout (2010).
% Edited and revised by Bernardo Diaz de Astarloa (2011-14).
% -------------------------------------------------------------------------

%% Settings 

global spec linear

% Parameters from EEJKT.
par = [0.1891  34.1008  38.5289  2.8870  2.5967  0.6033  0.8913  0.7724];

fprintf('Initializing... ')
% Set linear (1)/Cobb-Douglas(0) specification for signals
linear = 1;

ah        = ahf;                  % Success parameter, beta distribution (home)
bh        = bhf;                  % Failure parameter, beta distribution (home)
af        = ah;                   % Success parameter, beta distribution (foreign)
bf        = bh;                   % Failure parameter, beta distribution (foreign)

f         = par(1);               % search cost
F         = par(2);               % fixed cost of maint. a relationship
beta      = par(3);               % cost function parameter
scale_f   = par(4);               % export profits scale parameter
scale_h   = par(5);               % home profits scale parameter
rho_z     = par(6);               % root of product appeal shock
sig_eps_z = par(7);               % variance of seller-specific effect 
rho_phi   = par(8);               % root of productivity shock

specs      = zeros(3,2);          % ecost/psi specifications holder          
specs(1,:) = [3000 0];
specs(2,:) = [3000 0.3];
specs(3,:) = [0 0];

flags = zeros(3,1);

% Put parameters in structure and set additional parameters 
% (Creates structure 'mm' used as argument below)
setparams;

fprintf('Done.\n')
fprintf('\n')
    
% Set sunk cost/scrap value specification.
spec  = 2;
% for spec=0:2
    mm.ecost = specs(spec,1);               % Sunk cost of creating an establishment
    mm.psi   = specs(spec,2);               % Fraction recovered upon exit
    mm.scrap = mm.psi*mm.ecost;             % Scrap value
    
    if spec==1
        savefile = 'sol_bte_spec1';
    %    fprintf('Specification: sunk costs & no scrap value.\n')
    elseif spec==2
        savefile = 'sol_bte_spec2';
    %    fprintf('Specification: sunk costs & scrap value.\n')
    else
        savefile = 'sol_bte_spec0';
    %    fprintf('Specification: no sunk costs & no scrap value.\n')
    end

    
    %% Solve and simulate
    
    %fprintf('Solving the model... ')
    tic
    [lambda_f,lambda_h,pi_tilda_h,pi_tilda_f,p_phi,p_xf,p_xh,v_new,e_hazz,chi_c,chi_e,exitflag] = solve_bte(mm);
    time = toc/60;
    if exitflag == 1
        %fprintf('Done.\n')
        fprintf('\nBTE model solved: %2i.\n', time);
        fprintf('\n');
        
        %fprintf('Saving BTE solution... ');
        save(savefile,'lambda_h','lambda_f','e_hazz','v_new','chi_c','chi_e','pi_tilda_h','pi_tilda_f');
        %fprintf('Done.\n')
    else
        fprintf('There was a problem solving the model. Exiting now...\n');
    end
    flags(spec) = exitflag;
    
    % Simulate firm trajectories
    if exitflag==1
        %fprintf('\n');
        fprintf('Simulating trajectories... ');
        MM = trajec_bte(mm,pi_tilda_h,pi_tilda_f,lambda_f,lambda_h,e_hazz,chi_c,chi_e);
        fprintf('--------------------------------------------------------------\n')
    end
% end


%% Generate figures

% Comment this out for estimation.
%{
if prod(flags)==1
    fprintf('\n');
    fprintf('Generating figures... ');
    plot_compare(mm);
    plot_trajs;
    fprintf('Done.\n')
else
    fprintf('\n');
    fprintf('There was a problem solving some specification.\n');
    fprintf('No figures generated.\n');
end

clear all
close all

fprintf('\n');
fprintf('\n');
fprintf('The End.\n');
%}    


                       