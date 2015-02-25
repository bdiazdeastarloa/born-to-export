function [lambda_f,lambda_h,pi_tilda_h,pi_tilda_f,p_phi,p_xf,p_xh,v_new,e_hazz,chi_c,chi_e,exitflag] = solve_bte(mm)

% -------------------------------------------------------------------------
% Born to export 
%
% solve_bte: solves the model, computing policy and value functions.
%
% Written by James Tybout (2010).
% Edited and revised by Bernardo Diaz de Astarloa (2011-13).
% -------------------------------------------------------------------------


%% Assign parameter values

scale_f     = mm.scale_f;   
scale_h     = mm.scale_h;
mean_phi    = mm.mean_phi;     
rho_phi     = mm.rho_phi;      
sig_eps_phi = mm.sig_eps_phi; 
mean_z      = mm.mean_z;     
rho_z       = mm.rho_z;      
sig_eps_z   = mm.sig_eps_z; 
mean_xh     = mm.mean_xh;   
rho_xh      = mm.rho_xh;     
sig_eps_xh  = mm.sig_eps_xh; 

mean_xf     = mm.mean_xf;   
rho_xf      = mm.rho_xf;     
sig_eps_xf  = mm.sig_eps_xf; 

grid_length = mm.grid_length;  
%n_size      = mm.n_size;        
phi_size    = mm.phi_size;     
x_size      = mm.x_size;        
z_size      = mm.z_size;
%dim0        = mm.dim0;

%ahf         = mm.ahf;
%bhf         = mm.bhf;
%theta0      = mm.theta0;


%% Start solution routine

% Grids and transition matrices for exogenous variables:    
% match-specific shocks 
[Z,p_z,erg_pz] = tauchen(sig_eps_z, z_size, mean_z, grid_length, rho_z);         
% productivity
[Phi,p_phi,erg_phi] = tauchen(sig_eps_phi, phi_size, mean_phi, grid_length, rho_phi);
% home macro shocks
[X_h,p_xh,erg_xh] = tauchen(sig_eps_xh, x_size, mean_xh, grid_length, rho_xh); 
% foreign macro shocks
[X_f,p_xf,erg_xf] = tauchen(sig_eps_xf, x_size, mean_xf, grid_length, rho_xf); 
   
save('erg_dist','erg_pz','erg_phi','erg_xh','erg_xf');        
    
% Solve for expected payoffs to matches at home and abroad
pi_tilda_h = pi_tilda(Z, p_z, X_h, p_xh, Phi, p_phi, scale_h, erg_pz, mm);
pi_tilda_f = pi_tilda(Z, p_z, X_f, p_xf, Phi, p_phi, scale_f, erg_pz, mm);
         
% Solve for firms' search effort at home at and abraod
[lambda_h,~] = policy_h(pi_tilda_h,mm);
[lambda_f,e_hazz,v_new,chi_c,chi_e,exitflag] = policy_f_bte(pi_tilda_f,p_phi,p_xf,mm);
