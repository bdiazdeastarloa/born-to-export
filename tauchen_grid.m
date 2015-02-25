function [X,step_size] = tauchen_grid(sig_eps_in,size_in,mean_in,length_in,rho_in)  

% -------------------------------------------------------------------------
% Born to export 
%
% tauchen_grid: constructs shock grid and step size for shock X.
% Arguments are: 
% - sig_eps: std. dev. of innovation.
% - size_in: # of discretized cells.
% - mean_in: mean of process.
% - length_in: grid length.
% - rho_in: root of process.
%
% Written by James Tybout (2010).
% -------------------------------------------------------------------------

% Grid and transition matrix
sig_x     = sig_eps_in/sqrt(1-rho_in^2);
step_size = 2*length_in*sig_x/(size_in-1);                                 % width of a discretized interval
X         = (mean_in-length_in*sig_x):step_size:(mean_in+length_in*sig_x); % vector of possible shocks (logs)
