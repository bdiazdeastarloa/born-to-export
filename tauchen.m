function [X,p_x,erg_px] = tauchen(sig_eps_in,size_in,mean_in,length_in,rho_in) 

% -------------------------------------------------------------------------
% Born to export 
%
% tauchen: constructs shock grids and transition matrices for shock X.
% Arguments are: 
% - sig_eps: std. dev. of innovation.
% - size_in: # of discretized cells.
% - mean_in: mean of process.
% - length_in: grid length.
% - rho_in: root of process.
%
% Written by James Tybout (2010).
% Edited and revised by Bernardo Diaz de Astarloa (2011-13).
% -------------------------------------------------------------------------

% Grid
[X,step_size] = tauchen_grid(sig_eps_in,size_in,mean_in,length_in,rho_in);

% Transition densities matrix
p_x      = zeros(size_in,size_in);                                      
p_x(:,1) = normcdf(((X(1)+step_size/2)*ones(size_in,1)-rho_in*X')/sig_eps_in,0,1);
p_x(:,size_in) = ones(size_in,1)-normcdf(((X(size_in)-step_size/2)*ones(size_in,1)-rho_in*X')/sig_eps_in,0,1);
for j=2:(size_in-1)
    p_x(:,j)=normcdf(((X(j)+step_size/2)*ones(size_in,1)-rho_in*X')/sig_eps_in,0,1)...
          -normcdf(((X(j)-step_size/2)*ones(size_in,1)-rho_in*X')/sig_eps_in,0,1);
end;

% Ergodic distribution
erg_px = (1/size_in)*ones(size_in,1);
test = 1;
while test > 10^(-8)
    erg_px1 = p_x'*erg_px;
    test    = max(abs(erg_px1-erg_px));
    erg_px  = erg_px1;
end