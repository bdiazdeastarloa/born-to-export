function [pi_tilda]= pi_tilda(ZZ,p_z,XX,p_x,Phi,p_phi,scal,erg_pz,mm)

% -------------------------------------------------------------------------
% Born to export 
%
% pi_tilda: computes payoffs to a match.
%
% Written by James Tybout (2010).
% Edited and revised by Bernardo Diaz de Astarloa (2011-13).
% -------------------------------------------------------------------------

% Assign parameter values
beta_adj     = mm.beta_adj;
F            = mm.F;
eta          = mm.eta;
pi_tol       = mm.pi_tol;
phi_size     = mm.phi_size;        % # of productivity states     
x_size       = mm.x_size;          % # of macro states  
z_size       = mm.z_size;          % # of buyer states


% Expected match payoff (conditioned on buyer state z)
payoffs = kron(exp((eta-1).*Phi'),exp(XX));
        
pi_tilda_old = cell(z_size);
pi_tilda_new = cell(z_size);

for j = 1:z_size
    pi_tilda_old{j} = zeros(phi_size,x_size);
    pi_tilda_new{j} = zeros(phi_size,x_size);  
end

% Iterate on contraction
eps=1;
while eps > pi_tol   
    for j = 1:z_size
        c_val = zeros(phi_size,x_size);
        for k = 1:z_size
            c_val = c_val + p_z(j,k)*p_phi*pi_tilda_old{k}*p_x';
        end
        pi_tilda_new{j} = scal*exp(ZZ(j))*payoffs + max(0,beta_adj*c_val-F); 
    end
    eps = norm(cell2mat(pi_tilda_old)-cell2mat(pi_tilda_new))/norm(cell2mat(pi_tilda_new)) ;
    pi_tilda_old = pi_tilda_new;
end
         
% Expected profits (over buyer states)
pi_tilda = zeros(phi_size,x_size);
for j=1:z_size             
    pi_tilda = pi_tilda + pi_tilda_new{j}*erg_pz(j) ;
end

     