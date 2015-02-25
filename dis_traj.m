function [x_disc,indx] = dis_traj(SS,TT,X,x_size,x_shock,x_mean,rho_x) 

% -------------------------------------------------------------------------
% Born to export 
%
% dis_traj: generates trajectories of shocks and maps them to a grid.
%
% Written by James Tybout (2010).
% Edited and revised by Bernardo Diaz de Astarloa (2011-2013).
% -------------------------------------------------------------------------

indx   = zeros(SS,TT);
x_disc = zeros(SS,TT);
x_traj = zeros(SS,TT);

% Generate trajectories and map them onto discrete grid

x_traj(:,1) = x_mean + x_shock(:,1)/sqrt(1-rho_x^2); 
dif = abs(x_traj(:,1)*ones(1,x_size)-ones(SS,1)*X); 
for ii=1:SS
 [~,indx(ii,1)] = min(dif(ii,:));
end;
x_disc(:,1)=X(indx(:,1));   % map to grid

for t=2:TT
    x_traj(:,t)= x_mean*(1-rho_x) + rho_x*x_traj(:,t-1)+x_shock(:,t); 
    dif = abs(x_traj(:,t)*ones(1,x_size)-ones(SS,1)*X);
    for ii=1:SS
        [~,indx(ii,t)] = min(dif(ii,:));
    end;
    x_disc(:,t)=X(indx(:,t));  
end

