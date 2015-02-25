function [lambda,v_new]= policy_h(pi_tilda,mm)

% -------------------------------------------------------------------------
% Born to export 
%
% policy_h: solves for the home search intensity policy function, using
% FOCs.
%
% Written by James Tybout (2010).
% Edited and revised by Bernardo Diaz de Astarloa (2011-13).
% -------------------------------------------------------------------------


b        = mm.b;   
cost_fn  = mm.cost_fn;
x_size   = mm.x_size;
phi_size = mm.phi_size;
dim0     = mm.dim0;
dim1     = mm.dim1;
hazz     = mm.hazz;             % use hazz: no updating at home
v_new    = cell(dim0,dim1);
lambda   = cell(dim0,dim1);       

for typ = 1:dim0
    for i = 1:dim1
        lambda{i,typ} = zeros(phi_size,x_size);
    end
end

for i = 1:dim1                       % i indexes theta1 (home) values        
    for typ = 1:dim0                 % typ indexes theta0 values  
        
        dvdl = hazz(typ,i).*pi_tilda ;     
        dvdl = max(dvdl,10^(-5));
      
        % Solve FOC to get search policy
        lam = max(1 - (b./dvdl).^(0.5),0);   
        % Compute value function
        v_1 = -cost_fn(lam) + lam.*dvdl;   
        test = v_1 >= 0 ;               
    
        lambda{typ,i} = test.*lam;
        v_new{typ,i} = test.*v_1;
    end % typ
end % i
       