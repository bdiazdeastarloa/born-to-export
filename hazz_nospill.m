function [e_hazz]= hazz_nospill(a0,b0,mm)

% ------------------------------------------------------------------------- 
% Born to export 
%
% hazz_nospill
% Gives firms' posterior beliefs as a function of a, b, priors, 
% and known idiosyncratic productivity shock.
% Each cell contains values for all possible values of shocks, given (a,b).
%
% Written by James Tybout (2010).
% Edited and revised by Bernardo Diaz de Astarloa (2011-13).
% -------------------------------------------------------------------------

dim0   = mm.dim0;
dim2   = mm.dim2;
theta2 = mm.theta2;
amax   = mm.n_size;
bmax   = mm.n_size;
hazz   = mm.hazz;

e_hazz = cell(amax+1,bmax+1);

for a = 0:amax
    for b = 0:bmax
        % Get posterior beliefs;    
        jden = binopdf(a,a+b,hazz).*(ones(dim0,1)*betapdf(theta2,a0,b0)); 
        jden = jden./(sum(jden,2)*ones(1,dim2));
        e_hazz{a+1,b+1}  = sum(hazz.*jden,2);
    end
end

