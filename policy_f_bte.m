function [lambda,e_hazz,v_new,chi_c,chi_e,flag] = policy_f_bte(pi_tilda,p_phi,p_x,mm)

% -------------------------------------------------------------------------
% Born to export 
%
% policy_f_bte: solves for the foreign search intensity policy function, using
% FOCs.
%
% Written by James Tybout (2010).
% Edited and revised by Bernardo Diaz de Astarloa (2011-13).
% -------------------------------------------------------------------------

%% Assign parameters

v_tol     = mm.v_tol;
r         = mm.r;
b         = mm.b;
cost_fn   = mm.cost_fn;
ecost     = mm.ecost;
scrap     = mm.scrap;
n_size    = mm.n_size;
x_size    = mm.x_size;
phi_size  = mm.phi_size;
dim0      = mm.dim0;
af        = mm.af;
bf        = mm.bf;
    
v_old  = cell(n_size,n_size,dim0);
v_new  = cell(n_size,n_size,dim0);
lambda = cell(n_size,n_size,dim0);  
chi_c  = cell(n_size,n_size,dim0);  
chi_e  = cell(1,dim0);  


%% Compute value and policy functions
    
% Expected hazard rate 
% Each cell contains values for all possible theta0's, given a and b
e_hazz = hazz_nospill(af,bf,mm); 

for typ = 1:dim0
    for i = 1:n_size
        for j = 1:n_size
            t_hat = e_hazz{i,j}(typ);
            v_old{i,j,typ} = t_hat*pi_tilda;
            v_new{i,j,typ} = t_hat*pi_tilda;
            lambda{i,j,typ} = zeros(phi_size,x_size);
        end           
    end
end
el_count = (n_size^2)*phi_size*x_size*dim0;
m_old = zeros(el_count,1);

eps = 1;
iter = 1;
while eps > v_tol
  for typ=1:dim0;     
      for i = 1:n_size-1       % i indexes successes, plus 1        
         for j = 1:n_size-1    % j indexes failures, plus 1
             
             Ev11 = p_phi*max(v_old{i,j,typ},scrap)*p_x';     % no matches
             Ev12 = p_phi*max(v_old{i,j+1,typ},scrap)*p_x';   % no success
             Ev21 = p_phi*max(v_old{i+1,j+1,typ},scrap)*p_x'; % success

             t_hat = e_hazz{i,j}(typ);      
             dvdl = t_hat*(pi_tilda + (Ev21-Ev11)/(1+r)) + (1-t_hat)*(Ev12-Ev11)/(1+r);
             dvdl = max(dvdl,10^(-5));

             % Solve FOC to get search policy
             lam = max(1 - (b./dvdl).^(0.5),0) ;              
             % Compute value of searching        
             v_1 = -cost_fn(lam) + lam.*dvdl + Ev11/(1+r);
             % Compute value of not searching
             v_0 = Ev11/(1+r);                              
             test = v_1 >= v_0 ;

             lambda{i,j,typ} = test.*lam;
             v_new{i,j,typ} = test.*v_1 + (1-test).*v_0;
             % Continue or exit decision
             chi_c{i,j,typ} = v_new{i,j,typ} >= scrap;   
             % Entry decision
             if i == 1 && j == 1
                 chi_e{typ} = v_new{i,j,typ} >= ecost;   
             end
         end % j
         % Take care of boundaries
         lambda{i,n_size,typ} = lambda{i,n_size-1,typ};  
         v_new{i,n_size,typ}  = v_new{i,n_size-1,typ};
         chi_c{i,n_size,typ}  = v_new{i,n_size-1,typ} >= scrap; 
      end % i
      
      % Take care of boundaries (max number of successes)
      for j = 1:n_size     
          lambda{n_size,j,typ} = lambda{n_size-1,j,typ};
          v_new{n_size,j,typ}  = v_new{n_size-1,j,typ};
          chi_c{n_size,j,typ}  = v_new{n_size-1,j,typ} >= scrap; 
      end
  end % typ
  
  m_new = cell2mat(v_new);
  m_new = reshape(m_new,el_count,1);       
  eps   = norm(m_new-m_old)/norm(m_new);
  m_old = m_new;
  v_old = v_new;
            
  iter = iter + 1;
  if iter > 200
      fprintf('\n');
      fprintf('Max, number of iterations reached.\n'); 
      fprintf('Exiting solution routine at tol=%6.5f.\n',eps);
      flag = 0;
      break       
  end
  
flag = 1;  
end
