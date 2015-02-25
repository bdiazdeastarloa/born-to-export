% -------------------------------------------------------------------------
% Born to export 
%
% setparams: sets main parameters of the model.
%
% Written by James Tybout (2010).
% Edited and revised by Bernardo Diaz de Astarloa (2011-13).
% -------------------------------------------------------------------------

mm = struct();

% Technology parameters
mm.r       = 0.05;              % Rate of time preference
mm.delta   = delta;             % Match separation rate (inverse of an integer)
mm.b       = beta;              % Cost function parameter
mm.scale_f = scale_f;           % Export profits scale parameter
mm.scale_h = scale_h;           % Domestic profits scale parameter
mm.f       = f;                 % Fixed cost of searching
mm.eta     = 5;                 % Demand elasticity 
mm.yr1     = 0.5;               % Fraction of sales in first year
mm.dud     = 0.2;               % Discount for matches that fail after first period


% Theta distributions.
% In the BTE case home is irrelevant and the global parameter does not play
% a role, so the parameters that matter are af and bf. The specification of
% three different thetas comes from an old version of EEJKT.
mm.ahf     = ahf;               % Beta distribution parameter (theta0)
mm.bhf     = bhf;               % Beta distribution parameter (theta0)
mm.ah      = ah;                % Beta distribution parameter (theta1)
mm.bh      = bh;                % Beta distribution parameter (theta1)
mm.af      = af;                % Beta distribution parameter (theta2)
mm.bf      = bf;                % Beta distribution parameter (theta2)

mm.gamma1  = ah/(ah+bh);        % True mean home theta1 dist
mm.gamma2  = mm.gamma1;         % True mean foreign theta2 dist
mm.gamma2p = 0.5*mm.gamma2;     % Prior mean of foreign theta2 dist
mm.alpha   = 0;                 % Weight on home success rate
    
mm.mean_zh = 0;                 % Mean product appeal (home) 
%mm.sig_zh   = sig_z_h;          % Std. dev. product appeal 
mm.F       = F;                 % Fixed cost of maintaining a client


% Costs 
% Cost function
mm.cost_fn   = @(lambda)(mm.b.*(lambda./(1.-lambda)) + mm.f.*(lambda>0)); 

% Shock processes
mm.mean_phi    = 0;             % Mean productivity level
mm.rho_phi     = rho_phi;                    % Root, productivity shock
mm.sig_eps_phi = 0.6839*sqrt(1-rho_phi^2);   % Std. dev. of innovation in productivity shock
%mm.rho_phi     = 0.74;          % Root, productivity shock (CJ's estimate)
%mm.sig_eps_phi = 0.46;          % Std. dev. of innovation in productivity shock (CJ's estimate)  

mm.mean_z      = 0;             % Mean product appeal (foreign)
mm.sig_eps_z   = sig_eps_z;     % Std. dev. of innovation in product appeal shock
mm.rho_z       = rho_z;         % Root of product appeal shock 
    
mm.mean_xh     = 0;             % Mean home macro state 
mm.rho_xh      = 0.961;         % Root of home macro shock
mm.sig_eps_xh  = 0.081;         % Std. dev. of innovation in home macro shock
    
mm.mean_xf     = 0;             % Mean foreign macro state
mm.rho_xf      = 0.953;         % Root of foreign macro shock
mm.sig_eps_xf  = 0.052;         % Std. dev. of innovation in foreign macro shock


% Discretization of state-space
mm.grid_length = 2.5;           % # of std. dev. from mean for discretization
mm.n_size      = 50;            % max # of informative signals per firm
mm.z_size      = 30;            % # of discretized buyer states
mm.zh_size     = mm.z_size;     % # of discretized home buyer states
mm.phi_size    = 30;            % # of different discretized profit shocks
mm.x_size      = 15;            % # of different discretized macro shocks  
mm.lambda_size = 40;            % # of possible effort levels
mm.theta_size  = 51;            % # of possible market potential values
if mm.alpha==0
    mm.dim0    = 1;             % # of possible theta0 values (common)
else                            % makes no sense to set >0 if they don't play a role
    mm.dim0    = 5;             % (they enlarge matrices unnecesarily)
end
mm.dim1        = 100;           % # of possible theta1 values (home)
mm.dim2        = 100;           % # of possible theta2 values (foreign)
    
mm.theta0      = 1/mm.dim0:1/mm.dim0:1;
mm.theta1      = 1/mm.dim1:1/mm.dim1:1;
mm.theta2      = 1/mm.dim2:1/mm.dim2:1;
theta0vec      = repmat(mm.theta0',1,mm.dim2);
theta2vec      = repmat(mm.theta2,mm.dim0,1);

if linear==1
    mm.hazz = (theta0vec*mm.alpha)+(theta2vec*(1-mm.alpha)); 
else
    mm.hazz = (theta0'.^mm.alpha)*(theta2.^(1-mm.alpha));
end

% Solution parameters
mm.v_tol  = 1e-4;               % Tolerance, value function iteration
mm.pi_tol = 1e-4;               % Tolerance, profit function

% Implied parameters
mm.beta      = 1/(1+mm.r);             % Discount factor
mm.beta_adj  = (1-mm.delta)/(1+mm.r);  % Match profit stream discount factor adjusted for separation risk
mm.match_dur = 1/mm.delta;             % Expected duration of a match (must be an integer)
