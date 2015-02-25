function MM = trajec_bte(mm,pi_tilda_h,pi_tilda_f,lambda_f,lambda_h,e_hazz,chi_c,chi_e)

% -------------------------------------------------------------------------
% Born to export 
%
% trajec_bte: simulates trajectories for a sample of firms and computes
% moments.
%
% Written by James Tybout (2010).
% Edited and revised by Bernardo Diaz de Astarloa (2011-14).
% -------------------------------------------------------------------------

global linear

%% Settings & Parameters

% Simulation settings.
trunc = 50;                   % # of observations to be discarded as burn-in
esT   = 100;                  % # of ergodic state periods to be simulated   
TS    = trunc + esT;          % Number of time periods to be simulated
T     = 25;                   % Time start for estimation after burn-in 
S     = 2000;                 % # of firms to simulate

% Assign parameters from structure.
delta   = mm.delta;     
scale_f = mm.scale_f; 
scale_h = mm.scale_h;
eta     = mm.eta;
alpha   = mm.alpha;
yr1     = mm.yr1;
dud     = mm.dud;
F       = mm.F;

mean_z      = mm.mean_z;    
rho_z       = mm.rho_z;
sig_eps_z   = mm.sig_eps_z;
mean_phi    = mm.mean_phi;     
rho_phi     = mm.rho_phi;      
sig_eps_phi = mm.sig_eps_phi; 
mean_xh     = mm.mean_xh;   
rho_xh      = mm.rho_xh;     
sig_eps_xh  = mm.sig_eps_xh; 
mean_xf     = mm.mean_xf;   
rho_xf      = mm.rho_xf;     
sig_eps_xf  = mm.sig_eps_xf; 

dim0   = mm.dim0;           % # of possible theta0 (global firm effect) values
dim1   = mm.dim1;           % # of possible theta1 (home firm effect) values
dim2   = mm.dim2;           % # of possible theta2 (foreign firm effect) values
ahf    = mm.ahf;            % home and foreign true theta0 parameter
bhf    = mm.bhf;            % home and foreign true theta0 parameter
ah      = mm.ah;
bh      = mm.bh;
theta0 = mm.theta0;         % vector of global firm effects
theta1 = mm.theta1;         % vector of home market firm effects
theta2 = mm.theta2;         % vector of foreign market firm effects

grid_length = mm.grid_length;  
n_size      = mm.n_size;        
phi_size    = mm.phi_size;     
x_size      = mm.x_size;    
z_size      = mm.z_size;


%% Generate product appeal probabilities and map them onto discrete grid

% No need to draw these every time.
% Comments out the relevant commands in the section below.

%{
yy0  = rand(S,1);               % draws for theta0
yy1  = rand(S,1);               % draws for theta1
yy2  = rand(S,1);               % draws for theta2

% Shocks
eps1 = normrnd(zeros(S,TS),1);  % idiosyncratic productivity
eps2 = normrnd(zeros(1,TS),1);  % idiosyncratic home macro shocks
eps3 = normrnd(zeros(1,TS),1);  % foreign idiosyncratic productivity
eps4 = normrnd(zeros(S,TS),1);  % idiosyncratic match-specific shocks
eps5 = normrnd(zeros(S,TS),1);  % idiosyncratic match-specific shocksS

drw1 = rand(S,TS);              % random client matching, foreign
drw2 = rand(S,TS);              % random client separations, foreign
drw3 = rand(S,TS);              % match success, foreign 

drw4 = rand(S,TS);              % random client matching, home
drw5 = rand(S,TS);              % random client separations, home
drw6 = rand(S,TS);              % match success , home
 
save('drw_thetas','yy0','yy1','yy2');
save('drw_eps','eps1','eps2','eps3','eps4','eps5');
save('drw_other','drw1','drw2','drw3','drw4','drw5','drw6');
%}

load drw_thetas
load drw_eps
load drw_other


%%  Grids for exogenous variables   
      
Z_f = tauchen_grid(sig_eps_z,z_size,mean_z,grid_length,rho_z);     
Z_h = tauchen_grid(sig_eps_z,z_size,mean_z,grid_length,rho_z); 
Phi = tauchen_grid(sig_eps_phi,phi_size,mean_phi,grid_length,rho_phi);       
X_h = tauchen_grid(sig_eps_xh,x_size,mean_xh,grid_length,rho_xh); 
X_f = tauchen_grid(sig_eps_xf,x_size,mean_xf,grid_length,rho_xf);
   

%% Generate acceptance rates for home and foreign market

th0 = betaincinv(yy0,ahf,bhf);
th1 = betaincinv(yy1,ah,bh);
th2 = betaincinv(yy2,ah,bh); 
th_draw = cat(2,th0,th1,th2);

th    = cell(3,1);
th{1} = theta0;
th{2} = theta1;
th{3} = theta2;

in_type = zeros(S,3);               % type index to fit to grid
theta = zeros(S,3);                 % put all thetas together in a matrix
dim   = cat(2,dim0,dim1,dim2);      % put all theta dimensions in a matrix

% Fit random draws to discrete grid
for j = 1:3 % index global/home/abroad
    for s = 1:S % index firms
        dif            = abs(th_draw(s,j)*ones(1,dim(j))-th{j}); 
        [~,in_type(s,j)] = min(dif);
        theta(s,j)     = th{j}(in_type(s,j));
    end
end

if linear==1
    % Linear specification
    mu_h = (theta(:,1)*alpha)+(theta(:,2)*(1-alpha));
    mu_f = (theta(:,1)*alpha)+(theta(:,3)*(1-alpha));
else
    % Cobb-Douglas specification
    mu_h = (theta(:,1).^alpha).*(theta(:,2).^(1-alpha));
    mu_f = (theta(:,1)^alpha).*(theta(:,3).^(1-alpha));
end


%% Generate exogenous trajectories and map them onto discrete grid

v_phi_shock = sig_eps_phi*eps1;     % idiosyncratic productivity shocks
v_h_shock   = sig_eps_xh*eps2;      % idiosyncratic home macro shocks
v_f_shock   = sig_eps_xf*eps3;      % foreign idiosyncratic productivity shocks
v_z_shock   = sig_eps_z*eps4;       % idiosyncratic match-specific shocks
v_zc_shock  = sig_eps_z*eps5;       % idiosyncratic match-specific shocks

[phi_s,in_phi] = dis_traj(S,TS,Phi,phi_size,v_phi_shock,mean_phi,rho_phi);
[x_f,in_xf]    = dis_traj(1,TS,X_f,x_size,v_f_shock,mean_xf,rho_xf);
[x_h,in_xh]    = dis_traj(1,TS,X_h,x_size,v_h_shock,mean_xh,rho_xh);
[z_f,~]        = dis_traj(S,TS,Z_f,z_size,v_z_shock,mean_z,rho_z);
[z_h,~]        = dis_traj(S,TS,Z_h,z_size,v_zc_shock,mean_z,rho_z);
 

%% Map signals onto search intensities and impute export trajectories

maxc     = 10;                  % maximum number of clients considered
hazard_f = zeros(S,TS);          % stores firm-specific hazards through time
cont     = zeros(S,TS);          % continue/exit indicator
entry    = zeros(S,TS);          % entry indicator
export   = zeros(S,TS);          % firm-specific total exports
clients  = zeros(S,TS);          % firm-specific client counts
dur      = zeros(S,TS);          % duration of matches begun in current period
initex   = zeros(S,TS);          % first-period exports of match begun in t by s
newex    = zeros(1,TS);          % number of new exporters each period
transit  = zeros(maxc,maxc,TS);  % transition densities, number of clients
muf_hat  = zeros(S,TS);          % theta-hat: posterior mean success rate in foreign
maxsig   = zeros(S,1);          % indicator for whether firm has reached learning max

% 3 dimensional arrays are more efficient than cells;
excum   = zeros(TS,TS,S);         % client-specific export trajectories, given firm
exrev_m = zeros(TS,TS,S);         % sales trajectory for a newly added buyer
ecount  = zeros(TS,TS,S);         % trach # of active relationships
match   = zeros(S,TS);            % track # of matches

succ = mu_f*ones(1,TS) > drw3(:,1:TS);    %#ok success! This is a (SxT) matrix
track_succ = zeros(S,2);                  % track succeses (:,1) and failures (:,2)
track_act  = zeros(S,TS);                 % track active firms     
track_succ_cum = zeros(S,2,TS);
track_succ_lag = track_succ;

% Load policy fn, value fn, cont and exit indicators
% Comment this out if running standalone, and check which spec is being
% simulated (spec is global above).
% load sol_bte_spec0;

% -------------------------------------------------------------------------
% Simulation starts
% -------------------------------------------------------------------------
for t=1:TS    
    for ss = 1:S
        if maxsig(ss)==1;  % stop updating firms' records after n_size successes or failures
             track_succ(ss,:) = track_succ_lag(ss,:);
        else
           if track_succ(ss,1) >= n_size-1 || track_succ(ss,2) >= n_size-1,  % 'maxed out firm-specific parameter'
               track_succ(ss,1) = min(track_succ(ss,1),n_size-1);
               track_succ(ss,2) = min(track_succ(ss,2),n_size-1);
               maxsig(ss) = 1;
           else
               track_succ(ss,1) = track_succ_lag(ss,1) + succ(ss,t);         % update successes
               track_succ(ss,2) = track_succ_lag(ss,2) + 1 - succ(ss,t);     % update failures
           end
        end

        % Entry and exit indicators
        cont(ss,t)  = chi_c{track_succ(ss,1)+1,track_succ(ss,2)+1,in_type(ss,1)}(in_phi(ss,t),in_xf(t));
        entry(ss,t) = chi_e{in_type(ss,1)}(in_phi(ss,t),in_xf(t));
        
        % Track active firms
        if t==1,
            track_act(ss,1) = entry(ss,1);
        else
            track_act(ss,t) = (1-track_act(ss,t-1))*entry(ss,t) + track_act(ss,t-1)*cont(ss,t);
        end
        
        % Update successes taking into account active status
        track_succ(ss,:) = track_act(ss,t)*track_succ(ss,:);
        
        % Compute hazard of a match (search policy)
        % Note that with alpha=0, the lambda_f{.,.,typ} should be equal for
        % all typ's, since global theta0 does not count.
        hazard_f(ss,t) = track_act(ss,t)*...
          lambda_f{track_succ(ss,1)+1,track_succ(ss,2)+1,in_type(ss,1)}(in_phi(ss,t),in_xf(t));  
        
        % Check whether there's a match
        match(ss,t) = hazard_f(ss,t) > drw1(ss,t);

        % Calculate revenues for matched firms
        if match(ss,t) == 1,
            exrev_m(t,t:TS,ss) = eta*scale_f*exp((eta-1)*phi_s(ss,t:TS) + x_f(t:TS) + z_f(ss,t:TS));
            exrev_m(t,t,ss) = ((1-succ(ss,t))*dud + succ(ss,t)*yr1)*exrev_m(t,t,ss);
            initex(ss,t)  = exrev_m(t,t,ss);
            for tt = t+1:TS
            % Truncate export streams that began in t if match dies in tt (hazard delta)
            % or continuation value is less than the fixed cost of maintaining client 
                dum2 = pi_tilda_f(in_phi(ss,tt),in_xf(tt))>F;
                dum3 = drw2(ss,tt)> delta ;
                flag1 = dum2*dum3*succ(ss,t);
                if flag1==0,
                    exrev_m(t,tt:TS,ss)=0; 
                    dur(ss,t) = tt-t; % record duration of the match that started in t
                    break
                end
            end
        end
        
        % Store exports, active matches
        excum(t,t:TS,ss) = exrev_m(t,t:TS,ss);             % trajectory for exports to the client met at t (if any)                                        
        ecount(:,t,ss)  = excum(:,t,ss) > 0;               % vector of dummies for active relationships, firm ss
        clients(ss,t)   = sum(ecount(:,t,ss));             % # of clients in year t, firm ss
        export(ss,t)    = sum(excum(:,t,ss));              % value of exports in year t, firm ss
        
        % Update prior
        muf_hat(ss,t)   = e_hazz{track_succ(ss,1)+1,track_succ(ss,2)+1}(in_type(ss,1)); 
    end
       
    track_succ_lag = track_succ;
    track_succ_cum(:,:,t) = track_succ;
end 

% Check for empty industry.
if isempty(find(entry, 1))==1
    display('ERROR: No firm ever enters! No moments to compute.');
    emptyflag = 1;
else
    emptyflag = 0;
end

%% Compute aggregates, indicators and moments

if emptyflag==0
    % ---------------------------------------------------------------------
    % Generate Brooks tables (cohorts)
    % ---------------------------------------------------------------------

    %{ 
    cuts = 10;
    % [avdur,qtile] = quant_dur(dur,initex,10);
    % [mtotex,~,mavex,maxdur] = brooks(export,S,TS);
    % 
    % if spec==1
    %     save brooks_spec1 mtotex mavex maxdur;      % sunk cost & scrap value
    % elseif spec==2
    %     save brooks_spec2 mtotex mavex maxdur;      % sunk cost & no scrap value
    % else
    %     save brooks_spec0 mtotex mavex maxdur;      % no sunk cost & no scrap value
    % end
    %}

    
    % ---------------------------------------------------------------------
    % Compute aggregates
    % ---------------------------------------------------------------------
    % Successes and failures
    agg_aa  = sum(track_succ_cum(:,1,:));
    agg_bb  = sum(track_succ_cum(:,2,:));
    active  = sum(track_act(:,:));

    % Continuing, entering and exporting firms.
    exdum   = export > 0;
    totexp   = sum(export);   
    totexdum = sum(exdum);

    % Normalized aggregates.
    % totexp = sum(export)./sum(export(:,1));   % normalize initial exports = 1
    % totcli = sum(clients)./sum(clients(:,1)); % normalizeinitial total clients = 1
    % totexdum = sum(exdum)./sum(exdum(:,1));
    % totcli   = sum(clients);

    % Frequency counts on numbers of clients.    
    % [nr,nc] = size(clients);
    % clivec = reshape(clients,nr*nc,1);
    % dexp = clivec; 

    
    % ---------------------------------------------------------------------
    % Survival rates by age profiles.
    % ---------------------------------------------------------------------
    %{ 
    durcount   = zeros(7,1);        % # of exporters of age k
    % totexp_age = zeros(7,1);        % total exports in age k
    % for ss = 1:S
    %     first_t = find(track_act(ss,:),1);              % period of entry
    %     if isempty(first_t) == 0 && first_t ~= 1            % don't count the first cohort
    %         for k = 1:7
    %             if first_t + k-1 <= TS 
    %                 if track_act(ss,first_t+k-1)>0        % check firm is active
    %                     durcount(k) = durcount(k) + 1;
    %                     totexp_age(k) = totexp_age(k) + export(ss,first_t+k-1); 
    %                 end
    %             end
    %         end
    %     end
    % end
    %}

    
    % ---------------------------------------------------------------------
    % Create transition densities, by period.
    % ---------------------------------------------------------------------
    for t = 1:TS-1
        for i=1:S        
            ii  = clients(i,t)+1;     % ii = 1 correpsonds to zero clients in t
            jj  = clients(i,t+1)+1;   % jj = 1 corresponds to zero clients in t+1
            iii = min(ii,maxc);
            jjj = min(jj,maxc);
            transit(iii,jjj,t) = transit(iii,jjj,t)+1; 
        end
    end

    for t=1:TS-1
        newex(t) = sum(transit(1,:,t)) - transit(1,1,t);
        for i=1:maxc
            temp = sum(transit(i,:,t),2);
            temp = max(temp,1);
            transit(i,:,t)=transit(i,:,t)./temp; % convert transitions to cond. probs.
        end
    end

    % mtran = sum(transit,3)/(TS-1);
    % vtran = cat(2,mtran(2,2:3),mtran(3,2:3)); % picks off key transition probs. for calibration

    % ---------------------------------------------------------------------
    % Calculate domestic sales trajectories
    % ---------------------------------------------------------------------
    %{
    % match_h  = zeros(S,TS);
    % hazard_h = zeros(S,TS);
    % sale_h   = zeros(S,TS);
    % succ     = mu_h*ones(1,TS)>drw6(:,1:TS);      %#ok success!
    % %lambda_h = policy_h(pi_tilda_h,mm);
    % 
    % for ss=1:S
    %     sales_sc = zeros(TS,TS);
    %     homerev = zeros(1,T);            % matrix of sales for matches at home
    %     for t=1:T
    %         hazard_h(ss,t) = lambda_h{in_type(ss,1),in_type(ss,2)}(in_phi(ss,t),in_xf(t));
    %         if drw4(ss,t)<hazard_h(ss,t),
    %             match_h(ss,t:TS) = match_h(ss,t:TS) + succ(ss,t);
    %             homerev(t:TS) = eta*scale_h*exp((eta-1)*phi_s(ss,t:TS) + x_h(t:TS)+(z_h(ss,t:TS)));
    %             homerev(t) = ((1-succ(ss,t))*dud + succ(ss,t)*yr1)*homerev(t);        
    %         end
    %         for tt = t+1:TS
    %             % Check exit conditions (fixed costs, exog. death)
    %             dum4  = pi_tilda_h(in_phi(ss,tt),in_xh(tt))>F;
    %             dum5  = drw5(ss,tt)>delta;
    %             flag2 = dum4*dum5*succ(ss,t);
    %             % Firm ss exits, break here
    %             if flag2==0,
    %                 homerev(tt:TS)=0;
    %                 break
    %             end
    %         end
    %         sales_sc(t,t:TS)= homerev(t:TS);
    %     end
    %     sale_h(ss,:)= sum(sales_sc); % value of domestic sales by firm and year
    % end
    %}

    
    % ---------------------------------------------------------------------
    % Construct moments to match to data
    % ---------------------------------------------------------------------
    t0 = trunc+1;
    Tm = trunc+T;

    % Match separation hazard rates.
    survive = zeros(T-1,1);
    hazrate = zeros(T-1,1);
    for t = 1:T-1;
        survive(t) = sum(sum(dur(:,t0:Tm)==t));
    end
    for t = 1:T-1;
        if sum(survive(t:T-1)) > 30;             % min # of cases for inference
            hazrate(t) = survive(t)/sum(survive(t:T-1));
        end
    end

    % Alternative match separation rate: using Brooks tables.
    %{
    % matchcount = zeros(T,T);
    % matchhaz  = zeros(5,1);
    % counter    = zeros(5,1);
    % for ss = 1:S
    %     for t = 1:T
    %         for tt = 1:(T-t+1)
    %             check = find(dur(ss,t0+t)==tt,1);            
    %             if check>0
    %                 for ttt = 1:tt
    %                     matchcount(t+ttt-1,t) = matchcount(t+ttt-1,t)+1;
    %                 end
    %             end
    %         end
    %     end
    % end
    % for t = 1:T
    %     for tt = 1:5
    %         if tt<=T-t
    %             surv = matchcount(t+tt,t)/matchcount(t+tt-1,t);
    %             if isnan(surv) == 1
    %                 surv = 0;
    %             end
    %             matchhaz(tt) = matchhaz(tt) + surv;
    %             if surv>0
    %                 counter(tt) = counter(tt)+1;
    %             end
    %         end
    %     end
    % end
    % hazrate_alt = 1-matchhaz./counter;
    %}

    haz1 = hazrate(1);
    temp = min(5,find(hazrate,1,'last'));
    haz2 = mean(hazrate(2:temp));
    matchhaz = [haz1;haz2];

    % ML estimate of log-normal distribution parameters for sales.
    exptrunc = export(:,t0:Tm);
    xx_ind   = exptrunc>0;                                     
    xx       = reshape(exptrunc(xx_ind),[],1);
    if size(xx,1)>0
        [b,~] = lognfit(xx);
        xdist = [b(1);b(2)];
    else
        xdist = [100;100];
        display('ERROR: Nothing to regress moms.m exports distribution');
    end

else
    matchhaz = [100;100];
    xdist    = [100;100];
end

% Put moments together.
MM = cat(1,xdist,matchhaz);


%% Save trajectories for figures

% Comment this out for estimation.
% Check periods to select (i.e. take into account burn-in).
%{
xr_traj = x_f(1:T)- mean(x_f(1:T));    % normalized foreign demand shifter 
save traj_par xr_traj mu_h mu_f;

global spec
if spec==1
    save traj_spec1 hazard_f muf_hat export totexdum active totexp agg_aa agg_bb;      % sunk cost & scrap value
elseif spec==2
    save traj_spec2 hazard_f muf_hat export totexdum active totexp agg_aa agg_bb;      % sunk cost & no scrap value
else
    save traj_spec0 hazard_f muf_hat export totexdum active totexp agg_aa agg_bb;      % no sunk cost & no scrap value
end
%}
