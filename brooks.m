function [mtotex,mnumex,mavex,maxdur] = brooks(export,S,TT)

% -------------------------------------------------------------------------
% Born to export 
%
% brooks: generate Brooks tables (cohorts).
%
% Written by James Tybout (2010).
% -------------------------------------------------------------------------

start = zeros(S,TT);
cont = zeros(S,TT);
start(:,1) = export(:,1)>0;
start(:,2:TT) = (export(:,2:TT)>0).*(export(:,1:TT-1)==0);
cohort = ones(S,1)*(1:1:TT).*start;

for t = 2:TT;
cont(:,t) = (cont(:,t-1)+start(:,t-1).*cohort(:,t-1)).*(export(:,t)>0);
end;
cohort = cont + cohort;

maxdur = max(max((cohort>0).*(ones(S,1)*(1:1:TT)-cohort))) + 1; % duration, longest-lived cohort

totex = zeros(TT,maxdur);
numex = zeros(TT,maxdur);
avex  = zeros(TT,maxdur);
for t=1:TT-maxdur+1
    active = cohort(:,t:t+maxdur-1)==t;
    totex(t,:) = sum(active.*export(:,t:t+maxdur-1));
    numex(t,:) = sum(active);
    avex(t,:)  = totex(t,:)./max(numex(t,:),1) ;
end

burn = 10;
mtotex = mean(totex(burn:TT,:));
mnumex = mean(numex(burn:TT,:));
mavex  = mean(avex(burn:TT,:));