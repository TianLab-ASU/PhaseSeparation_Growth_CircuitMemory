function [t,X, N_dilute, N_droplet,concentration_dilute,concentration_droplet,concentration_ave,vt_rec]=ExactStoch_BatchGrowth_EndPoint(DefineReactions,tspan, IC)
%Implementation of Gillespie Exact Stochastic Algorithm.
% S = stoichiometry of C substrates in R reactions
% P = stoichiometry of products
% K = vector of reaction rates

global Vt Vt0 V_TF Omega Drop_Factor 

% INITIALIZE by calling the passed function to define the reactions:
if ~isa(DefineReactions, 'function_handle') ...
        disp('pass a function handle with your reaction defs'); return; end
X = IC;  % initialize the matrix of chemicals at time t
vt_rec=Vt; 
[S, P, K, N_dilute_SS, N_droplet_SS] = DefineReactions(0,X(end,:)); % defines reaction stoichiometry

N_dilute=N_dilute_SS;
N_droplet=N_droplet_SS;
V_droplet=N_droplet_SS*V_TF;
V_dilute=Vt-V_droplet; 
concentration_dilute=N_dilute_SS/V_dilute/Omega;
concentration_droplet=N_droplet_SS/V_droplet/Omega*Drop_Factor;
concentration_ave=X(end,:)/Vt/Omega; 

% [S, P, K] = DefineReactions(params); % defines reaction stoichiometry
% % check that size of inputs are correct:
% [R, C] = size(S);
% if size(S) ~= size(P); disp('reaction defs inconsistent'); return;end
% if length(IC) ~= C; disp('ICs inconsistent with reactions'); return; end
% if size(K,2) == 1; K = K'; end  % force column vector
% if length(K) ~= R; disp ('parameters inconsistent with reactions');return;end
if length(tspan) ~= 2; disp('Time vector: [0;maxtime] needed'); return;end
maxT = tspan(2);
tau = maxT/10000000; % initial check; time step is dynamic.
rand('state',sum(100000*clock));  %set the random number generator
t = 0;  %time
timeafterdivision=0;
% run loop until time is up or all chemicals are gone:

global Nmax1 lagphase growthrate
CellNum=1/64; 
DoublingTime=log(2)/growthrate;  

while (t<maxT && any(X(:)) ) 

    Vt=Vt0*(timeafterdivision/DoublingTime+1); 
 
    timeafterdivision=timeafterdivision+tau;
    if abs(timeafterdivision-DoublingTime)<DoublingTime/100 %cell division
        timeafterdivision=0;
        X(:)=ceil(X(:)/2);
        Vt=Vt/2; 
        CellNum=CellNum*2;
        GrowthRATE=growthrate*(1-CellNum/Nmax1)*(t>lagphase);   
        DoublingTime=log(2)/GrowthRATE; 
    end
 

    [S, P, K, N_dilute_SS, N_droplet_SS] = DefineReactions(t,X(:)); % defines reaction stoichiometry

    N_dilute=N_dilute_SS;
    N_droplet=N_droplet_SS; 
 
    vt_rec=Vt;
    V_droplet=N_droplet_SS*V_TF;
    V_dilute=Vt-V_droplet; 
    concentration_dilute=N_dilute_SS/V_dilute/Omega;
    concentration_droplet=N_droplet_SS/V_droplet/Omega*Drop_Factor;
    concentration_ave(:)=X(:)/Vt/Omega; 

    % check that size of inputs are correct:
    [R, C] = size(S);
    if size(S) ~= size(P); disp('reaction defs inconsistent'); return;end
    if length(IC) ~= C; disp('ICs inconsistent with reactions'); return; end
    if size(K,2) == 1; K = K'; end  % force column vector
    if length(K) ~= R; disp ('parameters inconsistent with reactions');return;end


    %step 1: Calculate a's (reaction rates given system state)
    a = K; 
    for r = 1:R
        for c = 1:C
            if S(r,c) == 1; a(r) = a(r)*X(c);
            elseif S(r,c) == 2
                a(r) = a(r)*X(c)*(X( c)-1)/2;
            elseif S(r,c) == 3
                a(r) = a(r)*X(c)*(X(c)-1)/2*(X(c)-2)/3;
            end
        end
    end
     
    a0 = sum(a); % a0 is the total rate of change of system
    if a0 == 0 % system can't change; finish and exit;
        X(:) = X(:);
        t =maxT;
        break;
    end

    %Step 2: calculate tau and r using random number generators
    % determine time of next reaction:
    p1  = rand;  tau = (1/a0)*log(1/p1);
    % determine which next reaction is:
    p2 = rand;
    for r=1:R
        if (sum(a(1:r)) >= p2*a0); break; end
    end
    %Step 3: carry out the reaction
    t = t + tau;  % t is time array; add last entry to it.
    %         nRC = nRC + 1       ; % nRC is number of reactions so far.
    X=X-S(r,:)+P(r,:);
    %         X(end,:)


end %end of while (t(end)<tspan(2) & any(X))
end % end of function

