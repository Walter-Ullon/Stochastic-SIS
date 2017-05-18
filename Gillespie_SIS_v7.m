%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Gillespie SIS v7   6/2016   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This version uses a statistical threshold to pulse the population once
% the threshold is met. Said threshold is comprised of the autocorrelation
% output supplied by the EWS analysis in R. Results from said analysis are
% saved in CSV files which are then read by this program and compared in
% order to ascertain whether threshold has been met. The resuldt from the
% perturbation trials are then saved into a CSV file for further analysis.


format long

for i=1:1000
    
    if i==1
        ExtinctionM = zeros(i,7);           % Matrix containing extinction information.
        thresh_counter = 0;                 % Counts the number of times gone extinct.
        e_counter = 0;
    end
    
SEED = i;
stream0=RandStream('mt19937ar','Seed',SEED);
RandStream.setGlobalStream(stream0);

FileName = ['Seed' num2str(SEED) 'Output.csv'];
ACFs = csvread(FileName,1,9);                           % read acf(1) from csv (vector)
time_index_ACFs = csvread(FileName,1,1,[1, 1, length(ACFs), 1]);     % read time index from csv

meanACFs = 0.82912;
minACFs  = 0.71253;
maxACFs = 0.90713;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Gillespie SIS  v7  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This version captures the extinction time, plus it counts the number of
% times the infected population has gone extinct.


%%%%%%%%%%%%%%%%%%
%%% Parameters %%%
%%%%%%%%%%%%%%%%%%
beta = 1000;                         % Contact rate.
gamma = 99.98;                       % Recovery rate.
mu = 0.02;                           % Birth/death rate (+/-).
R0 = beta/(gamma + mu);              % Reproductive number.

n = 10^7;                            %%%  Max number of iterations %%%

t = 0;                               %%%  Time  %%%

t_ext = 0;                          % Captures the extinction time.

%%%%%%%%%%%%%%%%%%%%
%%% Populations  %%%
%%%%%%%%%%%%%%%%%%%%
N = 40;                          %% Total Population (dN/dt = 0)
endemic_S = (N/R0);
endemic_I_analytic = N*(1 - (1/R0));

I_init = 2;                      %% Initial INFECTED population.
S_init = N - I_init;             %% Initial SUSCEPTIBLE population.


S = S_init;                      %%  SUSCEPTIBLE pop.
I = I_init;                      %%  INFECTED pop.


t_new(1) = 0;                   %% Time array (even steps).
S_new(1) = S_init;              %% SUSCEPTIBLE population array (time course).
I_new(1) = I_init;              %% INFECTED population array (time course).
        

%%%%%%%%%%%%%%
%%% STEP 1 %%%
%%%%%%%%%%%%%%
    
a = zeros(5, 1);              % Events array.
 
jump = 1;                     % Designated time step. %%
j = 1; 
k = 1;
thresh_met = false;           % determines if threshold is met
 
for index = 1:n
    %%%%%%%%%%%%%%
    %%% STEP 2 %%
    %%%%%%%%%%%%%%
    r1 = rand(1);
    r2 = rand(1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       Events        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a(1) = mu*N;                 % Birth ------->  (S, I) --> (S+1 , I)
    a(2) = mu*S;                 % Death(S) ---->  (S, I) --> (S-1 , I)
    a(3) = mu*I;                 % Death(I) ---->  (S, I) --> (S , I-1)
    a(4) = (beta*I*S/N);         % Infection --->  (S, I) --> (S-1 , I+1)
    a(5) = gamma*I;              % Recovery ---->  (S, I) --> (S+1 , I-1)
    
    a0 = sum(a);                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tau = log(1/r1)*(1/a0);          % Time increment.
    t = t + tau;
    
    st = r2*a0;                     % Stochastic time?

     if st <= a(1)
        S = S + 1;
    elseif st > a(1) && st <= (a(1) + a(2))
        S = S - 1;
    elseif st > (a(1) + a(2)) && st <= (a(1) + a(2) + a(3))
        I = I - 1;
    elseif st > (a(1) + a(2) + a(3)) && st <= (a(1) + a(2) + a(3) + a(4))
        S = S - 1;
        I = I + 1;
    else 
        I = I - 1;
        S = S + 1;
     end
     
     
     
    if t > jump
        t_new(j+1,1) = floor(t);
        S_new(j+1,1) = S;
        I_new(j+1,1) = I;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if thresh_met==false                      % ensures pulse happens only once
            if t_new(j) == time_index_ACFs(k)     % matches time index with csv time index
                if ACFs(k) >= maxACFs
                    r_factor = 0.5;               % removal factor      
                    p = floor(r_factor*I);        % Perturbation (scaled)
                    I = I - p;                            
                    N = S+I;                      % Update population post-pulse
                    thresh_met = true;            % theshold met!
                    info = [t_new(j), time_index_ACFs(k), ACFs(k), k, I, p+I, 0]; % save info to matrix
                end
                k = k+1;
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        if I<=0
             e_counter = e_counter + 1;
             t_new(j+1,1) = floor(t);
             S_new(j+1,1) = S;
             I_new(j+1,1) = I;
             t_ext = t_ext + t;
             info(7) = 1;
            break;
        end
        
       jump = jump + 1;
       j = j+1;
    end
end
        if thresh_met == true
             thresh_counter = thresh_counter+1;
             ExtinctionM(thresh_counter, 1) = SEED;
             ExtinctionM(thresh_counter, 2) = info(7);
             ExtinctionM(thresh_counter, 3) = info(2);
             ExtinctionM(thresh_counter, 5) = info(3);
             ExtinctionM(thresh_counter, 4) = floor(t_ext);
             ExtinctionM(thresh_counter, 6) = info(5);
             ExtinctionM(thresh_counter, 7) = info(6);
        end
        
        disp(SEED)
end

% csvwrite(['ExtinctionMatrix_' num2str(r_factor) '_LOW_.csv'], ExtinctionM);
% ML = csvread(['ExtinctionMatrix_' num2str(r_factor) '_LOW_.csv']);
% disp(['Removal factor = ' num2str(r_factor)]);
% disp(['Threshold met ' num2str(length(ML(:,1))) ' times (' num2str(length(ML(:,1))/i) ')']);
% disp(['Extinction reached ' num2str(sum(ML(:,2))) ' times']);
% disp(['Ratio = ' num2str(sum(ML(:,2))/length(ML(:,1)))]);


% csvwrite(['ExtinctionMatrix_' num2str(r_factor) '_HIGH_.csv'], ExtinctionM);
% MH = csvread(['ExtinctionMatrix_' num2str(r_factor) '_HIGH_.csv']);
% disp(['Removal factor = ' num2str(r_factor)]);
% disp(['Threshold met ' num2str(length(MH(:,1))) ' times (' num2str(length(MH(:,1))/i) ')']);
% disp(['Extinction reached ' num2str(sum(MH(:,2))) ' times']);
% disp(['Ratio = ' num2str(sum(MH(:,2))/length(MH(:,1)))]);


csvwrite(['ExtinctionMatrix_' num2str(r_factor) '_ABSHIGH_.csv'], ExtinctionM);
MH = csvread(['ExtinctionMatrix_' num2str(r_factor) '_ABSHIGH_.csv']);
disp(['Removal factor = ' num2str(r_factor)]);
disp(['Threshold met ' num2str(length(MH(:,1))) ' times (' num2str(length(MH(:,1))/i) ')']);
disp(['Extinction reached ' num2str(sum(MH(:,2))) ' times']);
disp(['Ratio = ' num2str(sum(MH(:,2))/length(MH(:,1)))]);


% csvwrite(['ExtinctionMatrix_' num2str(r_factor) '_ABSLOW_.csv'], ExtinctionM);
% ML = csvread(['ExtinctionMatrix_' num2str(r_factor) '_ABSLOW_.csv']);
% disp(['Removal factor = ' num2str(r_factor)]);
% disp(['Threshold met ' num2str(length(ML(:,1))) ' times (' num2str(length(ML(:,1))/i) ')']);
% disp(['Extinction reached ' num2str(sum(ML(:,2))) ' times']);
% disp(['Ratio = ' num2str(sum(ML(:,2))/length(ML(:,1)))]);
