clear
format long

stream0=RandStream('mt19937ar','Seed',6);
RandStream.setGlobalStream(stream0);

fileName = 'SIS_seed6-N30.csv';

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Gillespie SIS  v3  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This version captures the extinction time, plus it counts the number of
% times the infected population has gone extinct.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Steady states                           %
%  (N, 0) -------------------------------->  Disease free.  %
%  ((N/R0), N(1 - (1/R0))----------------->  Endemic.       %
%                                                           %
%     S_dot = -(beta*I(i)*S(i)/N) + mu*I(i) + gamma*I(i);   %
%     I_dot = (beta*I(i)*S(i)/N) - gamma*I(i) -  mu*I(i);   %
%                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
%%% Parameters %%%
%%%%%%%%%%%%%%%%%%
beta = 1000;                         % Contact rate.
gamma = 99.98;                       % Recovery rate.
mu = 0.02;                           % Birth/death rate (+/-).
R0 = beta/(gamma + mu);              % Reproductive number.

n = 10^10;                            %%%  Max number of iterations %%%

t = 0;                               %%%  Time  %%%

t_ext = 0;                          %% Captures the extinction time.
e_counter = 0;                      %% Counts the number of times gone extinct.

%%%%%%%%%%%%%%%%%%%%
%%% Populations  %%%
%%%%%%%%%%%%%%%%%%%%
N = 30;                          %% Total Population (dN/dt = 0)
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
    
a = zeros(5, 1);                %% Events array.
 
jump = 1;                     %% Designated time step. %%

j = 1; 
 
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
        t_new(j+1) = floor(t);
        S_new(j+1) = S;
        I_new(j+1) = I;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                          
%          if t_new(j)== 172
%          p = 28;                           % Perturbation
%          I = I - p;                            
%          N = S+I;                          % Update population post-pulse
%          end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        if I<=0
             e_counter = e_counter + 1;
             t_ext = t_ext + t;
             
             t_new(j+1) = floor(t);
             S_new(j+1) = S;
             I_new(j+1) = I;
            break;
        end
        
       jump = jump + 1;
       j = j+1;
    end

end

M = [t_new; I_new]';                          % matrix of time vs. pop (for CSV)
csvwrite(fileName, M);                        % write to CSV


hold on
plot(t_new, S_new);
plot(t_new, I_new);
% plot(365,33,'ro', 'MarkerSize',12, 'Linewidth', 1.5)
% plot(172,33,'bo', 'MarkerSize',12, 'Linewidth', 1.5)
% axis([0,800, 0, 50])
% title('Gillespie SIS v3');
xlabel('Time')                              % t-axis label
ylabel('Population')                        % N-axis label
legend('Susceptible','Infected', 'low resilience');

endemic_I_numerical = mean(I_new);
