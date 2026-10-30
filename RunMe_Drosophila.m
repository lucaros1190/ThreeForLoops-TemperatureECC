
% RunME ODE script - Drosophila suzukii with temperature-dependent rates

% Created by Luca Rossini on 13 September 2026
% Contributors: Luca Rossini
% Last update 21 September 2026
% e-mail: l.rossini@reading.ac.uk

% INSTRUCTIONS:
%
% RunMe_Drosophila.m: it is the runme file, the executable script to run
% simulations. Keep it minimal!
%
% Parameters.m: contains all the model parameters and can be used to set
% the initial conditions
% 
% Functions.m: contains, in a static class, all the functions involved
% needed by the problem. Use this file to add any additional functions!
%
% Plots.m: contains only plots to show the results
%
% TemperatureInput.xlsx: this file contains the temperature with daily
% values. The simulation time DEPENDS on this array. If you want to extend
% the simulation time, you should extend this array with more temperature
% values.
%
% Computational complexity disclaimer: This increasing the number of
% treatments or the temporal horizon of the simulation will increase
% exponentially the computational time! 


%% Compute the simulation time

tic

%% Clearing the workspace before the beginning of the simulation

clear
clc
close all


%% Start parallel process - This is necessary to reduce the computation 
 % time required by the problem.

if isempty(gcp('nocreate'))
    parpool;
end


%% Load the parameters and other inputs

run("Parameters.m")

    % Load the functions from Functions.m

Fun = Functions;


%% Pre-compute all temperature-dependent rates.
 % This operation is necessary to not compute the rates within the ODE
 % solver at every iteration. As the horizon is defined and we have a map
 % of the daily temperatures, we can precompute the rates and store them
 % into a database. Then, the solver will just pick the corresponding value
 % from that database.

Rates = Fun.PrecomputeRates(DailyTemp, FertPar, MortPar_Egg, MortPar_L1, ...
                            MortPar_L2, MortPar_L3, DevPar_Egg, DevPar_L1, ...
                            DevPar_L2, DevPar_L3, DevPar_P, DevPar_AM, ...
                            DevPar_AF);


%% Solve the equation - ODE system open loop
 % We solve the ODE system in the case of ZERO treatments, corresponding to
 % a free evolution of the system. This solution acts as a ground truth for
 % comparisons.

    % Define the time span

t_span = load("Parameters.mat", "t_span");
t_span = cell2mat(struct2cell(t_span));

    % Note: the initial condition can be set in 'Parameters.m'

TreatmentTimes = [];

[timeODE, solODE] = ode45(@(t,y) Fun.ODEfun(t, y, SR, Rates, FruitAbortion, ...
                                         FruitPredation, TreatmentTimes, ...
                                         M_tratt, TrDecay, T_cover), ...
                                         t_span, InitCond_ODE);

    % Compute the yield (objective function) associated with the case of
    % no treatments. 

FruitsHarvest_0 = max(0, solODE(end,1));

J0 = -F_price * FruitsHarvest_0;

    % Print the results of the no-treatment strategy

fprintf('\nNo-treatment strategy:\n');
fprintf('Fruits at harvest: %.2f\n', FruitsHarvest_0);
fprintf('J = %.4f\n', J0);


%% Optimisation - Closed loop.
 % This part of the code deals with the optimal treatment strategy. We
 % hereafter compute all the possible scenarios, to find the most
 % economically convenient.

    % Define some storage factors and keeps J_0 as first best global
    % objective function.

BestGlobal_J = J0;
BestGlobal_N = 0;
BestGlobal_ts = [];

    % Store the optimal strategy for each number of treatments

BestJ_all  = NaN(1, N_tratt_max);
BestTs_all = cell(1, N_tratt_max);

    % Computes multiple treatments (in order: 1, 2, ..., N_tratt_max)

for N_tratt = 1:N_tratt_max

    fprintf('\nTesting %d treatments...\n', N_tratt);

    BestN_J = Inf;
    BestN_ts = [];

        % Identify the first possible treatment

    latestFirst = T_harvest ...
                  - T_preharvest ...
                  - (N_tratt - 1) * T_min;

    % Sores the possible first treatment times

    FirstTimes = 0:latestFirst;

    nFirst = length(FirstTimes);

    % Storage for parallel results: the parallel process assigns different
    % combinations of treatment times to each worker, then it stores the
    % objective function and the timing of the optimal treatments

    J_values = Inf(1, nFirst);
    ts_values = cell(1, nFirst);

    % Makes the parallel process

    parfor idx = 1:nFirst

        i = FirstTimes(idx);

        [J_values(idx), ts_values{idx}] = ...
            Fun.CostFunEvaluator_VarTemp_NT(i, N_tratt, Fun, InitCond_ODE, ...
            T_harvest, T_min, T_preharvest, M_tratt, TrDecay, T_cover, ...
            P_tratt, F_price, FruitAbortion, FruitPredation, SR, Rates);

    end

    % Wraps everything up: Find the best result for this number of 
    % treatments

    [BestN_J, best_idx] = min(J_values);

    BestN_ts = ts_values{best_idx};
    
    % Store optimum for this number of treatments

    BestJ_all(N_tratt)  = BestN_J;
    BestTs_all{N_tratt} = BestN_ts;

    % Prints the results on the screen

    fprintf('Best J for N = %d: %.4f\n', N_tratt, BestN_J);

    fprintf('Treatment times: ');
    fprintf('%d ', BestN_ts);
    fprintf('\n');

    % Updates the objective functions to compare the optimal timing for the
    % current number of treatments with the other number of treatments. In
    % this way we identify which is the most convenient number of
    % treatments and their optimal timing of application.

    if BestN_J < BestGlobal_J

        BestGlobal_J = BestN_J;
        BestGlobal_N = N_tratt;
        BestGlobal_ts = BestN_ts;

    end

end


%% Simulate optimal strategy: This is for plotting and presentation of the 
 % results purposes

TreatmentTimes = BestGlobal_ts;

[timeOPT, solOPT] = ode45(@(t,y) Fun.ODEfun(t, y, SR, Rates, FruitAbortion, ...
                    FruitPredation, TreatmentTimes, M_tratt, TrDecay, ...
                    T_cover), [0 T_harvest], InitCond_ODE);

% Print results on the screen

fprintf('\n========================================\n');
fprintf('OPTIMAL STRATEGY\n');
fprintf('========================================\n');

fprintf('Number of treatments: %d\n', BestGlobal_N);

fprintf('Treatment times: ');
fprintf('%d ', BestGlobal_ts);
fprintf('\n');

fprintf('Fruits at harvest: %.2f\n', solOPT(end,1));

fprintf('Treatment cost: %.2f\n', ...
        BestGlobal_N * P_tratt);

fprintf('Fruit value: %.2f\n', ...
        F_price * solOPT(end,1));

fprintf('Profit: %.2f\n', ...
        -BestGlobal_J);

fprintf('========================================\n');



%% Eigenvalue analysis: compute the dominant eigenvalues over time to see 
 % when the system is switching from stability to instability and the other
 % way round. This is essential to understand the connections between the
 % treatments and their connection with temperature changes.

    % Fine time grid for computation purposes

timeEig = linspace(0, T_harvest, 2000);

nEig = length(timeEig);

    % Define storage variables

Eig_NoTreat = zeros(8, nEig);
Eig_Optimal = zeros(8, nEig);

LambdaMax_NoTreat = zeros(1, nEig);
LambdaMax_Optimal = zeros(1, nEig);


    % Compute eigenvalues and dominant real parts over the time grid

for j = 1:nEig

    t = timeEig(j);

        % Case with 0 treatments

    A0 = Fun.JacobianInsect( t, SR, Rates, [], M_tratt, TrDecay, T_cover);

    lambda0 = eig(A0);
    Eig_NoTreat(:,j) = lambda0;
    LambdaMax_NoTreat(j) = max(real(lambda0));

        % Case of the optimal treatment
   
    AT = Fun.JacobianInsect(t, SR, Rates, BestGlobal_ts, M_tratt, TrDecay, ...
                            T_cover);

    lambdaT = eig(AT);
    Eig_Optimal(:,j) = lambdaT;
    LambdaMax_Optimal(j) = max(real(lambdaT));

end



%% Plot the results calling the 'Plots.m' script

run('Plots.m')


%% End the calculation of the simulation time

toc

