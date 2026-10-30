
% Definition of the functions contained in the RunMe_Drosophila.m 
% file

% Created by Luca Rossini on 3 November 2023
% Contributors: Luca Rossini
% Last update 17 September 2026
% E-mail: l.rossini@reading.ac.uk


%% Definition of a static class where all the functions are stored

classdef Functions

    methods (Static) % Insert any function in between "methods" and "end"


    %% Function 1: Temperature reader

     % Reads the temperatures from the Excel file associated with the
     % dataset. This file will determine the time span of the simulation,
     % as well. If you wish to extend your simulation time span, you should
     % provide additional temperatures!

        function Temp = TempFunction(time, TempArray)
            
            % This if/else statement is needed because there is a problem
            % when time is 0, because it goes to read the array in the 0
            % position while Matlab starts to count from 1

            if time <= 0
                Temp = TempArray(1);
            else
                Temp = TempArray(time);
            end 
        end


    %% Function 2: Fertility rate function

     % Temperature-dependent egg production rate

        function B = FertRate(Par, Temp)

            alpha = Par(1); 
            gamma = Par(2); 
            lambda = Par(3); 
            delta = Par(4); 
            tau = Par(5);
            
            B = alpha * ((gamma + 1)/(pi * lambda^(2 * gamma + 2)) * ...
                   (lambda^2 - ((Temp - tau)^2 + delta^2))^gamma);
            
        end


    %% Function 3: Logan development rate function

     % Temperature-dependent development rate function. To see which stages
     % are described from this equation, check the 'Parameters.m' file.

        function G = Log(Par, Temp)
        
            G = Par(1) * (exp(Par(2) * Temp) - exp(Par(2) * Par(3) - ...
                ((Par(3) - Temp) / Par(4))));
    
            if G > 0
                G = Par(1) * (exp(Par(2) * Temp) - exp(Par(2) * Par(3) - ...
                ((Par(3) - Temp) / Par(4))));
            else
                G = 0;
            end
            
        end


    %% Function 4: Sharpe and De Michele development rate function

     % Temperature-dependent development rate function. To see which stages
     % are described from this equation, check the 'Parameters.m' file.

        function G = SDM(Par, Temp)

                T = Temp + 273.15; % Temperature in Kelvin for this function!
    
                if T <= 0
                    T = 1;
                else
                    T = Temp;
                end

                G = (T * (exp(Par(1) - (Par(2) / T)))) / ((1 + exp(Par(3) - ...
                    (Par(4) / T)) + exp(Par(5) - (Par(6) / T))));

                if G > 0
                    (T * (exp(Par(1) - (Par(2) / T)))) / ((1 + exp(Par(3) - ...
                    (Par(4) / T)) + exp(Par(5) - (Par(6) / T))));
                else
                    G = 0;
                end
                
        end


    %% Function 5: Lactin-1 Development rate function

     % Temperature-dependent development rate function. To see which stages
     % are described from this equation, check the 'Parameters.m' file.

        function G = LacOne(Par, Temp)

            G = exp(Par(1) * Temp) - exp(Par(1) * Par(2) - ( (Par(2) - ...
                 Temp) / Par(3) ) );

            if G > 0
                G = exp(Par(1) * Temp) - exp(Par(1) * Par(2) - ( (Par(2) - ...
                 Temp) / Par(3) ) );
            else
                G = 0;
            end

        end


    %% Function 6: Briere development rate function

     % Temperature-dependent development rate function. To see which stages
     % are described from this equation, check the 'Parameters.m' file.

        function G = Bri(Par, Temp)

            a = Par(1);
            T_L = Par(2);
            T_M = Par(3);
            m = Par(4);

            if Temp < T_L || Temp > T_M
                G = 0;
            else
                G = a * (Temp - T_L) * (T_M - Temp)^(1/m);
            end

        end


    %% Function 7: Mortality rate function

     % Temperature-dependent development rate function. To see which stages
     % are described from this equation, check the 'Parameters.m' file.

        function MP = MortRate(Par, Temp)

            MP = exp(Par(1) * Temp^2 + Par(2) * Temp + Par(3));

            if MP >= 0
                MP = exp(Par(1) * Temp^2 + Par(2) * Temp + Par(3));
            else
                MP = 0;
            end
            
        end
   
        
    %% Function 8: Pre-calculation of temperature-dependent rates

     % As the time horizon of the simulation is already known, it is more
     % efficient to pre-compute all the rates at the beginning of the
     % simulation and then pass the corresponding values to the ODE45
     % algorithm. The alternative is to compute the rates at each iteration
     % of ODE45, but it is extremely computationally expensive!
        
        function Rates = PrecomputeRates(TempArray, FertPar, MortPar_Egg, ...
                                         MortPar_L1, MortPar_L2, MortPar_L3, ...
                                         DevPar_Egg, DevPar_L1, DevPar_L2, ...
                                         DevPar_L3, DevPar_P, DevPar_AM, ...
                                         DevPar_AF)

            F = Functions;
            nDays = length(TempArray);

        % Prepare the databases to store the results

            Rates.Fert    = zeros(nDays,1);

            Rates.MortEgg = zeros(nDays,1);
            Rates.MortL1  = zeros(nDays,1);
            Rates.MortL2  = zeros(nDays,1);
            Rates.MortL3  = zeros(nDays,1);

            Rates.DevEgg  = zeros(nDays,1);
            Rates.DevL1   = zeros(nDays,1);
            Rates.DevL2   = zeros(nDays,1);
            Rates.DevL3   = zeros(nDays,1);
            Rates.DevP    = zeros(nDays,1);

            Rates.DevAM   = zeros(nDays,1);
            Rates.DevAF   = zeros(nDays,1);


        % Compute rates and store them into the corresponding storage
        % points

            for d = 1:nDays

                Temp = TempArray(d);

                Rates.Fert(d) = F.FertRate(FertPar, Temp);

                Rates.MortEgg(d) = F.MortRate(MortPar_Egg, Temp);
                Rates.MortL1(d)  = F.MortRate(MortPar_L1, Temp);
                Rates.MortL2(d)  = F.MortRate(MortPar_L2, Temp);
                Rates.MortL3(d)  = F.MortRate(MortPar_L3, Temp);

                Rates.DevEgg(d) = F.SDM(DevPar_Egg, Temp);
                Rates.DevL1(d)  = F.Bri(DevPar_L1, Temp);
                Rates.DevL2(d)  = F.Log(DevPar_L2, Temp);
                Rates.DevL3(d)  = F.LacOne(DevPar_L3, Temp);
                Rates.DevP(d)   = F.Bri(DevPar_P, Temp);

                Rates.DevAM(d)  = F.LacOne(DevPar_AM, Temp);
                Rates.DevAF(d)  = F.LacOne(DevPar_AF, Temp);

                end

            end


    %% Function 9: compute the Jacobian matrix to support eigenvalue 
     % analysis

     % This function computes the Jacobian of the system, but it does not
     % include the plant, as it is not relevant for the purpose of this
     % study. The Jacobian can be hand-computed from the system and
     % inserted in this function, should you need to adapt this code to
     % other cases.
        
        function A = JacobianInsect(t, SR, Rates, TreatmentTimes, ...
                                    M_tratt, TrDecay, T_cover)

            % Computes the daily temperature-dependent rates

            day = max(1, floor(t) + 1);
            day = min(day, length(Rates.Fert));

            Fert   = Rates.Fert(day);

            DevEgg = Rates.DevEgg(day);
            DevL1  = Rates.DevL1(day);
            DevL2  = Rates.DevL2(day);
            DevL3  = Rates.DevL3(day);
            DevP   = Rates.DevP(day);

            DevAM  = Rates.DevAM(day);
            DevAF  = Rates.DevAF(day);

            MortEgg = Rates.MortEgg(day);
            MortL1  = Rates.MortL1(day);
            MortL2  = Rates.MortL2(day);
            MortL3  = Rates.MortL3(day);

            MortP = 0.01;

            % Compute daily mortality associated with treatments

            M_temp = 0;

            for k = 1:length(TreatmentTimes)

                dt = t - TreatmentTimes(k);

                if dt >= 0 && dt < T_cover

                    M_temp = M_temp + M_tratt * TrDecay^dt;

                end

            end

            M_temp = min(M_temp, 1);

            % Compute the Jacobian matrix. The state order is the
            % following:

            % 1 = Egg
            % 2 = L1
            % 3 = L2
            % 4 = L3
            % 5 = Pupa
            % 6 = Adult male
            % 7 = Non-mated female
            % 8 = Mated female

            A = zeros(8,8);

                % Egg

            A(1,1) = -(DevEgg + MortEgg);
            A(1,8) = DevAF * Fert;

                % L1

            A(2,1) = DevEgg;
            A(2,2) = -(DevL1 + MortL1);

                % L2

            A(3,2) = DevL1;
            A(3,3) = -(DevL2 + MortL2);

                % L3

            A(4,3) = DevL2;
            A(4,4) = -(DevL3 + MortL3);

                % Pupa

            A(5,4) = DevL3;
            A(5,5) = -(DevP + MortP);

                % Adult males

            A(6,5) = (1 - SR) * DevP;
            A(6,6) = -(DevAM + M_temp);

                % Non-mated females

            A(7,5) = SR * DevP;
            A(7,7) = -(1 + M_temp);

                % Mated females

            A(8,7) = 1 - DevAF;
            A(8,8) = -(DevAF + M_temp);

        end


    %% Function 10: Definition of the ODE system

     % This function contains the definition of the ODE system and served
     % the ODE45 solver. Should you need to adapt the code, here is where
     % you should write your model!

        function dydt_ODE = ODEfun(t, y, SR, Rates, FruitAbortion, ...
                                   FruitPredation, TreatmentTimes, M_tratt, ...
                                   TrDecay, T_cover)

            % Takes the daily temperature-dependent rate previously
            % computed and stored

            day = max(1, floor(t) + 1);
            day = min(day, length(Rates.Fert));

            Fert   = Rates.Fert(day);

            DevEgg = Rates.DevEgg(day);
            DevL1  = Rates.DevL1(day);
            DevL2  = Rates.DevL2(day);
            DevL3  = Rates.DevL3(day);
            DevP   = Rates.DevP(day);

            DevAM  = Rates.DevAM(day);
            DevAF  = Rates.DevAF(day);

            MortEgg = Rates.MortEgg(day);
            MortL1 = Rates.MortL1(day);
            MortL2 = Rates.MortL2(day);
            MortL3 = Rates.MortL3(day);
            MortP = 0.01;

            % Compute the effect of the treatments at time t

            M_temp = 0;

            for k = 1:length(TreatmentTimes)

                dt = t - TreatmentTimes(k);

                if dt >= 0 && dt < T_cover

                    M_temp = M_temp + M_tratt * TrDecay^dt;

                end

            end

            % Prevent treatment mortality from exceeding 100%
            
            M_temp = min(M_temp, 1);

            % List of the state variables associated with the life stages

            % y(1) = Fruits
            % y(2) = Eggs
            % y(3) = L1
            % y(4) = L2
            % y(5) = L3
            % y(6) = Pupa
            % y(7) = Adult males
            % y(8) = Non-mated females
            % y(9) = Mated females

                % Fruits

            Fruits = - FruitAbortion * y(1) ...
                     - FruitPredation * DevAF * Fert * y(9);

                % Insect population

            Egg = DevAF * Fert * y(9) - DevEgg * y(2) - MortEgg * y(2);


            L1 = DevEgg * y(2) - DevL1 * y(3) - MortL1 * y(3);


            L2 = DevL1 * y(3) - DevL2 * y(4) - MortL2 * y(4);


            L3 = DevL2 * y(4) - DevL3 * y(5) - MortL3 * y(5);


            P = DevL3 * y(5) - DevP * y(6) - MortP * y(6);

                % Adults + treatments

            Am = (1 - SR) * DevP * y(6) - DevAM * y(7) - M_temp * y(7);


            Anmf = SR * DevP * y(6) - y(8) - M_temp * y(8);


            Amf = y(8) - DevAF * y(8) - DevAF * y(9) - M_temp * y(9);

                % Output array

            dydt_ODE = [Fruits; Egg; L1; L2; L3; P; Am; Anmf; Amf];

                end


    %% Function 11: Optimiser - Brute force approach

     % This function is the core of the optimisation problem, as it will
     % evaluate the different combinations of treatments and will take into
     % account the explored and unexplored combinations. The time
     % combinations that fall in any coverage time will be automatically
     % skipped, saving computational power.

        function [best_J, best_ts] = CostFunEvaluator_VarTemp_NT(i, N_tratt, ...
                                         Fun, InitCond_ODE, T_harvest, ...
                                         T_min, T_preharvest, M_tratt, ...
                                         TrDecay, T_cover, P_tratt, F_price, ...
                                         FruitAbortion, FruitPredation, SR, ...
                                         Rates)

            best_J = Inf;
            best_ts = [];

            current_comb = zeros(1, N_tratt);
            current_comb(1) = i;

            recursive_build(1, i + T_min);

            % Define an internal function to deal with recursive 
            % combination generator

            function recursive_build(n_selected, current_time)

                % Complete treatment combination

                if n_selected == N_tratt

                    TreatmentTimes = current_comb;

                % Calling the ODE solver and solves the ODE
                    [~, sol] = ode45(@(t,y) Fun.ODEfun(t, y, SR, Rates, ...
                                FruitAbortion, FruitPredation, ...
                                TreatmentTimes, M_tratt, TrDecay, T_cover), ...
                                [0 T_harvest], ...
                                InitCond_ODE);

                % Retrieve the number of fruits at harvest time

                    FruitsHarvest = max(0, sol(end,1));

                % Evaluate the economic objective
                    J = N_tratt * P_tratt - F_price * FruitsHarvest;

                % Store the optimum objective function J
                    
                    if J < best_J
                        best_J = J;
                        best_ts = TreatmentTimes;
                    end

                    return
                end

            % Generate combinations

            remainingTreatments = N_tratt - (n_selected + 1);
            
            latestTreatment = T_harvest - T_preharvest ...
                                  - remainingTreatments * T_min;

            for t = current_time:latestTreatment

                current_comb(n_selected + 1) = t;
                recursive_build(n_selected + 1, t + T_min);

            end

            end

        end


    end

end
