
% File containing all the parameters needed by 'RunMe_Drosophila.m'
% Change here the parameters and run directly 'RunMe_Drosophila.m'

% Created by Luca Rossini on 3 November 2023
% Last update 17 September 2026
% E-mail: l.rossini@reading.ac.uk


%% Absorb the temperature array

InputFileTemp = readmatrix('TemperatureInput.xlsx');
DailyTemp = InputFileTemp(:,3);


%% Define the time span

t_span = [0, length(DailyTemp)];

% Define the harvest time

T_harvest = length(DailyTemp);


%% Equation parameters: Plant parameters

    % Natural abortion rate of the fruits (mortality rate fruits)

FruitAbortion = 0.001;

    % Fruit abortion rate due to insects (insect predation on fruits)

FruitPredation = 1;


%% Equation parameters: Treatment parameters

    % Efficacy of the treatment

M_tratt = 0.9;

    % Decay rate of the treatment

TrDecay = 0.8;

    % Coverage time

T_cover = 25;


%% Equation parameters: Economic parameters

    % Price of a single treatment

P_tratt = 50;
    
    % Price of the fruits

F_price = 0.1;


%% Optimisation parameters

    % Maximum number of treatments

N_tratt_max = 3;

    % Minimum interval between two treatments

T_min = T_cover;

    % Minimum interval between last treatment and harvest

T_preharvest = T_cover;


%% Equation parameters: Insect parameters

    % Sex ratio

SR = 0.5;

    % Fertility

alpha = 37000.000000;
gamma = 17.761392;
lambda = 26.036039;
delta = 4.00000;
tau = 22.944989;

FertPar = [alpha gamma lambda delta tau];

    % Mortality Egg - Exponential

a_Mort_Egg = 0.012533;
b_Mort_Egg = -0.482365;
c_Mort_Egg = 2.256795;

MortPar_Egg = [a_Mort_Egg, b_Mort_Egg, c_Mort_Egg];

    % Mortality L1 - Exponential

a_Mort_L1 = 0.0143505;
b_Mort_L1 = -0.5538988;
c_Mort_L1 = 2.584712;

MortPar_L1 = [a_Mort_L1, b_Mort_L1, c_Mort_L1];

    % Mortality L2 - Exponential
    
a_Mort_L2 = 0.017029;
b_Mort_L2 = -0.658377;
c_Mort_L2 = 3.308932;

MortPar_L2 = [a_Mort_L2, b_Mort_L2, c_Mort_L2];

    % Mortality L3 - Exponential

a_Mort_L3 = 0.018840;
b_Mort_L3 = -0.727484;
c_Mort_L3 = 3.652471;

MortPar_L3 = [a_Mort_L3, b_Mort_L3, c_Mort_L3];

    % Mortality Pupa - Constant value

MortPar_P = 0.01;

    % Development - EGGS - Sharpe and De Michele
    
A_Egg = 4.782932;
B_Egg = -113.997964;
C_Egg = 20.442876;
D_Egg = 272.1352;
E_Egg = 6.756630;
F_Egg = -137.529823;

DevPar_Egg = [A_Egg, B_Egg, C_Egg, D_Egg, E_Egg, F_Egg];

    % Development - L1 - Briere

a_L1 = 0.001134;
TL_L1 = 3.148968;
TM_L1 = 32.023358;
m_L1 = 5.686684;

DevPar_L1 = [a_L1, TL_L1, TM_L1, m_L1];

    % Development - L2 - Logan

psi_L2 = 0.043368;
rho_L2 = 0.116276;
TM_L2 = 32.329243;
DeltaT_L2 = 2.217094;

DevPar_L2 = [psi_L2, rho_L2, TM_L2, DeltaT_L2];

    % Development - L3 - Lactin-1

rho_L3 = 0.185309;
TM_L3 = 32.973793;
DeltaT_L3 = 5.374509;

DevPar_L3 = [rho_L3, TM_L3, DeltaT_L3];

    % Development - Pupa - Briere 

a_P = 0.000046;
TL_P = 6.602687;
TM_P = 35.271650;
m_P = 1.00;

DevPar_P = [a_P, TL_P, TM_P, m_P];

    % Survival Ad Males - Lactin-1

a_AM = 0.106797;
TM_AM = 34.999382;
DeltaT_AM = 9.333894;
    
DevPar_AM = [a_AM, TM_AM, DeltaT_AM];

    % Survival Ad Females - Lactin-1

a_AF = 0.097196;
TM_AF = 35.000000;
DeltaT_AF = 10.243252;

DevPar_AF = [a_AF, TM_AF, DeltaT_AF];


%% Initial conditions - ODE

InitCond_ODE = [100000,  % Fruits
                0,       % Egg
                0,       % L1
                0,       % L2
                0,       % L3
                1000,    % Pupa
                0,       % Adult males
                0,       % Non mated females
                0];      % Mated females


% Create the .mat file containing all the parameters

delete("Parameters.mat")
save Parameters.mat t_span InitCond_ODE
