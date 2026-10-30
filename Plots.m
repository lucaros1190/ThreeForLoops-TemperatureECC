
% File that plots the results calculated in 'RunmeDDE_HalyomorphaVD.m'

% Created by Luca Rossini on 15 September 2026
% Last update 17 September 2026
% e-mail: l.rossini@reading.ac.uk


%% Compare Fruits and immature stages for 0, 1, ..., N optimal treatments

 % This plot is helpful to understand the effect of different treatments on
 % the population and on the number of fruits.

    % Definition of the colours associated with 0, 1, ..., N treatments

colors = lines(N_tratt_max + 1);

    % States to plot:
    % 1 = Fruits
    % 2 = Eggs
    % 3 = L1
    % 4 = L2
    % 5 = L3
    % 6 = Pupa

StateIndex = 1:6;

StateNames = {'Fruits', 'Eggs', 'Larva 1', 'Larva 2', 'Larva 3', 'Pupae'};

    % Carry out the loop over states

for s = 1:length(StateIndex)

    state = StateIndex(s);

    figure('Position', [100 100 1200 650])
    hold on

    % Manage left axis - Population

    yyaxis left

    % Plot solution 0 treatments

    plot(timeODE, solODE(:,state), 'LineWidth', 2, 'Color', colors(1,:), ...
         'DisplayName', '0 treatments');

    % Plot solution 1, 2, ..., N optimal treatments

    for N_tratt = 1:N_tratt_max

        % Optimal treatment times for this N

        TreatmentTimes_N = BestTs_all{N_tratt};

        % Simulate optimal strategy

        [timeN, solN] = ode45(@(t,y) Fun.ODEfun(t, y, SR, Rates, ...
            FruitAbortion, FruitPredation, TreatmentTimes_N, M_tratt, ...
            TrDecay, T_cover), [0 T_harvest], InitCond_ODE);

        % Colour corresponding to this strategy

        thisColor = colors(N_tratt + 1,:);

        % Plot state

        plot(timeN, solN(:,state), 'LineWidth', 2, 'Color', thisColor, ...
             'DisplayName', sprintf('%d treatment(s)', N_tratt));


        % Treatment times

        for k = 1:length(TreatmentTimes_N)

            xline(TreatmentTimes_N(k), '--', 'LineWidth', 1.2, 'Color', ...
                  thisColor, 'HandleVisibility', 'off');

        end

    end

        % Left axis settings

    if state == 1

        ylabel('Number of fruits', 'FontSize', 16, 'FontWeight', 'bold');

    else

        ylabel('Population density', 'FontSize', 16, 'FontWeight', 'bold');

    end

    % Manage right axis - Temperature

    yyaxis right

    plot(0:length(DailyTemp)-1, DailyTemp, 'k-', 'LineWidth', 3, ...
         'DisplayName', 'Temperature');

    ylabel('Temperature (°C)', 'FontSize', 16, 'FontWeight', 'bold');

    ylim([min(DailyTemp)-2, max(DailyTemp)+2])

    % Figure settings

    xlabel('Time (days)', 'FontSize', 16, 'FontWeight', 'bold');

    title(sprintf('%s - Effect of optimal treatment strategies', ...
          StateNames{s}), 'FontSize', 18, 'FontWeight', 'bold');

    legend('Location', 'best', 'FontSize', 12);

    grid on
    box on

    ax = gca;
    ax.FontSize = 14;
    ax.LineWidth = 1.2;

    xlim([0 T_harvest])

end


%% Compare adult population for 0, 1, ..., N optimal treatments

% This plot is helpful to show the results of the efficacy of the
% treatments

figure('Position', [100 100 1200 650])
hold on

% Define colours for 0, 1, ..., N treatments

colors = lines(N_tratt_max + 1);

% Define left axis - Adult population

yyaxis left

% Plot 0 treatments

Adults0 = solODE(:,7) + solODE(:,8) + solODE(:,9);

plot(timeODE, Adults0, 'LineWidth', 2, 'Color', colors(1,:), 'DisplayName', ...
    '0 treatments');

% Plot 1, 2, ..., N optimal treatments

for N_tratt = 1:N_tratt_max

    TreatmentTimes_N = BestTs_all{N_tratt};

    % Simulate optimal strategy

    [timeN, solN] = ode45(@(t,y) Fun.ODEfun(t, y, SR, Rates, FruitAbortion, ...
        FruitPredation, TreatmentTimes_N, M_tratt, TrDecay, T_cover), ...
        [0 T_harvest], InitCond_ODE);

    % Total adults

    AdultsN = solN(:,7) + solN(:,8) + solN(:,9);

    % Colour corresponding to this strategy

    thisColor = colors(N_tratt + 1,:);

    % Population

    plot(timeN, AdultsN, 'LineWidth', 2, 'Color', thisColor, 'DisplayName', ...
        sprintf('%d treatment(s)', N_tratt));

    % Treatment times

    for k = 1:length(TreatmentTimes_N)

        xline(TreatmentTimes_N(k), '--', ...
            'LineWidth', 1.2, ...
            'Color', thisColor, ...
            'HandleVisibility', 'off');
    end

end

ylabel('Adult population', 'FontSize', 16, 'FontWeight', 'bold')

% Define right axis - Temperature

yyaxis right

plot(0:length(DailyTemp)-1, DailyTemp, 'k-', 'LineWidth', 3, 'DisplayName', ...
     'Temperature');

ylabel('Temperature (°C)', 'FontSize', 16, 'FontWeight', 'bold')

ylim([min(DailyTemp)-2, max(DailyTemp)+2])

% Figure settings

xlabel('Time (days)', 'FontSize', 16, 'FontWeight', 'bold')

title('Effect of optimal treatment strategies on adults', 'FontSize', 18, ...
      'FontWeight', 'bold')

legend('Location', 'best', 'FontSize', 12)

ax = gca;
ax.FontSize = 14;
ax.LineWidth = 1.2;

xlim([0 T_harvest])


%% Plot the dominant instantaneous eigenvalue

figure('Position', [100 100 1200 650])
hold on

    % Compute the dominant eigenvalue

h0 = plot(timeEig, LambdaMax_NoTreat, 'LineWidth', 2.5, 'DisplayName', ...
          'No treatment');

hT = plot(timeEig, LambdaMax_Optimal, 'LineWidth', 2.5, 'DisplayName', ...
          'Optimal treatment');

    % Stability threshold

yline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    % Optimal treatment times

TreatmentColor = hT.Color;

for k = 1:length(BestGlobal_ts)

    xline(BestGlobal_ts(k), '--', 'LineWidth', 1.5, 'Color', ...
          TreatmentColor, 'HandleVisibility', 'off');

end

    % Figure settings

xlabel('Time (days)', 'FontSize', 16, 'FontWeight', 'bold');

ylabel('max Re(\lambda)', 'FontSize', 16, 'FontWeight', 'bold');

title('Instantaneous dominant eigenvalue', 'FontSize', 18, 'FontWeight', ...
      'bold');

legend('Location', 'best', 'FontSize', 12);

grid on
box on

ax = gca;
ax.FontSize = 14;
ax.LineWidth = 1.2;

xlim([0 T_harvest])