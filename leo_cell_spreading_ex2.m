%Leo Harjanto
% Experiment 2 - Graphing all conditions on same plot
clear all

% Model Parameters
n_m = 75;          
F_m = -2;          
v_u = -120;        
n_c = 75;          
k_on = 1;          
k_off = 0.1;       
F_b = -2;          
k_clutch = 5;      

% Figure Setup
figure;
hold on;
xlabel("Relative cell area");
ylabel("Actin Retrograde Flow [nm/s]");
title("Movement vs Cell Spreading");

% Substrate stiffness values
k_substrate_values = [0.1, 1, 10, 100];
cell_areas = 0.5:0.05:4; 

% Colors for different plots
colors = {'r', 'g', 'b', 'k'};

% Loop through substrate stiffness values
for k_idx = 1:length(k_substrate_values)
    k_substrate = k_substrate_values(k_idx);
    x_substrate = 0;
    t = 0;
    events_simulated = 1e5;
    v_filaments_t = zeros(1, events_simulated);
    dts = zeros(1, events_simulated);
    average_vs = zeros(1, length(cell_areas));

    i = 1;
    for cell_area = cell_areas
        new_n_c = round(n_c * cell_area);
        x_clutches = zeros(1, new_n_c);
        clutch_states = zeros(1, new_n_c);

        for ii = 1:events_simulated
            % Calculate filament velocity
            v_filament = v_u * (1 - ((k_substrate * x_substrate) / (n_m * F_m)));

            % Calculate clutch forces
            F_clutch = k_clutch * (x_clutches - x_substrate);
            
            % Simulate next event
            [is_bind, state_idx, dt] = simulate_next_event(clutch_states, F_clutch, k_on, F_b, k_off);

            % Update clutch state
            clutch_states(state_idx) = is_bind;
            if ~is_bind
                x_clutches(state_idx) = x_substrate;
            end

            engaged_clutch_idx = find(clutch_states == 1);
            disengaged_clutch_idx = find(clutch_states == 0);

            % Update substrate position
            x_substrate = (k_clutch * sum(x_clutches(engaged_clutch_idx))) / ...
                          (k_substrate + length(engaged_clutch_idx) * k_clutch);

            % Update clutch displacements
            x_clutches(engaged_clutch_idx) = x_clutches(engaged_clutch_idx) + v_filament * dt;
            x_clutches(disengaged_clutch_idx) = x_substrate;
            
            t = t + dt;
            v_filaments_t(ii) = v_filament;
            dts(ii) = dt;
        end
        average_vs(i) = sum(v_filaments_t .* dts) / sum(dts);
        i = i + 1;
    end
    
    % Plot results
    plot(cell_areas, -average_vs, '-', 'Color', colors{k_idx});
end

legend('k_{substrate}=0.1', 'k_{substrate}=1', 'k_{substrate}=10', 'k_{substrate}=100');
hold off;

% Function Definition Moved to the End
function [is_bind, state_idx, dt] = simulate_next_event(clutch_states, F_clutch, k_on, F_b, k_off)
    unengaged_clutch_indices = find(clutch_states == 0);
    engaged_clutch_indices = find(clutch_states == 1);

    if isempty(unengaged_clutch_indices)
        t_bind = inf;
    else
        t_bind = -log(rand(1, length(unengaged_clutch_indices))) / k_on;
    end

    if isempty(engaged_clutch_indices)
        t_unbind = inf;
    else
        t_unbind = -log(rand(1, length(engaged_clutch_indices))) ./ ...
                   (k_off * exp(F_clutch(engaged_clutch_indices) / F_b));
    end

    [min_t_unbind, unbind_idx] = min(t_unbind);
    [min_t_bind, bind_idx] = min(t_bind);

    if min_t_bind < min_t_unbind
        is_bind = 1;
        state_idx = unengaged_clutch_indices(bind_idx);
        dt = min_t_bind;
    else
        is_bind = 0;
        state_idx = engaged_clutch_indices(unbind_idx);
        dt = min_t_unbind;
    end
end
