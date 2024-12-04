function [is_bind, state_idx, dt] = simulate_next_event(clutch_states, F_clutch, k_on, F_b, k_off)
    % Gillepsie algorithm implementation
    unengaged_clutch_indices = find(clutch_states==0);
    engaged_clutch_indices = find(clutch_states==1);

    if isempty(unengaged_clutch_indices)
        t_bind=inf;
    else
        t_bind=-log(rand(1, length(unengaged_clutch_indices)))/k_on;
    end

    if isempty(engaged_clutch_indices)
        t_unbind=inf;
    else
        t_unbind=-log(rand(1, length(engaged_clutch_indices)))./(k_off * exp(F_clutch(engaged_clutch_indices) / F_b));
    end

    [min_t_unbind, unbind_idx] = min(t_unbind);
    [min_t_bind, bind_idx] = min(t_bind);
    
    if min_t_bind < min_t_unbind
        % Binding event happens first
        is_bind = 1;
        state_idx = unengaged_clutch_indices(bind_idx);
        dt = min_t_bind;
    else
        % Unbinding event happens first
        is_bind = 0;
        state_idx = engaged_clutch_indices(unbind_idx);
        dt = min_t_unbind;
    end
end