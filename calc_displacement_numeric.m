% function [clutch_displacement, x_substrate_new, v_substrate_new] = calc_displacement(eta, k_substrate, k_clutch, n_eng, v_u, F_stall, t_, dt, x_substrate_i, v_substrate_i)
%     a = eta;
%     b = (k_substrate + k_clutch * n_eng + (k_clutch * n_eng * v_u)/(F_stall));
%     c = (k_clutch * n_eng * v_u * k_substrate) / F_stall;
%     d = k_clutch * n_eng * v_u;
% 
%     syms x(t);
%     Dx(t) = diff(x, t);
%     ode = a * diff(x, t, 2) + b * diff(x, t) + c * x == d;
%     cond1 = x(t_) == x_substrate_i;
%     cond2 = Dx(t_) == v_substrate_i;
%     xSol(t) = dsolve(ode, [cond1 cond2]);
%     DxSol(t) = diff(xSol, t);
%     v_f(t) = v_u * (1 - (eta * DxSol(t) + k_substrate * xSol(t))/(F_stall));
%     clutch_displacement = double(int(v_f(t), t_, t_+dt));
%     x_substrate_new = double(xSol(t_+dt));
%     v_substrate_new = double(DxSol(t_+dt));
% end

function [clutch_displacement, x_substrate_new, v_substrate_new] = calc_displacement_numeric(eta, k_substrate, k_clutch, n_eng, v_u, F_stall, t_, dt, x_substrate_i, v_substrate_i)
    % Parameters for ODE
    a = eta;
    b = (k_substrate + k_clutch * n_eng + (k_clutch * n_eng * v_u)/(F_stall));
    c = (k_clutch * n_eng * v_u * k_substrate) / F_stall;
    d = k_clutch * n_eng * v_u;

    % Define the system of first-order ODEs
    ode_func = @(t, y) [y(2); (-b * y(2) - c * y(1) + d) / a];
    
    % Initial conditions: [x(0), Dx(0)]
    y0 = [x_substrate_i, v_substrate_i];
    
    % Time span for solving
    t_span = [t_, t_ + dt];
    
    % Solve using ode45
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    [t_sol, y_sol] = ode45(ode_func, t_span, y0, options);
    
    % Extract displacement and velocity
    x_substrate_new = y_sol(end, 1);
    v_substrate_new = y_sol(end, 2);
    
    % Calculate velocity function v_f numerically
    v_f = v_u * (1 - (eta * y_sol(:, 2) + k_substrate * y_sol(:, 1)) / F_stall);
    
    % Numerical integration of v_f using trapezoidal rule
    clutch_displacement = trapz(t_sol, v_f);
end