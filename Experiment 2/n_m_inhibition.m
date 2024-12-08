function [nm_active, nm_states, treatment_day, treatment_time, C, C0] = n_m_inhibition(nm,nm_states,nm_min,Dose,treatment_day,treatment_time,dt,t,C0_)
    %Pharmacokinetics specs
    V = 1;                          % [mL] Volume units
    D_t = 10;                       % [Days] Degradation time
    a_t = 2/24;                     % [Days] Absorbtion time
    Kd = 1/D_t;                     % [1/Days]Degradation rate
    Ka = 1/a_t;                     % [1/Days] Absorbtion rate
    Kb_max =  0.97;                 % [1/Days] Max binding rate of drug (97% in Mavactem)
    Kb_min = 0;                     % [1/Days] Minimum binding rate of drug
    Kr_max = 1;                     % [1/Days] Max unbinding rate of drug/Binding rate of ATP
    Kr_min = 0.03;                  % [1/Days] Minimum unbinding rate of drug/Binding rate of ATP
    target_ss_concentration = 500;  % [ng/mL] Terapeutic target steady state concentration
                                    % (ng/mL)
    C0 = C0_;
    C = drug_C(Ka,Kd,V,Dose,C0,t-treatment_time);
    if floor(t) >= treatment_day 
        C0 = C;
        treatment_day = treatment_day + 1;
        treatment_time = t;
    end
    
    if C <= target_ss_concentration && t > 0
        Kb = Kb_max*C/target_ss_concentration;
        Kr = Kr_max - Kb;
    else
        Kb = Kb_max;
        Kr = Kb_min;
    end

    for j = 1:nm
        if Kb > rand()
            nm_state(j) = 0;
        elseif Kr > rand()
            nm_state(j) = 1;
        end
    end
    if sum(nm_state) < nm_min
        nm_active = nm_min;
    else
        nm_active = sum(nm_state);
    end

end