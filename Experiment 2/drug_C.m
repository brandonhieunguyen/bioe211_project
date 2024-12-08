function C = drug_C(Ka,Kd,V,Dose,C0,t)
    C = ((Ka*Dose)/(V*(Ka - (Kd/V))))*(exp(-(Kd/V)*t)-exp(-Ka*t))+C0*exp(-(Kd/V)*t);
end