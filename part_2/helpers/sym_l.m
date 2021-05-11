function dy = sym_nl(t, y, CAin, Fc, Ca0, CAin0, T0, Fc0)
% y(1) - Ca
% y(2) - T
% dy(1) - dCa/dt
% dy(2) - dT/dt
Ca = y(1);
T = y(2);

dy = zeros(2, 1);
global g m K mol kmol cal min ro ro_c cp c_pc k E_R h a b V F Fin Tin Tcin;

dy(1) = -((Ca - Ca0) * (F + V * k * exp(-E_R / T0)) + Ca0 * F ...
        - CAin0 * Fin - Fin * (CAin - CAin0) + Ca0 * V * k * exp(-E_R / T0) ...
        + (Ca0 * E_R * V * k * exp(-E_R / T0) * (T - T0)) / T0^2) / V;
dy(2) = (a * (Fc - Fc0) * (T0 - Tcin) ...
        * ((exp(log(Fc0) * (b + 1)) * ((a * b * exp(b * log(Fc0))) ...
        / (2 * Fc0 * cp * ro) + 1)) ...
        / (Fc0 + (a * exp(b * log(Fc0))) / (2 * cp * ro))^2 ...
        - (exp(log(Fc0) * (b + 1)) * (b + 1)) ...
        / (Fc0 * (Fc0 + (a * exp(b * log(Fc0))) / (2 * cp * ro)))) ...
        - (T - T0) * ((a * exp(log(Fc0) * (b + 1))) ...
        / (Fc0 + (a * exp(b * log(Fc0))) / (2 * cp * ro)) ...
        + F * cp * ro - (Ca0 * E_R * V * h * k * exp(-E_R / T0)) / T0^2) ...
        - (a * exp(log(Fc0) * (b + 1)) * (T0 - Tcin)) ...
        / (Fc0 + (a * exp(b * log(Fc0))) / (2 * cp * ro)) ...
        - F * T0 * cp * ro + Fin * Tin * cp * ro ...
        + Ca0 * V * h * k * exp(-E_R / T0) ...
        + V * h * k * exp(-E_R / T0) * (Ca - Ca0)) / (V * cp * ro);
end