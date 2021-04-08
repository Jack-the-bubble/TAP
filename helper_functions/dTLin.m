function [dT] = dTLin(Ca, T)
    global CAin CAin0 ro cp k E_R h a b V Fin Fc Fc0 Tin Ca0 Tcin T0 F;
    dT = (a * (Fc - Fc0) * (T0 - Tcin) * ((exp(log(Fc0) * (b + 1)) * ((a * b * exp(b * log(Fc0))) / (2 * Fc0 * cp * ro) + 1)) / (Fc0 + (a * exp(b * log(Fc0))) / (2 * cp * ro))^2 - (exp(log(Fc0) * (b + 1)) * (b + 1)) / (Fc0 * (Fc0 + (a * exp(b * log(Fc0))) / (2 * cp * ro)))) - (T - T0) * ((a * exp(log(Fc0) * (b + 1))) / (Fc0 + (a * exp(b * log(Fc0))) / (2 * cp * ro)) + F * cp * ro - (Ca0 * E_R * V * h * k * exp(-E_R / T0)) / T0^2) - (a * exp(log(Fc0) * (b + 1)) * (T0 - Tcin)) / (Fc0 + (a * exp(b * log(Fc0))) / (2 * cp * ro)) - F * T0 * cp * ro + Fin * Tin * cp * ro + Ca0 * V * h * k * exp(-E_R / T0) + V * h * k * exp(-E_R / T0) * (Ca - Ca0)) / (V * cp * ro);
end
