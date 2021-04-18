function [dT] = dT(Ca, T)
    global CAin0 CAin ro cp k E_R h a b V Fin Fc Fc0 Tin Ca0 Tcin T0 F;
    dT = (Fin * ro * cp * Tin - Fin * ro * cp * T + V * h * k * exp(-E_R / T) * Ca - (a * Fc^(b + 1) / (Fc + (a * Fc^b / (2 * ro * cp)))) * (T - Tcin)) / (V * ro * cp);
end
