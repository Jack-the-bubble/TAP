function [dCa] = dCaLin(Ca, T)
    global CAin0 CAin ro cp k E_R h a b V Fin Fc Fc0 Tin Ca0 Tcin T0 F;
    dCa = -((Ca - Ca0) * (F + V * k * exp(-E_R / T0)) + Ca0 * F - CAin0 * Fin - Fin * (CAin - CAin0) + Ca0 * V * k * exp(-E_R / T0) + (Ca0 * E_R * V * k * exp(-E_R / T0) * (T - T0)) / T0^2) / V;
end
