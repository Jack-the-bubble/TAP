function [dCa] = dCa(Ca, T)
    global CAin0 CAin ro cp k E_R h a b V Fin Fc Fc0 Tin Ca0 Tcin T0 F;
    dCa = (Fin * CAin - Fin * Ca - V * k * exp(-E_R / T) * Ca) / V;
end
