function dy = sym_nl(t, y, CAin, Fc)
% y(1) - Ca
% y(2) - T
% dy(1) - dCa/dt
% dy(2) - dT/dt
Ca = y(1);
T = y(2);

dy = zeros(2, 1);
global g m K mol kmol cal min ro ro_c cp c_pc k E_R h a b V F Fin Tin Tcin;

dy(1) = (Fin * CAin - Fin * Ca - V * k * exp(-E_R / T) * Ca) / V;
dy(2) = (Fin * ro * cp * Tin - Fin * ro * cp * T + V * h * k * exp(-E_R / T) * Ca - (a * Fc^(b + 1) / (Fc + (a * Fc^b / (2 * ro * cp)))) * (T - Tcin)) / (V * ro * cp);
% dy(1) = (Fin*CAin)/V - (F * y(1))/V - k*exp(E_R/y(2)) *y(1);
% dy(2) = (Fin*Tin)/V - (F * y(2))/V + h*k*exp(E_R/y(2))*y(1)/(ro*cp)-...
%     (y(2)-Tcin)*((a*Fc^(b+1))/(Fc+(a*Fc^b)/(2*ro_c*c_pc)))/(V*ro*cp);
end
