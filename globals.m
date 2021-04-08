global g m K kmol cal min ro cp k E_R h a b ro cp k E_R h a b V Fin CAin Fc Tin Tcin Ca T F V0 Fin0 CAin0 Fc0 Tin0 Tcin0 Ca0 F0 T0;

% units
g = 1;
m = 1;
K = 1;
kmol = 1;
cal = 1;
min = 1;

% const
ro = 1e6 * g / m^3;
cp = 1 * cal / (g * K);
k = 1e10 * 1 / min;
% k0 =k;
E_R = 8330.1 * 1 / K;
h = 130 * 1e6 * cal / kmol;
a = 0.516 * 1e6 * cal / (K * m^3);
b = 0.5;

% operating point
V = 1 * m^3;
Fin = 1 * m^3 / min;
CAin = 2 * kmol / m^3;
Fc = 15 * m^3 / min;
F = Fin;
Tin = 343 * K;
Tcin = 310 * K;
Ca = 1.7895 * kmol / m^3;
T = 331.0083 * K;

% operating point under zero conditions
V0 = V;
Fin0 = Fin;
CAin0 = CAin;
Fc0 = Fc;
F0 = Fin0;
Tin0 = Tin;
Tcin0 = Tcin;
Ca0 = Ca;
T0 = T;
