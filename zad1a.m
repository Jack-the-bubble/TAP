clear all;
close all;
syms A B C D s Vt(t) Tint(t) Tcint(t) CAint(t) Fct(t) Ft(t) Fint(t) Cat(t) Tt(t) Ca g m K mol kmol cal roc cpc ro cp k k0 E_R h a b ro cp k E_R h a b V F Fin CAin Fc Tin Tcin T V0 Fin0 CAin0 Fc0 Tin0 Tcin0 Ca0 T0 Fin0 F0;

syms Vs Tins Tcins CAins Fcs Fs Cas Ts Fins

res = fsolve(@system, [1.79, 331]);

eq1 = (Fin * CAin - F * Ca - V * k0 * exp(-E_R / T) * Ca);
eq2 = ((Fin * ro * cp * Tin) - (F * ro * cp * T) + (V * h * k0 * exp(-E_R / T) * Ca) - ((a * (Fc^(b + 1)) * (T - Tcin)) / (Fc + ((a * Fc^b) / (2 * roc * cpc))))); 

lin1 = taylor(eq1, [CAin Fc Ca T], [CAin0 Fc0 Ca0 T0], 'Order', 2);
lin2 = taylor(eq2, [CAin Fc Ca T], [CAin0 Fc0 Ca0 T0], 'Order', 2);

lin1 = lin1 / V;
lin2 = lin2 / (V * ro * cp);
latex(lin1);
latex(lin2);

% linearization
eq1 = diff(Ca, t) == lin1;
eq2 = diff(Tt, t) == lin2;

A = sym('A%d%d', [2 2]);
B = sym('B%d%d', [2 2]);
C = [1 0; 0 1];
D = zeros(2, 2);

A = aV(A, 1, 1, lin1, Ca);
A = aV(A, 1, 2, lin1, T);
A = aV(A, 2, 1, lin2, Ca);
A = aV(A, 2, 2, lin2, T);

B = aV(B, 1, 1, lin1, CAin);
B = aV(B, 1, 2, lin1, Fc);
B = aV(B, 2, 1, lin2, CAin);
B = aV(B, 2, 2, lin2, Fc);

ro = 1e6;
roc = ro;
cp = 1;
cpc = cp;
k = 1e10;
k0 = k;
E_R = 8330.1;
h = 130 * 1e6;
a = 0.516 * 1e6;
b = 0.5;
V = 1;
Fin = 1;
F = Fin;
CAin = 2;
Fc = 15;
Tin = 343;
Tcin = 310;
Ca = 1.7895;
T = 331.0083;
T0 = T;
Fc0 = Fc;
Ca0 = Ca;
CAin0 = CAin;

A = [-(F + V * k0 * exp(-E_R / T0)) / V, -(Ca0 * E_R * k0 * exp(-E_R / T0)) / T0^2; (h * k0 * exp(-E_R / T0)) / (cp * ro), -((a * exp(log(Fc0) * (b + 1))) / (Fc0 + (a * exp(b * log(Fc0))) / (2 * cpc * roc)) + F * cp * ro - (Ca0 * E_R * V * h * k0 * exp(-E_R / T0)) / T0^2) / (V * cp * ro)];
B = [Fin / V, 0; 0, (a * (T0 - Tcin) * ((exp(log(Fc0) * (b + 1)) * ((a * b * exp(b * log(Fc0))) / (2 * Fc0 * cpc * roc) + 1)) / (Fc0 + (a * exp(b * log(Fc0))) / (2 * cpc * roc))^2 - (exp(log(Fc0) * (b + 1)) * (b + 1)) / (Fc0 * (Fc0 + (a * exp(b * log(Fc0))) / (2 * cpc * roc))))) / (V * cp * ro)];

SS = ss(A, B, C, D);
output = tf(SS);
disp('Model w postaci transmitancji')
output

disp('Model w postaci zmiennych stanu')
SS

opt = stepDataOptions('StepAmplitude', 0)


step(output)
step_output = step(output);

function res = system(CaT)
    globals
    Ca_sym = CaT(1);
    T_sym = CaT(2);

    res(1) = (Fin * CAin) / V - F * Ca_sym / V - k * exp(-E_R / T_sym) * Ca_sym;
    res(2) = Fin * Tin / V - F * T_sym / V + h * k * exp(-E_R / T_sym) * Ca_sym / (ro * cp) - (a * (Fc)^(b + 1)) / (Fc + (a * (Fc)^(b)) / (2 * ro * cp)) * (T_sym - Tcin) / (V * ro * cp);
end
