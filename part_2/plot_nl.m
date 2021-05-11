clear all;
close all;

global CAin Fc Ca T Tcin0 Tin;

globals

Tsim = 1e1;



figure(1)
% out1 - Ca
% out2 - T

hold on

Fc_vec = [Fc Fc]
Cain_vec = [CAin CAin+1];

% skok CAin (in1)
for i=1:length(Cain_vec)

[t_nl, y_nl] = ode45(@(t, y) sym_nl(t, y, Cain_vec(i), Fc_vec(i)), [0 Tsim], [Ca T]);
[t_l, y_l] = ode45(@(t, y) sym_l(t, y, Cain_vec(i), Fc_vec(i), Ca, CAin, T, Fc), ...
    [0 Tsim], [Ca T]);

subplot(2, 2, 1); plot(t_l, y_l(:, 1)); %in1 out1
subplot(2, 2, 3); plot(t_l, y_l(:, 2)); %in1 out2
end



% symulacja skoku Fc (in2)
Fc_vec = [Fc Fc+1]
Cain_vec = [CAin CAin];
for i=1:length(Cain_vec)

[t_nl, y_nl] = ode45(@(t, y) sym_nl(t, y, Cain_vec(i), Fc_vec(i)), [0 Tsim], [Ca T]);
[t_l, y_l] = ode45(@(t, y) sym_l(t, y, Cain_vec(i), Fc_vec(i), Ca, CAin, T, Fc), ...
    [0 Tsim], [Ca T]);

subplot(2, 2, 2); plot(t_l, y_l(:, 1)); % in2 out1
subplot(2, 2, 4); plot(t_l, y_l(:, 2)); % in2 out2
end



% state system
T_L = T0;
C_a_L = Ca0;
C_Ain_L = CAin0;
F_c_L = Fc0;
T_cin = Tcin0;
T_in = Tin;

% moje parametry z pliku sym_l.m, ale raczej te poniżej są lepsze
% a11 = -(F/V+k*exp(-E_R/T0));
% a12 = -(Ca0*E_R*k*exp(-E_R/T0)/T0^2);
% a21 = h*k*exp(-E_R/T0)/V;
% a22 = -(F/V+ (a*exp((b+1)*log(Fc0))/(Fc0+(a*exp(b*log(Fc0)))/(2*cp*ro))) -(Ca0*E_R*V*h*k*exp(-E_R/T0))/(T0^2))/(V*cp*ro);
% b1 = Fin/V;
% b2 = a*(T0-Tcin)* ((exp(log(Fc0) * (b + 1)) * ((a * b * exp(b * log(Fc0))) ...
%         / (2 * Fc0 * cp * ro) + 1)) ...
%         / (Fc0 + (a * exp(b * log(Fc0))) / (2 * cp * ro))^2 ...
%         - (exp(log(Fc0) * (b + 1)) * (b + 1)) ...
%         / (Fc0 * (Fc0 + (a * exp(b * log(Fc0))) / (2 * cp * ro))))/(V*cp*ro)


% old, probably good
% a11=-(F/V+k*exp(-(E_R/T_L))); %ok
% a12=-(C_a_L*k*E_R*exp((-E_R/T_L))/T_L^2); %ok
% a21=((h*k)/(ro*cp*exp(E_R/T_L))); %ok
% a22=-F/V + (h*k*E_R*C_a_L*exp(-E_R/T_L))/(ro*cp*T_L^2)-[a*(F_c_L^(b+1)/((V*ro*cp)*(F_c_L+(a*F_c_L^b/(2*ro*cp)))))]; %do sprawdzenia %ok 
a11=-(F/V+k*exp(-(E_R/T_L))); %ok
a12=-(C_a_L*k*E_R*exp((-E_R/T_L))/T_L^2); %ok
a21=((h*k)/(ro*cp*exp(E_R/T_L))); %ok
a22=-F/V + (h*k*E_R*C_a_L*exp(-E_R/T_L))/(ro*cp*T_L^2)-[a*(F_c_L^(b+1)/((V*ro*cp)*(F_c_L+(a*F_c_L^b/(2*ro*cp)))))]; %do sprawdzenia %ok 

b1=F/V; %ok
b2=-((a*(T_L-T_cin)/(V*ro*cp))*(((b+1)*F_c_L^b)/((a*F_c_L^b)/(2*ro*cp)+F_c_L)-(F_c_L^(b+1)*((a*b*F_c_L^(b-1))/(2*ro*cp)+1))/((a*F_c_L^b)/(2*ro*cp)+F_c_L)^2)); %ok

c1=1;
c2=1;
A=[a11 a12; a21 a22];
B=[b1 0; 0 b2];
C=[1 0; 0 1];
D=[0 0; 0 0];




%sta³e procesu - wartoœci potrzebne w modelu zlinearyzowanym w postaci
%równañ stanu
% w1=F*(C_Ain_L-C_a_L)/V-(C_a_L*k*E_R*exp((-E_R/T_L)))+a11*C_a_L+a12*T_L;
% w2=a21*C_a_L+a22*T_L+(F*(T_in-T_L)/V)+[a*(F_c_L^(b+1)/((V*ro*cp)*(F_c_L+(a*F_c_L^b/(2*ro*cp)))))];

%wyznaczanie transmitancji modelu zlinearyzowanego
sys_c=tf(ss(A,B,C,D));
sys_d=c2d(sys_c,0.001,'zoh');
sys11=sys_c(1,1);
sys12=sys_c(1,2);
sys21=sys_c(2,1);
sys22=sys_c(2,2);
sys11d=sys_d(1,1);
sys12d=sys_d(1,2);
sys21d=sys_d(2,1);
sys22d=sys_d(2,2);

opts = stepDataOptions('StepAmplitude', 1)
% y = step(sys12, opts);
% plot(y)

figure(2)
step(sys_c, opts)


