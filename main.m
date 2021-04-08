close all;
clear all;

global g m K mol kmol cal min ro cp k E_R h a b ro cp k E_R h a b V Fin CAin Fc Tin Tcin Ca T;
globals
step = 0.01;

CAin_vect = [1.7, 1.8, 1.9, 1.95, 2.0] .* kmol;
Fc_vect = [14, 15, 15.5, 16, 17];

Ca = 1.7895 * kmol / m^3;
T = 331.0083 * K;

CAin = 2; Fc = 15;

%% task b)

figure(2)
PlotModel('CAin', 'Ca', CAin_vect, @dCaLin, @dTLin, step, 'Lin zaleznosc Ca od czasu - skok CAin', '-');
PlotModel('CAin', 'Ca', CAin_vect, @dCa, @dT, step, 'NLin zaleznosc Ca od czasu - skok CAin', '-.');

plotLegend('Ca', CAin_vect)
hold off

CAin = 2; Fc = 15;

figure(3)
PlotModel('Fc', 'Ca', Fc_vect, @dCaLin, @dTLin, step, 'Lin zaleznosc Ca od czasu - skok Fc', '-');
PlotModel('Fc', 'Ca', Fc_vect, @dCa, @dT, step, 'NLin zaleznosc Ca od czasu - skok Fc', '-.');
plotLegend('Fc', Fc_vect)
hold off

CAin = 2; Fc = 15;

figure(4)
PlotModel('CAin', 'T', CAin_vect, @dCaLin, @dTLin, step, 'Lin zaleznosc T od czasu - skok CAin', '-');
PlotModel('CAin', 'T', CAin_vect, @dCa, @dT, step, 'NLin zaleznosc T od czasu - skok CAin', '-.');
plotLegend('Ca', CAin_vect)
hold off

CAin = 2; Fc = 15;

figure(5)
PlotModel('Fc', 'T', Fc_vect, @dCaLin, @dTLin, step, 'Lin zaleznosc T od czasu - skok Fc', '-');
PlotModel('Fc', 'T', Fc_vect, @dCa, @dT, step, 'NLin zaleznosc T od czasu - skok Fc', '-.');
plotLegend('Fc', Fc_vect)
hold off
CAin = 2; Fc = 15;

%% linear model
Ts = 10;
figure(6)
PlotModel('CAin', 'Ca', CAin_vect, @dCaLin, @dTLin, step, 'M Lin zaleznosc Ca od czasu - skok CAin', '-');
PlotModelDiscrete('CAin', 'Ca', Ts, CAin_vect, @dCaLin, @dTLin, step, 'M Dis zaleznosc Ca od czasu - skok CAin', '-');
plotLegend('Ca', CAin_vect)
hold off

figure(7)
PlotModel('CAin', 'T', CAin_vect, @dCaLin, @dTLin, step, 'M Lin zaleznosc T od czasu - skok CAin', '-');
PlotModelDiscrete('CAin', 'T', Ts, CAin_vect, @dCaLin, @dTLin, step, 'M Dis zaleznosc T od czasu - skok CAin', '-');
plotLegend('Fc', Fc_vect)
hold off
CAin = 2;

figure(8)
PlotModel('Fc', 'Ca', Fc_vect, @dCaLin, @dTLin, step, 'M Lin zaleznosc Ca od czasu - skok Fc', '-');
PlotModelDiscrete('Fc', 'Ca', Ts, Fc_vect, @dCaLin, @dTLin, step, 'Zaleznosc Ca od czasu - skok Fc', '-');
plotLegend('Ca', CAin_vect)
hold off

figure(9)
PlotModel('Fc', 'T', Fc_vect, @dCaLin, @dTLin, step, 'M Lin zaleznosc T od czasu - skok Fc', '-');
PlotModelDiscrete('Fc', 'T', Ts, Fc_vect, @dCaLin, @dTLin, step, 'M Dis zaleznosc T od czasu - skok Fc', '-');
plotLegend('Fc', Fc_vect)
hold off
Fc = 15;

%% discrete model

Ts = 10;
figure(9)
PlotModelDiscrete('CAin', 'Ca', Ts, CAin_vect, @dCaLin, @dTLin, step, 'M Dis zaleznosc Ca od czasu - skok CAin', '-');
hold off

figure(10)
PlotModelDiscrete('CAin', 'T', Ts, CAin_vect, @dCaLin, @dTLin, step, 'M Dis zaleznosc T od czasu - skok CAin', '-');
hold off

figure(11)
PlotModelDiscrete('Fc', 'Ca', Ts, Fc_vect, @dCaLin, @dTLin, step, 'M Dis zaleznosc Ca od czasu - skok Fc', '-');
hold off

figure(12)
PlotModelDiscrete('Fc', 'T', Ts, Fc_vect, @dCaLin, @dTLin, step, 'M Dis zaleznosc T od czasu - skok Fc', '-');
hold off

%% state equations discrete model
Ts=5;
[dis_ss, dis_tf] = get_discrete(Ts);
[Az,Bz,Cz,Dz]=ssdata(dis_ss);


PlotStateEqDiscrete('CAin', 'Ca', CAin_vect, Fc_vect, Az, Bz, Cz, Dz, 'M Dis zaleznosc T od czasu - skok Fc', '-')


