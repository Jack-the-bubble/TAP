close all;
clear all;

global g m K mol kmol cal min ro cp k E_R h a b ro cp k E_R h a b V Fin CAin Fc Tin Tcin Ca T;
globals
st = 0.01;

CAin_vect = [1.6,1.7, 1.8, 1.9, 1.95, 2.0,2.2];
Fc_vect = [13,14, 15, 15.5, 16, 17,19];

Ca = 1.7845;
T = 331.0083;

CAin = 2; Fc = 15;

Ts = 10;
figure(1)
PlotModel('CAin', 'Ca', CAin_vect, @dCaLin, @dTLin, st, 'M Lin zaleznosc Ca od czasu - skok CAin', '-');
PlotModelDiscrete('CAin', 'Ca', Ts, CAin_vect, @dCaLin, @dTLin, st, 'M Dis zaleznosc Ca od czasu - skok CAin', '-');
plotLegend('CAin', CAin_vect)
hold off

figure(2)
PlotModel('CAin', 'T', CAin_vect, @dCaLin, @dTLin, st, 'M Lin zaleznosc T od czasu - skok CAin', '-');
PlotModelDiscrete('CAin', 'T', Ts, CAin_vect, @dCaLin, @dTLin, st, 'M Dis zaleznosc T od czasu - skok CAin', '-');
plotLegend('CAin', CAin_vect)
hold off
CAin = 2;

figure(3)
PlotModel('Fc', 'Ca', Fc_vect, @dCaLin, @dTLin, st, 'M Lin zaleznosc Ca od czasu - skok Fc', '-');
PlotModelDiscrete('Fc', 'Ca', Ts, Fc_vect, @dCaLin, @dTLin, st, 'Zaleznosc Ca od czasu - skok Fc', '-');
plotLegend('Fc', Fc_vect)
hold off

figure(4)
PlotModel('Fc', 'T', Fc_vect, @dCaLin, @dTLin, st, 'M Lin zaleznosc T od czasu - skok Fc', '-');
PlotModelDiscrete('Fc', 'T', Ts, Fc_vect, @dCaLin, @dTLin, st, 'M Dis zaleznosc T od czasu - skok Fc', '-');
plotLegend('Fc', Fc_vect)
hold off