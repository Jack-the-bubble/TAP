close all;
clear all;

global g m K mol kmol cal min ro cp k E_R h a b ro cp k E_R h a b V Fin CAin Fc Tin Tcin Ca T;
globals
st = 0.01;

CAin_vect = [1.6,1.7, 1.8, 1.9, 1.95, 2.0,2.2] .* kmol;
Fc_vect = [13,14, 15, 15.5, 16, 17,19];

Ca = 1.7845;
T = 331.0083;

CAin = 2; Fc = 15;


figure(1)
PlotModel('CAin', 'Ca', CAin_vect, @dCaLin, @dTLin, st, 'Lin zaleznosc Ca od czasu - skok CAin', '-');
PlotModel('CAin', 'Ca', CAin_vect, @dCa, @dT, st, 'NLin zaleznosc Ca od czasu - skok CAin', '-.');

plotLegend('CAin', CAin_vect)
hold off

CAin = 2; Fc = 15;

figure(2)
PlotModel('Fc', 'Ca', Fc_vect, @dCaLin, @dTLin, st, 'Lin zaleznosc Ca od czasu - skok Fc', '-');
PlotModel('Fc', 'Ca', Fc_vect, @dCa, @dT, st, 'NLin zaleznosc Ca od czasu - skok Fc', '-.');
plotLegend('Fc', Fc_vect)
hold off

CAin = 2; Fc = 15;

figure(3)
PlotModel('CAin', 'T', CAin_vect, @dCaLin, @dTLin, st, 'Lin zaleznosc T od czasu - skok CAin', '-');
PlotModel('CAin', 'T', CAin_vect, @dCa, @dT, st, 'NLin zaleznosc T od czasu - skok CAin', '-.');
plotLegend('CAin', CAin_vect)
hold off

CAin = 2; Fc = 15;


figure(4)
PlotModel('Fc', 'T', Fc_vect, @dCaLin, @dTLin, st, 'Lin zaleznosc T od czasu - skok Fc', '-');
PlotModel('Fc', 'T', Fc_vect, @dCa, @dT, st, 'NLin zaleznosc T od czasu - skok Fc', '-.');
plotLegend('Fc', Fc_vect)
hold off
CAin = 2; Fc = 15;

% compare ss with continous model
figure(5)
PlotModel('CAin', 'Ca', CAin_vect, @dCaLin, @dTLin, st, 'Lin zaleznosc Ca od czasu - skok CAin', '-');
PlotModelCont('CAin', 'Ca', CAin_vect, @dCa, @dT, st, 'NLin zaleznosc Ca od czasu - skok CAin', '-.');

plotLegend('CAin', CAin_vect)
hold off

CAin = 2; Fc = 15;

figure(6)
PlotModel('Fc', 'Ca', Fc_vect, @dCaLin, @dTLin, st, 'Lin zaleznosc Ca od czasu - skok Fc', '-');
PlotModel('Fc', 'Ca', Fc_vect, @dCa, @dT, st, 'NLin zaleznosc Ca od czasu - skok Fc', '-.');
plotLegend('Fc', Fc_vect)
hold off

CAin = 2; Fc = 15;

figure(7)
PlotModel('CAin', 'T', CAin_vect, @dCaLin, @dTLin, st, 'Lin zaleznosc T od czasu - skok CAin', '-');
PlotModel('CAin', 'T', CAin_vect, @dCa, @dT, st, 'NLin zaleznosc T od czasu - skok CAin', '-.');
plotLegend('CAin', CAin_vect)
hold off

CAin = 2; Fc = 15;

figure(8)
PlotModel('Fc', 'T', Fc_vect, @dCaLin, @dTLin, st, 'Lin zaleznosc T od czasu - skok Fc', '-');
PlotModel('Fc', 'T', Fc_vect, @dCa, @dT, st, 'NLin zaleznosc T od czasu - skok Fc', '-.');
plotLegend('Fc', Fc_vect)
hold off
CAin = 2; Fc = 15;


