T_sim = 30;

Ca0 = 1.7895;
T0 = 331.0083;
Cain0 = 2;
Fc0 = 15;

tsim = transpose(1:T_sim);
T_zad = timeseries(zeros(T_sim, 1), tsim);
C_a_zad = timeseries(zeros(T_sim, 1), tsim);

T_zad.Data = T_zad.Data + T0;
C_a_zad.Data = C_a_zad.Data + Ca0;

T_zad.Data(10:30) = T_zad.Data(10:30) + 100; 

%% PID bez odsprzegania

sim('PID_bez_odsprz', T_sim)

figure(1);
subplot(2, 1, 1);
hold on;
title({'PID bez odsprzęgania'});
stairs(T_zad.Time, T_zad.Data, 'Color', 'blue');
plot(T_out.Time, T_out.Data, 'Color', 'red');
xlabel('czas symulacji') ;
ylabel('T');
legend('T zad','T');


subplot(2, 1, 2);
hold on;
stairs(C_a_zad.Time, C_a_zad.Data, 'Color', 'blue');
plot(C_a_out.Time, C_a_out.Data, 'Color', 'red');
xlabel('czas symulacji') ;
ylabel('Ca');
legend('Ca zad','Ca');


%% PID z odsprzeganiem

sim('PID_z_odsprz', T_sim)

figure(3);
subplot(2, 1, 1);
hold on;
title({'PID z odsprzęganiem'});
stairs(T_zad.Time, T_zad.Data, 'Color', 'blue');
plot(T_out.Time, T_out.Data, 'Color', 'red');
xlabel('czas symulacji') ;
ylabel('T');
legend('T zad','T');

subplot(2, 1, 2);
hold on;
stairs(C_a_zad.Time, C_a_zad.Data, 'Color', 'blue');
plot(C_a_out.Time, C_a_out.Data, 'Color', 'red');
xlabel('czas symulacji') ;
ylabel('Ca');
legend('Ca zad','Ca');

%%