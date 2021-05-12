clear all;
close all;

%% ciągły model przestrzeni stanu, wyznaczony w pliku plot_nl.m
A=[-1.118 -0.016;
   15.285 -0.792];

B=[1.0  0;
   0    -1.394];

C=[1 0;
   0 1];

D=[0 0;
    0 0];

% punkt pracy obliczony w czesci 1 projektu
Ca0 = 1.7895;
T0 = 331.0083;
Cain0 = 2;
Fc0 = 15;

%% dyskretny model przestrzeni stanu
Tp=0.01;
ss_c = ss(A, B, C, D);
ss_d = c2d(ss_c, Tp, 'zoh');

[Ad, Bd, Cd, Dd] = ssdata(ss_d);

%% wyznaczenie transmitancji
Tstep=10;
t=0:Tp:(Tstep-Tp);


sys_c=tf(ss(A,B,C,D));
sys_d=tf(ss(Ad,Bd,Cd,Dd,Tp));

%% konfiguracja liczby wyjść i sygnałów sterujących, sygnał wartości zadanej, parametry symulacji
ny=2;
nu=2;
N=200;
Nu=2;
Dyn=length(t)/2;

maxU1 = 0.5;
minU1 = -1.9;
minU2 = -5;

Tsim=30;
t_sim=0:Tp:(Tsim-Tp);

U_past=zeros((Dyn-1)*nu,1);
U=zeros(nu,length(t_sim));

offset = ones(1, length(t_sim));
Y=zeros(2,length(t_sim));
Yzad=zeros(2,length(t_sim));

Ca_zad_progi = [100 600 1200 1800 2400];
T_zad_progi = [200 700 1300 1900 2500];

Ca_zad = [0.4 -0.2 0.6 -0.5 -1.7];
T_zad = [3 -0.5 20 -4 -40];



for i=1:length(Ca_zad)
    Yzad(1, Ca_zad_progi(i):end) = Ca_zad(i);
    Yzad(2, T_zad_progi(i):end) = T_zad(i);
end

%% odpowiedz skokowa i macierze potrzebne do uzyskania prawa sterowania

Y_step=step(sys_d,t);
S=zeros(ny,Dyn*nu);

for i=1.0:1.0:Dyn
    for j=1:1:nu
        for k=1:1:ny
        S(k,nu*(i-1)+j)=Y_step(i,k,j);
        end
    end
end

M=zeros(N*ny,nu*Nu);
Mp=zeros(N*ny,nu*(Dyn-1));

for i=1:1:Nu
    for j=1:1:N
        for k=1:1:nu
            for l=1:1:ny
                if  ( (j-1)*nu+k-(i-1)*nu <= 0 )
                    M(ny*(j-1)+l,k+(i-1)*nu )=0;
                else    
                    M(ny*(j-1)+l,k+(i-1)*nu )=S(l,(j-1)*nu+k-(i-1)*nu);
                end
        
            end
        end
    end
end

for i=1:1:N   
    for j=1:1:(Dyn-1)
        for k=1:1:nu
            for l=1:1:ny
                if(k+(i+j-1)*nu<=Dyn*nu)
                    Mp(ny*(i-1)+l,nu*(j-1)+k)=S(l,k+(i+j-1)*nu)-S(l,k+(j-1)*nu);
                else
                    Mp(ny*(i-1)+l,nu*(j-1)+k)=S(l,k+Dyn*(nu-1))-S(l,k+(j-1)*nu); 
                end
            end
        end
    end    
end    

Lambda=eye(Nu*nu,Nu*nu);
Fi=eye(N*ny,N*ny);

K = (M'*Fi*M+Lambda);
K= (K^-1)*M'*Fi;

Ku=K(1:nu,1:N*ny);

Ke=zeros(nu,2);
for i=1:1:N
    Ke=Ke+K(1:nu,(2*(i-1)+1):(2*(i-1)+2));
end


 
%% DMC w wersji analitycznej

for i=1:1:(length(t_sim)-1)
    
    e=Yzad(:,i)-Y(:,i);
    dU=Ke*e-Ku*Mp*U_past;
    
    for f=1:1:(Dyn-2)
        for j=1:1:nu
            U_past(nu*(Dyn-2)-nu*(f-1)+j)=U_past(nu*(Dyn-2)-nu*(f)+j);
        end
    end
    
    Uprev=U(:, i);
    U(:, i+1)=U(:, i)+dU;
    
    % ograniczenie przez rzutowanie
    if U(1,i+1)>maxU1
        U(1,i+1)=maxU1;
    end
    
    if U(1,i+1)<minU1
        U(1,i+1)=minU1;
    end
    
     if U(2,i+1)<minU2
        U(2,i+1)=minU2;
    end
        
    U_past(1:nu)=U(:, i+1)-Uprev;    
    Y(:,i+1)=Ad*Y(:,i)+Bd*(U(:, i+1));
end

%% wykresy

figure(1)
Yzad(1,:) = Yzad(1, :)+offset*Ca0;
plot(t_sim,Yzad(1,:), 'Color', 'blue');
hold on

Y(1,:) = Y(1, :)+offset*Ca0;
stairs(t_sim,Y(1,:), 'Color', 'red');

title({'DMC analityczny wyjscie Ca'});
xlabel('czas symulacji') ;
ylabel('Ca=y1');
legend('y\_zad','y');


figure(2)
Yzad(2,:) = Yzad(2, :)+offset*T0;
plot(t_sim,Yzad(2,:), 'Color', 'blue');

hold on
Y(2,:) = Y(2, :)+offset*T0;
stairs(t_sim,Y(2,:), 'Color', 'red');
title({'DMC analityczny wyjscie T'});
xlabel('czas symulacji');
ylabel('T=y2');
legend('T\_zad','T');

figure(3)
subplot(2, 1, 1)
U(1, :) = U(1, :)+offset*Cain0;
plot(t_sim,U(1,:), 'Color', 'green');

title({'DMC analityczny przebieg sterowania'});
xlabel('czas symulacji') ;
ylabel('Cain=u1');
legend('Cain');

subplot(2, 1, 2)
U(2, :) = U(2, :)+offset*Fc0;
plot(t_sim,U(2,:), 'Color', 'magenta');

xlabel('czas symulacji') ;
ylabel('Fc=u2');
legend('Fc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%DMC_numerical simualation%%%%%%%%%%%%%%%
Alpha=1.0*eye(Nu*nu,Nu*nu);

Fi=1*eye(N*ny,N*ny);


Xzad_n=zeros(N*ny,1);
Xn = zeros(N*ny,1);

Xmax = ones(N*ny,1)*120;
Xmin = ones(N*ny,1)*-120;

Uhist=zeros((Dyn-1)*nu,1);
dU=zeros(nu*Nu,1);



H=2*(M'*Fi*M+Alpha);

I = eye(nu,nu);
J = zeros(Nu*nu,Nu*nu);
for i=1:1:Nu
    for j=1:1:Nu
        if(j-i<=0)
        
        for f=1:1:nu
            for k=1:1:nu
                J(nu*(i-1)+f,nu*(j-1)+k)=I(f,k);
%                 J(4*(i-1)+f,4*(j-1)+k)=I(f,k);
            end
        end
        
        end
    end
end

Aqp = [-J; J; -M; M];



Tsim=30;
t_sim=0:Tp:(Tsim-Tp);

dUmax=ones(nu*Nu,1)*200;
dUmax(1,1)=0.3;

dUmin=ones(nu*Nu,1)*-200;
dUmin(2,1)=-3.0;

Umin=ones(nu*Nu,1)*-100;
Umin(2,1)=-3;

Umax=ones(nu*Nu,1)*100;
Umax(1,1)=0.3;


Xout=zeros(ny, length(t_sim));
Xzad_out=zeros(ny, length(t_sim));
U=zeros(nu,length(t_sim));
Uutil=zeros(Nu*nu,length(t_sim));


for i=1:1:(length(t_sim)-1)
   
    if(i==100)
       for c=1:1:N
           Xzad_n((c-1)*ny+1,1)=0.2;
       end
    end
    
     X0=Xn+Mp*Uhist;
     
     F=(-2*M'*Fi*(Xzad_n-X0));
   
     Uutil_del=zeros(Nu*nu,1); 
     if(i>1)
         Uutil_del=Uutil(:,i-1);
     end
     
     b = [-Umin+Uutil_del;Umax-Uutil_del;-Xmin+X0;Xmax-X0;];   
     
     dU=quadprog(H,F,Aqp,b,[],[],dUmin,dUmax);
          
     if(i>1)
     U(:,i)=U(:,(i-1))+dU(1:nu,1); 
     Uutil(:,i)=Uutil(:,(i-1))+dU(:,1);
     else
     U(:,i)=dU(1:nu,1); 
     Uutil(:,i)=dU(:,1);    
     end
     
     %Disturbance
     %U(3,i)=0;
     %U(4,i)=-16;
     
     for f=1:1:(Dyn-2)
         for j=1:1:nu
             Uhist(nu*(Dyn-2)-nu*(f-1)+j,1)=Uhist(nu*(Dyn-2)-nu*(f)+j,1);
         end
     end

    Uhist(1:nu,1)=dU(1:nu,1);

    if(i>1) 
    Xn(1:ny)=Xout(:,(i-1));
    end 
    
    Ualg=U(:,i);
    
    for j=1:1:(N-1)
        
            
        
            if( (j<=Nu) && (j>1))
                Ualg=Ualg+dU(nu*(j-1)+1:nu*(j),1);
            end
            
        Xn( (j)*ny+1:(j+1)*ny,1)=Ad*Xn( (j-1)*ny+1:(j)*ny,1)+Bd*Ualg;
    end        

    Xout(:,i)=Xn(1*ny+1:2*ny,1);
    Xzad_out(:,i)=Xzad_n(1:ny,1);
end

figure(4)
plot(t_sim,Yzad(1,:));
hold on
stairs(t_sim,Xout(1,:));
title({'DMC numerical Ca control'})
xlabel('t') % x-axis label
ylabel('Ca=y1') % y-axis label
legend('Ca','Ca zad')
print('DMC_num_Ca','-djpeg','-r 300')

figure(5)
plot(t_sim,Yzad(2,:));
hold on
stairs(t_sim,Xout(2,:));
title({'DMC numerical T control'})
xlabel('t') % x-axis label
ylabel('T=y2') % y-axis label
legend('T','T zad')
print('DMC_num_T','-djpeg','-r 300')


U(1, :) = U(1, :)+offset*Cain0;
U(2, :) = U(2, :)+offset*Fc0;

figure(6)
subplot(2, 1, 1)
plot(t_sim,U(1,:), 'Color', 'green');

title({'DMC numeryczny przebieg sterowania'});
xlabel('czas symulacji') ;
ylabel('Cain=u1');
legend('Cain');

subplot(2, 1, 2)
plot(t_sim,U(2,:), 'Color', 'magenta');

xlabel('czas symulacji') ;
ylabel('Fc=u2');
legend('Fc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Comparing plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(6)
% subplot(2,2,1)
% stairs(t_sim,Xout(1,:));
% hold on
% plot(t_sim,Xzad_out(1,:));
% title({'DMC numerical Ca control'})
% xlabel('t') % x-axis label
% ylabel('Ca=y1') % y-axis label
% legend('Ca','Ca zad')