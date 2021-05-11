clear all;
close all;

%%%%%%%%%%%%%%Continious model%%%%%%%%%%%%%%%%
A=[-1.1176 -0.0160;
   15.2857 -0.7933];

B=[1.0  0;
   0    -1.3935];

C=[1 0;
   0 1];

D=[0 0;
    0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Discrete_model%%%%%%%%%%%%%%%%%%%%
Tp=0.01;

Ad=[0.9889     -0.0001585; 
    0.1514     0.9921];

Bd=[0.009944   1.108e-06   -7.95e-07   -1.49e-06;
    0.0007594  -0.01388    0.00996      0.01866  ];

Cd=[1  0;
    0  1];

Dd=[0 0 0 0;
    0 0 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Discrete continious transfer model%%%%%%%%%
Tfinal=6;
t=0:Tp:(Tfinal-Tp);
Dyn=length(t);

sys_c=tf(ss(A,B,C,D));
sys_d=tf(ss(Ad,Bd,Cd,Dd,0.01));

ny=2;
nu=4;


Y=step(sys_d,t);
S=zeros(ny,Dyn*nu);

for i=1.0:1.0:Dyn
    for j=1:1:nu
        for k=1:1:ny
        S(k,nu*(i-1)+j)=Y(i,k,j);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%Params for DMC%%%%%%%%%%%%%%%%%%%%%%%%%
N=200;
Nu=2;
M=zeros(N*ny,nu*Nu);
Mp=zeros(N*ny,nu*(Dyn-1));

%%%M matrix%%%
for f=1:1:Nu
for i=1:1:N
    
    for j=1:1:nu
        
        for k=1:1:ny
           
        if  ( (i-1)*nu+j-(f-1)*nu <= 0 )
             M(ny*(i-1)+k,j+(f-1)*nu )=0;
        else    
             M(ny*(i-1)+k,j+(f-1)*nu )=S(k,(i-1)*nu+j-(f-1)*nu);
        end
        
        end
        
        
    
    end
     
end
end
%%%Mp matrix%%%

for i=1:1:N
   
    for j=1:1:(Dyn-1)
    
        for k=1:1:nu
            
            for f=1:1:ny
                if(k+(i+j-1)*nu<=Dyn*nu)
                Mp(ny*(i-1)+f,nu*(j-1)+k)=S(f,k+(i+j-1)*nu)-S(f,k+(j-1)*nu);
                else
                Mp(ny*(i-1)+f,nu*(j-1)+k)=S(f,k+Dyn*(nu-1))-S(f,k+(j-1)*nu); 
                end
            end
        end
  
    end    
    
end    
%%%K matrix%%%
Alpha=eye(Nu*nu,Nu*nu);
Fi=eye(N*ny,N*ny);

K=inv((M'*Fi*M+Alpha))*M'*Fi;

Ke=zeros(nu,2);
% Ke=zeros(4,2);
for i=1:1:N
Ke=Ke+K(1:nu,(2*(i-1)+1):(2*(i-1)+2));
end

Ku=K(1:nu,1:N*ny);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%DMC_analitical simualation%%%%%%%%%%%%%%
% Uhist=zeros((Dyn-1)*nu,1);
% 
% 
% U=zeros(4,1);
% 
% Tsim=8;
% t_sim=0:Tp:(Tsim-Tp);
% 
% X=zeros(2,length(t_sim));
% 
% Xzad=zeros(2,length(t_sim));
% 
% %%Disturbance%%
% Tin_dist=0.0;
% TCin_dist=0.0;%-16.0;
% maxU1 = 0.5;
% minU2 = -5;
% 
% for i=1:1:(length(t_sim)-1)
%     
%     
%     if(i>=100)
%         Xzad(:,i)=[0.4;0.0];
%     end    
%     
%     e=Xzad(:,i)-X(:,i);
%     dU=Ke*e-Ku*Mp*Uhist;
%     
%     for f=1:1:(Dyn-2)
%         for j=1:1:4
%             Uhist(nu*(Dyn-2)-nu*(f-1)+j)=Uhist(nu*(Dyn-2)-nu*(f)+j);
%         end
%     end
%     
%     Uprev=U;
%     U=U+dU;
%     
%     %%%Seting uncontrolled u%%%
%     U(3,1)=0+Tin_dist;
%     U(4,1)=0+TCin_dist;
%     
%     %%%Limits for input%%%
%     if U(1,1)>maxU1
%         U(1,1)=maxU1;
%     end
%     
%      if U(2,1)<minU2
%         U(2,1)=minU2;
%     end
%     
%         
%     Uhist(1:4)=U-Uprev;
%      
%     
%     X(:,i+1)=Ad*X(:,i)+Bd*(U);
%     
% 
% end
% 
% figure(1)
% stairs(t_sim,X(1,:));
% hold on
% plot(t_sim,Xzad(1,:));
% title({'DMC Ca control'})
% xlabel('t') % x-axis label
% ylabel('Ca=y1') % y-axis label
% legend('Ca','Ca zad')
% print('DMC_anal_Ca','-djpeg','-r 300')
% 
% figure(2)
% stairs(t_sim,X(2,:));
% hold on
% plot(t_sim,Xzad(2,:));
% title({'DMC T control'})
% xlabel('t') % x-axis label
% ylabel('T=y2') % y-axis label
% legend('T','T zad')
% print('DMC_anal_T','-djpeg','-r 300')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%DMC_numerical simualation%%%%%%%%%%%%%%%
% 
% Xzad_n=zeros(N*ny,1);
% Xn = zeros(N*ny,1);
% 
% Xmax = ones(N*ny,1)*200;
% Xmin = ones(N*ny,1)*-200;
% 
% Uhist=zeros((Dyn-1)*nu,1);
% dU=zeros(nu*Nu,1);
% 
% Umin=ones(nu*Nu,1)*-100;
% Umax=ones(nu*Nu,1)*100;
% 
% 
% H=2*(M'*Fi*M+Alpha);
% 
% I = eye(nu,nu);
% J = zeros(Nu*nu,Nu*nu);
% for i=1:1:Nu
%     for j=1:1:Nu
%         for f=1:1:nu
%             for k=1:1:nu
%                 J(4*(i-1)+f,4*(j-1)+k)=I(f,k);
%             end
%         end
%     end
% end
% 
% Aqp = [-J; J; -M; M];
% 
% 
% 
% Tsim=8;
% t_sim=0:Tp:(Tsim-Tp);
% 
% dUmax=ones(nu*Nu,1)*1;
% dUmin=ones(nu*Nu,1)*-1;
% 
% Xout=zeros(ny, length(t_sim));
% Xzad_out=zeros(ny, length(t_sim));
% U=zeros(nu,length(t_sim));
% Uutil=zeros(Nu*nu,length(t_sim));
% 
% 
% for i=1:1:(length(t_sim)-1)
%    
%     if(i==100)
%        for c=1:1:N
%            Xzad_n((c-1)*ny+1,1)=0.2;
%        end
%     end
%     
%      X0=Xn+Mp*Uhist;
%      
%      f=(-2*M'*Fi*(Xzad_n-Xn));
%    
%      Uutil_del=zeros(Nu*nu,1); 
%      if(i>1)
%          Uutil_del=Uutil(:,i-1);
%      end
%      
%      b = [-Umin+Uutil_del;Umax-Uutil_del;-Xmin+X0;Xmax-X0;];   
%      
%      dU=quadprog(H,f,Aqp,b,[],[],dUmin,dUmax);
%           
%      U(:,i)=U(:,i)+dU(1:nu,1); 
%      Uutil(:,i)=Uutil(:,i)+dU(:,1);
%      
%      Uhist(1:nu)=U(:,i);
%      
%      for f=1:1:(Dyn-2)
%          for j=1:1:4
%              Uhist(nu*(Dyn-2)-nu*(f-1)+j)=Uhist(nu*(Dyn-2)-nu*(f)+j);
%          end
%      end
% 
%     for j=1:1:(N-1)
%         
%             Ualg=U(:,i);
%         
%             if( (j<=Nu) && (j>1))
%                 Ualg=Ualg+dU(nu*(j-1)+1:nu*(j),1);
%             end
%             
%         Xn( (j)*ny+1:(j+1)*ny,1)=Ad*Xn( (j-1)*ny+1:(j)*ny,1)+Bd*Ualg;
%     end        
% 
%     Xout(:,i)=Xn(ny+1:2*ny,1);
%     Xzad_out(:,i)=Xzad_n(1:ny,1);
% end
% 
% figure(1)
% stairs(t_sim,Xout(1,:));
% hold on
% plot(t_sim,Xzad_out(1,:));
% title({'DMC numerical Ca control'})
% xlabel('t') % x-axis label
% ylabel('Ca=y1') % y-axis label
% legend('Ca','Ca zad')
% print('DMC_num_Ca','-djpeg','-r 300')
% 
% figure(2)
% stairs(t_sim,Xout(2,:));
% hold on
% plot(t_sim,Xzad_out(2,:));
% title({'DMC numerical T control'})
% xlabel('t') % x-axis label
% ylabel('T=y2') % y-axis label
% legend('T','T zad')
% print('DMC_num_T','-djpeg','-r 300')
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



Tsim=8;
t_sim=0:Tp:(Tsim-Tp);

dUmax=ones(nu*Nu,1)*200;
dUmax(1,1)=0.3;

dUmax(3,1)=0;
% dUmax(7,1)=0;

dUmax(4,1)=0;
% dUmax(8,1)=0;


dUmin=ones(nu*Nu,1)*-200;
dUmin(2,1)=-3.0;

dUmin(3,1)=0;
% dUmin(7,1)=0;

dUmin(4,1)=0;
% dUmin(8,1)=0;


Umin=ones(nu*Nu,1)*-100;
Umin(2,1)=-3;

Umin(3,1)=0;
% Umin(7,1)=0;

Umin(4,1)=0;
% Umin(8,1)=0;



Umax=ones(nu*Nu,1)*100;
Umax(1,1)=0.3;

Umax(3,1)=0;
% Umax(7,1)=0;

Umax(4,1)=0;
% Umax(8,1)=0;



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

figure(1)
stairs(t_sim,Xout(1,:));
hold on
plot(t_sim,Xzad_out(1,:));
title({'DMC numerical Ca control'})
xlabel('t') % x-axis label
ylabel('Ca=y1') % y-axis label
legend('Ca','Ca zad')
print('DMC_num_Ca','-djpeg','-r 300')

figure(2)
stairs(t_sim,Xout(2,:));
hold on
plot(t_sim,Xzad_out(2,:));
title({'DMC numerical T control'})
xlabel('t') % x-axis label
ylabel('T=y2') % y-axis label
legend('T','T zad')
print('DMC_num_T','-djpeg','-r 300')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%