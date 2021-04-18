function [y,t] = rk4Discrete(dCa,dT,Ca,T,sp,Ts)
t=0:sp:10;
y(:,1) = [Ca T];
acc = [0 0];
for i=1:(length(t)-1)
    k11=dCa(Ca,T);
    k12=dT(Ca,T);

    k21=dCa(Ca+0.5*sp,T+0.5*sp*k11);
    k22=dT(Ca+0.5*sp,T+0.5*sp*k12);

    k31=dCa(Ca+0.5*sp,T+0.5*sp*k21);
    k32=dT(Ca+0.5*sp,T+0.5*sp*k22);

    k41=dCa(Ca+sp,T+sp*k31);
    k42=dT(Ca+sp,T+sp*k32);

    Ca=Ca+(sp/6)*(k11+k41+2*(k21+k31));
    T=T+(sp/6)*(k12+k42+2*(k22+k32));
    
    if(mod(i-1,Ts) == 0)
        acc = [Ca T];
    end
    y(:,i+1)=acc;
end
end