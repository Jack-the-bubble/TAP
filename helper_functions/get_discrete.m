function [cont_ss, discr_ss, discr_tf] = get_discrete(Ts)
%GET_DISCRETE Wyznacz modele dyskretne obiektu reaktora przeplywowego
    globals;
    global  ro cp k E_R h a b V Fin Tin F Fc0 Ca0 T0;
    
    A = [ -(F + V*k*exp(-E_R/T0))/V, -(Ca0*E_R*k*exp(-E_R/T0))/T0^2;
        (h*k*exp(-E_R/T0))/(cp*ro), -((a*exp(log(Fc0)*(b + 1)))/(Fc0 ...
            + (a*exp(b*log(Fc0)))/(2*cp*ro)) + F*cp*ro - (Ca0*E_R*V*h*k*exp(-E_R/T0))/T0^2)/(V*cp*ro)];
    B = [ Fin/V, 0;
        0, (a*(T0 - Tin)*((exp(log(Fc0)*(b + 1))*((a*b*exp(b*log(Fc0)))/(2*Fc0*cp*ro) + 1)) ... 
            /(Fc0 + (a*exp(b*log(Fc0)))/(2*cp*ro))^2 - (exp(log(Fc0)*(b + 1))*(b + 1))/(Fc0*(Fc0 ...
            + (a*exp(b*log(Fc0)))/(2*cp*ro)))))/(V*cp*ro)];
    C = [1 0; 0 1];
    D = zeros(2,2); % Układ z samymi wspolczynnikami rzeczywistymi - macierz zerowa
    
%   continuous model
    cont_ss = ss(A,B,C,D);

%   discrete model
    method = 'zoh'
    discr_ss = c2d(cont_ss, Ts, method);
    discr_tf = tf(discr_ss); 
    
end

