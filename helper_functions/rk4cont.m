function [y, t] = rk4cont(data, sp)
    t = 0:sp:10;
%     Ca = CaI;
%     T = TI;
%     y(:, 1) = [Ca T];
  u = [data * ones(size(t));zeros(size(t))];
A = [ -1.118, -0.01601; 15.3, -0.7927 ];

B = [ 1, 0; 0, -1.394 ];

C = [ 1, 0; 0, 1 ];

D = [ 0, 0; 0, 0 ];

sys = ss(A,B,C,D);
y = lsim(sys, u, t,[1.7895, 331.0083]); 
%     for i = 1:(length(t) - 1)
%         k11 = dCa(Ca, T);
%         k12 = dT(Ca, T);
% 
%         k21 = dCa(Ca + 0.5 * sp * k11, T + 0.5 * sp * k11);
%         k22 = dT(Ca + 0.5 * sp * k12, T + 0.5 * sp * k12);
% 
%         k31 = dCa(Ca + 0.5 * sp * k21, T + 0.5 * sp * k21);
%         k32 = dT(Ca + 0.5 * sp * k22, T + 0.5 * sp * k22);
% 
%         k41 = dCa(Ca + sp * k31, T + sp * k31);
%         k42 = dT(Ca + sp * k32, T + sp * k32);
% 
%         Ca = Ca + (sp / 6) * (k11 + k41 + 2 * (k21 + k31));
%         T = T + (sp / 6) * (k12 + k42 + 2 * (k22 + k32));
%         y(:, i + 1) = [Ca T];
%         %Ca
%     end

end




