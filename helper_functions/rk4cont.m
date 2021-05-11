function [y, t] = rk4cont(data, sp)
    t = 0:sp:10;

  u = [data * ones(size(t));zeros(size(t))];
A = [ -1.118, -0.01601; 15.3, -0.7927 ];

B = [ 1, 0; 0, -1.394 ];

C = [ 1, 0; 0, 1 ];

D = [ 0, 0; 0, 0 ];

sys = ss(A,B,C,D);
y = lsim(sys, u, t,[1.7895, 331.0083]); 

end




