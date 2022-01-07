
% Control system model
Ts = 2;
Gp_c = tf(1, [15, 1], 'InputDelay', 2)

Gp_d = c2d(Gp_c, Ts, 'zoh')
% pzmap(Gp)
% pzmap(Gp_d)


% Reference model
ks=12; Tn=12;
Gref = tf(ks, [Tn/2.2, 1])
Gref_d = c2d(Gref, Ts, 'zoh')
% pzmap(Gref)
% pzmap(Gref_d)


% Polynomials
z = tf('z');

[B, A] = tfdata(Gp_d, 'v');
A0 = 1;
[Bm, Am] = tfdata(Gref_d, 'v');

F = 1; % ?

a1 = A(2);
am1 = Am(2);
g0 = am1-a1;

s0 = g0;
r1 = B(2);
t0 = Bm(2);

% G = Am*A0 - A*F;
% 
% S = G;
% R = B * F;
% T = A0 * Bm;
% step(Gp_c)
% hold on
% step(Gp_d)


