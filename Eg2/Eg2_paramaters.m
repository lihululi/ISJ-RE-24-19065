function [Ad,Bd,Cd,F,Gd,kappa,T,mu,dd,c1,c2,Attack_signal,Channel_attacked,x00] = a1paramaters()

Tt1 = 0.33; Tt2 = 0.35;
Tg1 = 0.1; Tg2 = 0.1;
R1 = 0.2; R2 = 0.22;
H1 = 3.5; H2 = 2.5;
D1 = 3; D2 = 4;
P12 = 8; P21 = 8;

A11 = [0,1,0,0;
       -P12/(2*H1), -D1/(2*H1), 1/(2*H1), 0;
       0, 0, -1/Tt1, 1/Tt1;
       0, -1/(R1*Tg1), 0, -1/Tg1];

A22 = [0,1,0,0;
       -P21/(2*H2), -D2/(2*H2), 1/(2*H2), 0;
       0, 0, -1/Tt2, 1/Tt2;
       0, -1/(R2*Tg2), 0, -1/Tg2];

A12 = [0,0,0,0;
    P12/(2*H1),0,0,0;
    0,0,0,0;
    0,0,0,0];

A21 = [0,0,0,0;
    P21/(2*H2),0,0,0;
    0,0,0,0;
    0,0,0,0];

Ac = [A11, A12;
     A21, A22];

B1 = [0;0;0;1/Tg1];
B2 = B1;

Bc = blkdiag(B1,B2);

C1 = 0.02*[5, 1, 0, 0;
          5, 0, 2, 1];
C2 = 0.02*[5, 1, 0, 0;
          5, 0, 2, 1];

Cc = blkdiag(C1,C2);

G1 = [0;-1/(2*H1);0;0];
G2 = G1;

Gc = blkdiag(G1,G2);

D = zeros(4,2);

dt = 4;

S = ss(Ac,Bc,Cc,D);
XX = c2d(S,dt);
Ad = XX.A;
Bd = XX.B;
Cd = XX.C;

F = [0.1,0,0,0,0.1,0,0,0];

S = ss(Ac,Gc,Cc,D);
XX = c2d(S,dt);
Gd = XX.B;

T = 28;
mu = 0.01;
c1 = 1;
c2 = 9;
dd = 0.25;

R = eye(13,13);
kappa = 1.1;

Attack_signal1 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
Attack_signal2 = [0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0];
Attack_signal3 = [0,0,1,1,1,1,1,1,0,0,0,1,1,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0];
Attack_signal4 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0];

Attack_signal = [Attack_signal1;Attack_signal2;Attack_signal3;Attack_signal4];

Channel_attacked = [1];
for i = 1:28
    if Attack_signal(1,i+1) == 1
        Channel_attacked = [Channel_attacked,1];
    elseif Attack_signal(2,i+1) == 1
        Channel_attacked = [Channel_attacked,2];
    elseif Attack_signal(3,i+1) == 1
        Channel_attacked = [Channel_attacked,3];
    elseif Attack_signal(4,i+1) == 1
        Channel_attacked = [Channel_attacked,4];
    else
        Channel_attacked = [Channel_attacked,0];
    end
end

x00 = 0.5*[0.8,0.9,0.7,1.1,1,0.6,0.7,1]';
end