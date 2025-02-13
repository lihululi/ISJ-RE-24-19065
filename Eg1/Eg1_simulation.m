clc
clear
close all

[A,B,C,F,G,kappa,T,mu,dd,c1,c2,Attack_signal,Channel_attacked,x00,R] = Eg1_parameters();

bar_I1 = [eye(2,2);zeros(2,2)];
bar_I2 = [zeros(2,2);eye(2,2)];

Phi1 = [1,0;0,0];
Phi2 = [0,0;0,1];
Psi1 = [1,0;0,0];
Psi2 = [0,0;0,1];

rho = 0;
Phi_hat_1 = (1-rho)*bar_I1*Phi1*bar_I1';
Phi_bar_1 = (1-rho)*bar_I1*(eye(2,2)-Phi1)*bar_I1';

Phi_hat_2 = (1-rho)*bar_I1*Phi2*bar_I1';
Phi_bar_2 = (1-rho)*bar_I1*(eye(2,2)-Phi2)*bar_I1';

rho = 1;
Phi_hat_3 = rho*bar_I2*Psi1*bar_I2';
Phi_bar_3 = rho*bar_I2*(eye(2,2)-Psi1)*bar_I2';

Phi_hat_4 = rho*bar_I2*Psi2*bar_I2';
Phi_bar_4 = rho*bar_I2*(eye(2,2)-Psi2)*bar_I2';

% obtained by Eg1_Algorithm3.m
XXX = [-0.831634990075268	-0.0838664615798684	-0.00124103137883746	-0.0177165330126700	0.00192785107902779	1.24383353286016e-06	0.00188145039663906	-0.000725910300711874	1.54331735512823e-05	-0.000306664200517386	0.000393339508809670	0.000351454213316612	0.0300000000000000	0.0434910319767313	0.0302412605412358	0.0340430406117073	-0.223928507370471	-8.75351357839080e-09	-0.0119614799922210	-0.465185338431976	-0.452900784611884	-0.00430343943243012	0.0164434859486846	-0.693314982741930	0.00138784421593914	0.00787192608846375];
Ac1 = XXX(1);
Ac2 = XXX(2);
Ac3 = XXX(3);
Ac4 = XXX(4);

Bc1 = [XXX(5:6),0,0];
Bc2 = [XXX(7:8),0,0];
Bc3 = [0,0,XXX(9:10)];
Bc4 = [0,0,XXX(11:12)];

Dc1 = XXX(13);
Dc2 = XXX(14);
Dc3 = XXX(15);
Dc4 = XXX(16);

Ec1 = [XXX(17:18),0,0];
Ec2 = [XXX(19:20),0,0];
Ec3 = [0,0,XXX(21:22)];
Ec4 = [0,0,XXX(23:24)];

Ac = [Ac1,Ac2,Ac3,Ac4];
Bc = [Bc1,Bc2,Bc3,Bc4];
Dc = [Dc1,Dc2,Dc3,Dc4];
Ec = [Ec1,Ec2,Ec3,Ec4];

barQ = diag([XXX(25),XXX(26)]);

x0 = x00;
y0 = C*x0;
baryt = zeros(4,1);

State = [x0];
Y = [baryt];

disturbance = [];
for k = 1:T+1
    disturbance = [disturbance,1.5*(sin(k)^2)/(1+k)];
end

for k = 1:T
    xk_1 = A*x0 + G*disturbance(1,k);
    xk = xk_1;
    x0 = xk;
    if mod(k,4)+1 < 3
        if mod(k,4)+1 == 1
            Phi = Phi_hat_1;
            barPhi = Phi_bar_1;
            channelt = 1;
        else
            Phi = Phi_hat_2;
            barPhi = Phi_bar_2;
            channelt = 2;
        end
    else
        y_qta_3 = y0(3) - baryt(3);
        tod3 = y_qta_3'*barQ(1,1)*y_qta_3;
        y_qta_4 = y0(4) - baryt(4);
        tod4 = y_qta_4'*barQ(2,2)*y_qta_4;
        if tod3 >= tod4
            Phi = Phi_hat_3;
            barPhi = Phi_bar_3;
            channelt = 3;
        else
            Phi = Phi_hat_4;
            barPhi = Phi_bar_3;
            channelt = 4;
        end
    end
    baryt = Phi*y0 + barPhi*baryt;
    y0 = C*x0;
    
    State = [State,x0];
    Y = [Y,baryt];
end

xx = [State(:,T+1);0;Y(:,T+1)];

x0 = [State(:,1);0;Y(:,1)];
x0'*R*x0

xx'*R*xx

x0 = x00;
v0 = 0;
u0 = 0;

channel0 = 1;
y0 = C*x0;
baryt = zeros(4,1);

C_State = [x0];
U = [u0];
V = [v0];
C_Y = [baryt];

Channel_selected = [channel0];

l = 2;
ny = 4;

attack = 0;
for k = 0:T-1
    if mod(k,4)+1 < 3 
        for step = 0:1
            if Attack_signal(mod(k+attack+step,l)+1, k+2) == 0
                attack = attack + step;
                break;
            end
        end
        if mod(k+attack,l)+1 == 1
            Phi = Phi_hat_1;
            barPhi = Phi_bar_1;
            channelt = 1;
        else
            Phi = Phi_hat_2;
            barPhi = Phi_bar_2;
            channelt = 2;
        end
    else
        Mk = 0;
        minW = inf;
        maxW = -1;

        Ae = [-1];
        Se = [10000];
        for i = l+1:ny
            y_qta = y0(i) - baryt(i);
            tod = y_qta'*barQ(i-l,i-l)*y_qta;
            if Attack_signal(i, k+2) == 0
                Se = [Se, tod];
            else
                Ae = [Ae, tod];
            end
        end

        if min(Se) > max(Ae)
            Mk = 0;
        else
            Mk = max(Ae) - min(Se) + 1e-9;
        end
        
        WTOD = [];
        for i = l+1:ny
            y_qta = y0(i) - baryt(i);
            tod = y_qta'*barQ(i-l,i-l)*y_qta;
            if Attack_signal(i, k+2) == 0
                WTOD = [WTOD,tod+Mk];
            else
                WTOD = [WTOD,tod];
            end
        end

        [t,index] = sort(WTOD, 'descend');
        
        if index(1) == 1
            Phi = Phi_hat_3;
            barPhi = Phi_bar_3;
            channelt = 3;
        else
            Phi = Phi_hat_4;
            barPhi = Phi_bar_4;
            channelt = 4;
        end
    end
    baryt = Phi*y0 + barPhi*baryt;
    for i = 1:4
        if Attack_signal(i, k+2) == 1
            baryt(i,1) = 0;
        end
    end

    if channelt == 1
        Dc1 = Dc(1,1);
        Ec1 = Ec(1,1:4);
    elseif channelt == 2
        Dc1 = Dc(1,2);
        Ec1 = Ec(1,5:8);
    elseif channelt == 3
        Dc1 = Dc(1,3);
        Ec1 = Ec(1,9:12);
    else
        Dc1 = Dc(1,4);
        Ec1 = Ec(1,13:16);
    end

    if channel0 == 1
        Ac1 = Ac(1,1);
        Bc1 = Bc(1,1:4);
    elseif channel0 == 2
        Ac1 = Ac(1,2);
        Bc1 = Bc(1,5:8);
    elseif channel0 == 3
        Ac1 = Ac(1,3);
        Bc1 = Bc(1,9:12);
    else
        Ac1 = Ac(1,4);
        Bc1 = Bc(1,13:16);
    end

    u0 = Dc1*v0 + Ec1*baryt;
    v0 = Ac1*v0 + Bc1*baryt;
    x0 = A*x0 + B*u0 + G*disturbance(1,k+1);
    y0 = C*x0;
    
    channel0 = channelt;
    Channel_selected = [Channel_selected,channel0];

    C_State = [C_State,x0];
    U = [U,u0];
    V = [V,v0];
    C_Y = [C_Y,baryt];
end

Tn = T;
T = (0:T);
set(groot, 'defaultAxesColorOrder', lines);

figure(1)
hold on;
box on;
stairs(T,Channel_selected(:,1:Tn+1),'LineWidth',1.5); 
stairs(T,Channel_attacked(:,1:Tn+1), '--','LineWidth',1.5); 
ylim([0 4.3]);
xlim([0 Tn]);
xlabel('$k$', 'Interpreter', 'latex','FontSize', 15);
ylabel('$\alpha(k)$', 'Interpreter', 'latex','FontSize', 15);
xticks(0:4:Tn);
title('The selected sensor node and the attacked sensor node');


%% FTB Figure
cNorm = [];
Norm = [];

for k = 1:Tn+1
    cxk = C_State(:,k);
    cvk = V(:,k);
    cybark = C_Y(:,k);
    cetak = [cxk;cvk;cybark];

    cnorm = cetak'*R*cetak;
    
    xk = State(:,k);
    ybark = Y(:,k);
    etak = [xk;0;ybark];
    norm = etak'*R*etak;

    cNorm = [cNorm, cnorm];
    Norm = [Norm, norm];
end

figure(2);
z = zeros(1,Tn+1);
plot3(cNorm,T,z,'r','Linewidth',1.5);

hold on;
plot3(Norm,T,z,'b','LineWidth',1.5);
axis([-15 15 0 Tn -15 15])%设置坐标范围

h=0:0.5:Tn;%所画图形的z坐标范围

t=0:pi/Tn:2*pi;

x=zeros(length(h),length(t));

y=zeros(length(h),length(t));

z=zeros(length(h),length(t));

for i = 1:length(h)
    x(i,:) = 9*cos(t);
end
for ii = 1:length(h)
    y(ii,:) = 9*sin(t);
end
for iI = 1:length(h)
    z(iI,:) = h(iI);
end

mesh(x,z,y)
alpha(0.3);
set(gca,'XDir','reverse');
ylabel('k', 'Interpreter', 'latex', 'FontSize', 15);
legend('closed-loop', 'open-loop', '$c_2^2$', 'Interpreter', 'latex', 'FontSize', 15, 'Location', 'Best');
grid on;

sum = 0;
for k = 1:Tn
    sum = disturbance(1,k)'*disturbance(1,k) + sum;
end
