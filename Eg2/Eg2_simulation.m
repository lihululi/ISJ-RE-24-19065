clc
clear
close all

[A,B,C,F,G,kappa,T,mu,dd,c1,c2,Attack_signal,Channel_attacked,x00] = Eg2_paramaters();

omi_s3 = 0.5;
omi_s4 = 0.5;
omi_d3 = 0.5;
omi_d4 = 0.5;

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

% obtained by Eg2_Algorithm3.m
XXX = [-0.861000665218970	0.930875606330750	0.918085615157278	0.628781803672140	-0.251182349093371	0.00954940066406957	0.00140618281486311	-0.427231083915386	-0.187684604862179	0.00323877577568321	0.00177429701186810	-0.287361560736314	0.0874202196159915	0.101166663437518	0.0914735064561277	0.130000000000000	0.130000000000000	0.0751495463315811	0.130000000000000	0.0700000000000000	-13.4294278567931	-13.5977933967535	-0.0476363165059146	0.0518831986788200	-0.0380230784824705	-0.00158711842950169	-16.9178435298751	-17.0564368738765	-14.8658441844881	-15.0910636050351	0.0454865825631624	-0.0299700521682063	-0.0470941324794785	-0.0162485483693548	-18.1113205949319	-18.4251190593226	0.221017027732336	0.0163207922236799];
Ac1 = XXX(1);
Ac2 = XXX(2);
Ac3 = XXX(3);
Ac4 = XXX(4);

Bc1 = [XXX(5:6),0,0];
Bc2 = [XXX(7:8),0,0];
Bc3 = [0,0,XXX(9:10)];
Bc4 = [0,0,XXX(11:12)];

Dc1 = XXX(13:14);
Dc2 = XXX(15:16);
Dc3 = XXX(17:18);
Dc4 = XXX(19:20);

Ec1 = [XXX(21:22)',XXX(23:24)',zeros(2,2)];
Ec2 = [XXX(25:26)',XXX(27:28)',zeros(2,2)];
Ec3 = [zeros(2,2),XXX(29:30)',XXX(31:32)'];
Ec4 = [zeros(2,2),XXX(33:34)',XXX(35:36)'];

Ac = [Ac1,Ac2,Ac3,Ac4];
Bc = [Bc1,Bc2,Bc3,Bc4];
Dc = [Dc1,Dc2,Dc3,Dc4];
Ec = [Ec1,Ec2,Ec3,Ec4];

barQ = diag([XXX(37),XXX(38)]);

x0 = x00;
y0 = C*x0;
baryt = zeros(4,1);
State = [x0];
Y = [baryt];

disturbance = [];
for k = 1:T+1
    disturbance = [disturbance,[(sin(k))/((1+k)*pi);(cos(k))/((1+k)*pi)]];
end

x0 = x00;
v0 = 0;
u0 = zeros(2,1);
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
        Dc1 = Dc(:,1);
        Ec1 = Ec(:,1:4);
    elseif channelt == 2
        Dc1 = Dc(:,2);
        Ec1 = Ec(:,5:8);
    elseif channelt == 3
        Dc1 = Dc(:,3);
        Ec1 = Ec(:,9:12);
    else
        Dc1 = Dc(:,4);
        Ec1 = Ec(:,13:16);
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
    x0 = A*x0 + B*u0 + G*disturbance(:,k+1);
    y0 = C*x0;
    
    channel0 = channelt;
    Channel_selected = [Channel_selected,channel0];

    C_State = [C_State,x0];
    U = [U,u0];
    V = [V,v0];
    C_Y = [C_Y,baryt];
end

% 画图
Tn = T;
T = (0:T);
set(groot, 'defaultAxesColorOrder', lines);
figure(1)
box on;
hold on;
cx1 = C_State(1,:);
cx2 = C_State(2,:);
cx3 = C_State(3,:);
cx4 = C_State(4,:);
x1 = State(1,:);
x2 = State(2,:);
x3 = State(3,:);
x4 = State(4,:);

plot(T, cx1, 'r-', 'LineWidth', 1.1);
plot(T, cx2, 'b-', 'LineWidth', 1.1);
plot(T, cx3, 'c-', 'LineWidth', 1.1);
plot(T, cx4, 'g-', 'LineWidth', 1.1);

plot(T, x1, 'r--', 'LineWidth', 1.1);
plot(T, x2, 'b--', 'LineWidth', 1.1);
plot(T, x3, 'c--', 'LineWidth', 1.1);
plot(T, x4, 'g--', 'LineWidth', 1.1);
xlabel('$k$', 'Interpreter', 'latex','FontSize', 15);
xlim([0,Tn]);

legend('$cx_1$','$cx_2$','$cx_3$', '$cx_4$', '$x_1$', '$x_2$', '$x_3$', '$x_4$', 'Interpreter', 'latex','FontSize', 15);
title('system state of area1');
hold off;

figure(2)
box on;
hold on;
cx1 = C_State(5,:);
cx2 = C_State(6,:);
cx3 = C_State(7,:);
cx4 = C_State(8,:);
x1 = State(5,:);
x2 = State(6,:);
x3 = State(7,:);
x4 = State(8,:);

plot(T, cx1, 'r-', 'LineWidth', 1.1);
plot(T, cx2, 'b-', 'LineWidth', 1.1);
plot(T, cx3, 'c-', 'LineWidth', 1.1);
plot(T, cx4, 'g-', 'LineWidth', 1.1);

plot(T, x1, 'r--', 'LineWidth', 1.1);
plot(T, x2, 'b--', 'LineWidth', 1.1);
plot(T, x3, 'c--', 'LineWidth', 1.1);
plot(T, x4, 'g--', 'LineWidth', 1.1);
xlabel('$k$', 'Interpreter', 'latex','FontSize', 15);
xlim([0,Tn]);

legend('$cx_1$','$cx_2$','$cx_3$', '$cx_4$', '$x_1$', '$x_2$', '$x_3$', '$x_4$', 'Interpreter', 'latex','FontSize', 15);
title('system state of area2');
hold off;

figure(3)
box on;
hold on;
x1 = C_State(1,:);
x2 = C_State(2,:);
x3 = C_State(3,:);
x4 = C_State(4,:);
x5 = C_State(5,:);
x6 = C_State(6,:);
x7 = C_State(7,:);
x8 = C_State(8,:);
plot(T, x1, 'r-', 'LineWidth', 1.5);
plot(T, x2, 'b-', 'LineWidth', 1.5);
plot(T, x3, 'c-', 'LineWidth', 1.5);
plot(T, x4, 'g-', 'LineWidth', 1.5);

plot(T, x5, 'r--', 'LineWidth', 1.5);
plot(T, x6, 'b--', 'LineWidth', 1.5);
plot(T, x7, 'c--', 'LineWidth', 1.5);
plot(T, x8, 'g--', 'LineWidth', 1.5);
xlabel('$k$', 'Interpreter', 'latex','FontSize', 26);
xlim([0,Tn]);
set(gca, 'FontSize', 26);
legend('$x_1$','$x_2$','$x_3$', '$x_4$', '$x_5$', '$x_6$', '$x_7$', '$x_8$', 'Interpreter', 'latex','FontSize', 26);
% title('Two areas power system state trajectories under jump-like FRP and DoS attacks');
hold off;