clc
clear
close all

[A,B,C,F,G,kappa,T,mu,dd,c1,c2,Attack_signal,Channel_attacked,x00,R] = Eg1_paramaters();

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

%% closed loop
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

disturbance = [];
for k = 1:T+1
    disturbance = [disturbance,1.5*(sin(k)^2)/(1+k)];
end

for k = 0:T-1
    % traditional FRP
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

%% Compare
Tn = T;
T = (0:T);
set(groot, 'defaultAxesColorOrder', lines);

State = C_State;

% obtained by Eg1_A2_simulation.m
C_State = [0.0244500000000000	0.0306834701836807	0.0359294449641372	0.0267840270194724	0.0197101874180875	0.0266899296294753	0.0305083414929824	0.0244241914952572	0.0195658612484347	0.0232307283278252	0.0257722773954880	0.0211204710115724	0.0159701327366221	0.0175745201612985	0.0192293521783411	0.0158940618325856	0.0126390270778255	0.0136620181540592	0.0147109863759724	0.0115749560061394	0.00994781628957806	0.0103771558126032	0.0110432916924506	0.00996380950132100	0.00892916590564945	0.00919039103435894	0.00958596788555470	0.00793332681264590	0.00752223475709261
0.00815000000000000	0.0188357276377606	0.0332942913331743	0.0303894091816855	0.0195685770282610	0.0235923711500462	0.0259984691308824	0.0159512680520015	0.00999063044492301	0.0143872493731837	0.0145039935971487	0.00871945995052641	0.00579130870998991	0.00740337731829946	0.00805800264049517	0.00611044077693294	0.00370569226509659	0.00432955817852070	0.00609661256593059	0.00368282044954953	0.00193039592747749	0.00325046659231561	0.00435184903222602	0.00270000148508023	0.00214689576095219	0.00298145889331910	0.00271906563200826	0.00136902691022943	0.00173816116577336
0.00815000000000000	0.0625687340082784	0.0757852912199998	0.0288739280627454	0.0234482490408425	0.0383758212828677	0.0216262903806845	0.0118612825573874	0.0173438679485439	0.0134042353688876	0.0114622856916149	0.0135343361514925	0.00569761640603896	0.00568162717942813	0.0134023345624035	0.00757767359438190	0.00166091800983084	0.00913085241258117	0.00954888100537605	0.00207465020080101	0.00550726416921290	0.00778902865625865	0.00413841281598508	0.00531413016836663	0.00657952556392610	0.00340835108108187	0.00492360972454130	0.00584882382485980	0.00286985320869070
0.154850000000000	0.164056940367361	0.167014598873170	0.135406046450571	0.109102603575127	0.115234657186598	0.113434532453594	0.0908668448057905	0.0747696770913449	0.0781031590721566	0.0774313574582088	0.0623261472515182	0.0485376867214070	0.0494299380867637	0.0497782511935293	0.0398844694458201	0.0309302236255268	0.0314939288719556	0.0318600420443737	0.0231493362783355	0.0184219286653930	0.0185084493404522	0.0185405344632208	0.0146408304680712	0.0114915426826960	0.0111552346357341	0.0107538047807980	0.00636089638217709	0.00499110150697874];
cx1 = C_State(1,:);
cx2 = C_State(2,:);
cx3 = C_State(3,:);
cx4 = C_State(4,:);
x1 = State(1,:);
x2 = State(2,:);
x3 = State(3,:);
x4 = State(4,:);

set(groot, 'defaultAxesColorOrder', lines);
figure(1)
subplot(4,1,1)
box on;
hold on;
plot(T, x1, 'r--', 'LineWidth', 2);
plot(T, cx1, 'r-', 'LineWidth', 2);
legend('$x_1$ of FRP with Dos attacks', '$x_1$ of the jump-like FRP with Dos attacks', 'Interpreter', 'latex','FontSize', 18, 'Location', 'Northeast');
xlim([0,28]);
ylim([-0.02,0.23]);
xticks(0:4:28);
xlabel('$k$', 'Interpreter', 'latex','FontSize', 24);
set(gca, 'FontSize', 24);
hold off;

subplot(4,1,2)
box on;
hold on;
plot(T, x2, 'b--', 'LineWidth', 2);
plot(T, cx2, 'b-', 'LineWidth', 2);
legend('$x_2$ of FRP with Dos attacks', '$x_2$ of the jump-like FRP with Dos attacks', 'Interpreter', 'latex','FontSize', 18, 'Location', 'Northeast');
xlim([0,28]);
ylim([0,0.05]);
xticks(0:4:28);
xlabel('$k$', 'Interpreter', 'latex','FontSize', 24);
set(gca, 'FontSize', 24);
hold off;

subplot(4,1,3)
box on;
hold on;
plot(T, x3, 'c--', 'LineWidth', 2);
plot(T, cx3, 'c-', 'LineWidth', 2);
legend('$x_3$ of FRP with Dos attacks', '$x_3$ of the jump-like FRP with Dos attacks', 'Interpreter', 'latex','FontSize', 18, 'Location', 'Northeast');
xlim([0,28]);
ylim([-0.02,0.08]);
xticks(0:4:28);
yticks([0,0.05,0.08]);
xlabel('$k$', 'Interpreter', 'latex','FontSize', 24);
set(gca, 'FontSize', 24);
hold off;

subplot(4,1,4)
box on;
hold on;
plot(T, x4, 'm--', 'LineWidth', 2);
plot(T, cx4, 'm-', 'LineWidth', 2);
legend('$x_4$ of FRP with Dos attacks', '$x_4$ of the jump-like FRP with Dos attacks', 'Interpreter', 'latex','FontSize', 18, 'Location', 'Northeast');
xlim([0,28]);
ylim([0,0.25]);
xticks(0:4:28);
xlabel('$k$', 'Interpreter', 'latex','FontSize', 24);
set(gca, 'FontSize', 24);
hold off;
