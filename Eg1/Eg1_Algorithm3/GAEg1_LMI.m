function [flag, gamma, R] = GAEg1_LMI(XXX)
XXX = XXX';

[A,B,C,F,G,kappa,T,mu,dd,c1,c2,Attack_signal,Channel_attacked,x00] = Eg1_paramaters();

omi_s3 = 0.35;
omi_s4 = 0.65;
omi_d3 = 0.35;
omi_d4 = 0.65;

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

barQ = diag([XXX(25),XXX(26)]);

R = sdpvar(9,9);
gamma = sdpvar(1);
P1 = sdpvar(9,9);
P2 = sdpvar(9,9);
P3 = sdpvar(9,9);
P4 = sdpvar(9,9);

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
 
barA1 = [A+B*Ec1*Phi_hat_1*C, B*Dc1, B*Ec1*Phi_bar_1;
           Bc1*Phi_hat_1*C, Ac1, Bc1*Phi_bar_1;
           Phi_hat_1*C, zeros(4,1), Phi_bar_1];
barA2 = [A+B*Ec2*Phi_hat_2*C, B*Dc2, B*Ec2*Phi_bar_2;
           Bc2*Phi_hat_2*C, Ac2, Bc2*Phi_bar_2;
           Phi_hat_2*C, zeros(4,1), Phi_bar_2];
barA3 = [A+B*Ec3*Phi_hat_3*C, B*Dc3, B*Ec3*Phi_bar_3;
           Bc3*Phi_hat_3*C, Ac3, Bc3*Phi_bar_3;
           Phi_hat_3*C, zeros(4,1), Phi_bar_3];
barA4 = [A+B*Ec4*Phi_hat_4*C, B*Dc4, B*Ec4*Phi_bar_4;
           Bc4*Phi_hat_4*C, Ac4, Bc4*Phi_bar_4;
           Phi_hat_4*C, zeros(4,1), Phi_bar_4];

% barC
barC = [G',zeros(1,1),zeros(1,4)]';

% barF i
barF = [F,zeros(1,1),zeros(1,4)];

% mathscr_C
M_C = bar_I2'*[C,zeros(4,1),-eye(4,4)];

% arr_Q
arr_Q_3 = barQ - barQ*Psi1;
arr_Q_4 = barQ - barQ*Psi2;

% Gamma_1
Gamma_1_1 = blkdiag(-(1+mu)*P1, -gamma/(1+mu)^T);
Gamma_1_2 = blkdiag(-(1+mu)*P2, -gamma/(1+mu)^T);

Gamma_1_3 = blkdiag(-(1+mu)*(P3 + M_C'*arr_Q_3*M_C), -gamma/(1+mu)^T);
Gamma_1_4 = blkdiag(-(1+mu)*(P4 + M_C'*arr_Q_4*M_C), -gamma/(1+mu)^T);

% Gamma_2 i j alpha
Gamma_2_1 = [barA1,barC];
Gamma_2_2 = [barA2,barC];
Gamma_2_3 = [barA3,barC];
Gamma_2_4 = [barA4,barC];

% Gamma_3 i
Gamma_3 = [barF,zeros(1,1)];

% Omega D S
Omega_D_3 = P3 + M_C'*arr_Q_3*M_C;
Omega_D_4 = P4 + M_C'*arr_Q_4*M_C;

Omega_S_1 = P1;
Omega_S_2 = P2;

%% 约束
% Xi
% Case 1
Lmat_1_33 = Gamma_1_3 + Gamma_3'*Gamma_3 + Gamma_2_3'*(P3 + omi_d3*M_C'*arr_Q_3*M_C + omi_d4*M_C'*arr_Q_4*M_C)*Gamma_2_3;
Lmat_1_34 = Gamma_1_3 + Gamma_3'*Gamma_3 + Gamma_2_3'*(P4 + omi_d3*M_C'*arr_Q_3*M_C + omi_d4*M_C'*arr_Q_4*M_C)*Gamma_2_3;

Lmat_1_43 = Gamma_1_4 + Gamma_3'*Gamma_3 + Gamma_2_4'*(P3 + omi_d3*M_C'*arr_Q_3*M_C + omi_d4*M_C'*arr_Q_4*M_C)*Gamma_2_4;
Lmat_1_44 = Gamma_1_4 + Gamma_3'*Gamma_3 + Gamma_2_4'*(P4 + omi_d3*M_C'*arr_Q_3*M_C + omi_d4*M_C'*arr_Q_4*M_C)*Gamma_2_4;

% Case 2
Lmat_2_31 = Gamma_1_3 + Gamma_3'*Gamma_3 + Gamma_2_3'*P1*Gamma_2_3;
Lmat_2_32 = Gamma_1_3 + Gamma_3'*Gamma_3 + Gamma_2_3'*P2*Gamma_2_3;

Lmat_2_41 = Gamma_1_4 + Gamma_3'*Gamma_3 + Gamma_2_4'*P1*Gamma_2_4;
Lmat_2_42 = Gamma_1_4 + Gamma_3'*Gamma_3 + Gamma_2_4'*P2*Gamma_2_4;

% Dynamic sagment FTB
Lmat_FTB_D1_3 = R-Omega_D_3;
Lmat_FTB_D1_4 = R-Omega_D_4;

Lmat_FTB_D2_3 = Omega_D_3 - kappa*R;
Lmat_FTB_D2_4 = Omega_D_4 - kappa*R;

% Case 3
Lmat_3_11 = Gamma_1_1 + Gamma_3'*Gamma_3 + Gamma_2_1'*P1*Gamma_2_1;
Lmat_3_12 = Gamma_1_1 + Gamma_3'*Gamma_3 + Gamma_2_1'*P2*Gamma_2_1;

Lmat_3_21 = Gamma_1_2 + Gamma_3'*Gamma_3 + Gamma_2_2'*P1*Gamma_2_2;
Lmat_3_22 = Gamma_1_2 + Gamma_3'*Gamma_3 + Gamma_2_2'*P2*Gamma_2_2;

% Case 4
Lmat_4_13 = Gamma_1_1 + Gamma_3'*Gamma_3 + Gamma_2_1'*(P3 + omi_s3*M_C'*arr_Q_3*M_C + omi_s4*M_C'*arr_Q_4*M_C)*Gamma_2_1;
Lmat_4_14 = Gamma_1_1 + Gamma_3'*Gamma_3 + Gamma_2_1'*(P4 + omi_s3*M_C'*arr_Q_3*M_C + omi_s4*M_C'*arr_Q_4*M_C)*Gamma_2_1;

Lmat_4_23 = Gamma_1_2 + Gamma_3'*Gamma_3 + Gamma_2_2'*(P3 + omi_s3*M_C'*arr_Q_3*M_C + omi_s4*M_C'*arr_Q_4*M_C)*Gamma_2_2;
Lmat_4_24 = Gamma_1_2 + Gamma_3'*Gamma_3 + Gamma_2_2'*(P4 + omi_s3*M_C'*arr_Q_3*M_C + omi_s4*M_C'*arr_Q_4*M_C)*Gamma_2_2;

% Static sagment FTB
Lmat_FTB_S1_1 = R-Omega_S_1;
Lmat_FTB_S1_2 = R-Omega_S_2;

Lmat_FTB_S2_1 = Omega_S_1 - kappa*R;
Lmat_FTB_S2_2 = Omega_S_2 - kappa*R;

L_FTB = (1+mu)^T*kappa*c1 + gamma*dd^2 - c2;

Lmat = [Lmat_1_33<=0;Lmat_1_34<=0;Lmat_1_43<=0;Lmat_1_44<=0;
        Lmat_2_31<=0;Lmat_2_32<=0;Lmat_2_41<=0;Lmat_2_42<=0;
        Lmat_3_11<=0;Lmat_3_12<=0;Lmat_3_21<=0;Lmat_3_22<=0;
        Lmat_4_13<=0;Lmat_4_14<=0;Lmat_4_23<=0;Lmat_4_24<=0;
        Lmat_FTB_D1_3<=0;Lmat_FTB_D1_4<=0;
        Lmat_FTB_D2_3<=0;Lmat_FTB_D2_4<=0;
        Lmat_FTB_S1_1<=0;Lmat_FTB_S2_1<=0;
        Lmat_FTB_S1_2<=0;Lmat_FTB_S2_2<=0;L_FTB<=0;
        P1>=0;P2>=0;P3>=0;P4>=0;barQ>=0;gamma>=0;R>=0];

ops = sdpsettings('verbose', 0, 'debug', 1);
sol = optimize(Lmat, gamma, ops);

%% 验证
for i = 1
    P1 = value(P1);
    P2 = value(P2);
    P3 = value(P3);
    P4 = value(P4);
    R = value(R);
    kappa = value(kappa);
    gamma = value(gamma);
    
    % FRP
    bar_I1 = [eye(2,2);zeros(2,2)];
    bar_I2 = [zeros(2,2);eye(2,2)];
    
    Phi1 = [1,0;0,0];
    Phi2 = [0,0;0,1];
    Psi1 = [1,0;0,0];
    Psi2 = [0,0;0,1];
    
    % static segment
    rho = 0;
    Phi_hat_1 = (1-rho)*bar_I1*Phi1*bar_I1';
    Phi_bar_1 = (1-rho)*bar_I1*(eye(2,2)-Phi1)*bar_I1';
    
    Phi_hat_2 = (1-rho)*bar_I1*Phi2*bar_I1';
    Phi_bar_2 = (1-rho)*bar_I1*(eye(2,2)-Phi2)*bar_I1';
    
    % dynamic segment
    rho = 1;
    Phi_hat_3 = rho*bar_I2*Psi1*bar_I2';
    Phi_bar_3 = rho*bar_I2*(eye(2,2)-Psi1)*bar_I2';
    
    Phi_hat_4 = rho*bar_I2*Psi2*bar_I2';
    Phi_bar_4 = rho*bar_I2*(eye(2,2)-Psi2)*bar_I2';
    
    % barA alpha
    barA1 = [A+B*Ec1*Phi_hat_1*C, B*Dc1, B*Ec1*Phi_bar_1;
               Bc1*Phi_hat_1*C, Ac1, Bc1*Phi_bar_1;
               Phi_hat_1*C, zeros(4,1), Phi_bar_1];
    barA2 = [A+B*Ec2*Phi_hat_2*C, B*Dc2, B*Ec2*Phi_bar_2;
               Bc2*Phi_hat_2*C, Ac2, Bc2*Phi_bar_2;
               Phi_hat_2*C, zeros(4,1), Phi_bar_2];
    barA3 = [A+B*Ec3*Phi_hat_3*C, B*Dc3, B*Ec3*Phi_bar_3;
               Bc3*Phi_hat_3*C, Ac3, Bc3*Phi_bar_3;
               Phi_hat_3*C, zeros(4,1), Phi_bar_3];
    barA4 = [A+B*Ec4*Phi_hat_4*C, B*Dc4, B*Ec4*Phi_bar_4;
               Bc4*Phi_hat_4*C, Ac4, Bc4*Phi_bar_4;
               Phi_hat_4*C, zeros(4,1), Phi_bar_4];
    
    % barC
    barC = [G',zeros(1,1),zeros(1,4)]';
    
    % barF i
    barF = [F,zeros(1,1),zeros(1,4)];
    
    % mathscr_C
    M_C = bar_I2'*[C,zeros(4,1),-eye(4,4)];
    
    % arr_Q
    arr_Q_3 = barQ - barQ*Psi1;
    arr_Q_4 = barQ - barQ*Psi2;
    
    % Gamma_1
    Gamma_1_1 = blkdiag(-(1+mu)*P1, -gamma/(1+mu)^T);
    Gamma_1_2 = blkdiag(-(1+mu)*P2, -gamma/(1+mu)^T);
    
    Gamma_1_3 = blkdiag(-(1+mu)*(P3 + M_C'*arr_Q_3*M_C), -gamma/(1+mu)^T);
    Gamma_1_4 = blkdiag(-(1+mu)*(P4 + M_C'*arr_Q_4*M_C), -gamma/(1+mu)^T);
    
    % Gamma_2 i j alpha
    Gamma_2_1 = [barA1,barC];
    Gamma_2_2 = [barA2,barC];
    Gamma_2_3 = [barA3,barC];
    Gamma_2_4 = [barA4,barC];
    
    % Gamma_3 i
    Gamma_3 = [barF,zeros(1,1)];
    
    % Omega D S
    Omega_D_3 = P3 + M_C'*arr_Q_3*M_C;
    Omega_D_4 = P4 + M_C'*arr_Q_4*M_C;
    
    Omega_S_1 = P1;
    Omega_S_2 = P2;
    
    %% 约束
    % Xi
    % Case 1
    Lmat_1_33 = Gamma_1_3 + Gamma_3'*Gamma_3 + Gamma_2_3'*(P3 + omi_d3*M_C'*arr_Q_3*M_C + omi_d4*M_C'*arr_Q_4*M_C)*Gamma_2_3;
    Lmat_1_34 = Gamma_1_3 + Gamma_3'*Gamma_3 + Gamma_2_3'*(P4 + omi_d3*M_C'*arr_Q_3*M_C + omi_d4*M_C'*arr_Q_4*M_C)*Gamma_2_3;
    
    Lmat_1_43 = Gamma_1_4 + Gamma_3'*Gamma_3 + Gamma_2_4'*(P3 + omi_d3*M_C'*arr_Q_3*M_C + omi_d4*M_C'*arr_Q_4*M_C)*Gamma_2_4;
    Lmat_1_44 = Gamma_1_4 + Gamma_3'*Gamma_3 + Gamma_2_4'*(P4 + omi_d3*M_C'*arr_Q_3*M_C + omi_d4*M_C'*arr_Q_4*M_C)*Gamma_2_4;
    
    % Case 2
    Lmat_2_31 = Gamma_1_3 + Gamma_3'*Gamma_3 + Gamma_2_3'*P1*Gamma_2_3;
    Lmat_2_32 = Gamma_1_3 + Gamma_3'*Gamma_3 + Gamma_2_3'*P2*Gamma_2_3;
    
    Lmat_2_41 = Gamma_1_4 + Gamma_3'*Gamma_3 + Gamma_2_4'*P1*Gamma_2_4;
    Lmat_2_42 = Gamma_1_4 + Gamma_3'*Gamma_3 + Gamma_2_4'*P2*Gamma_2_4;
    
    % Dynamic sagment FTB
    Lmat_FTB_D1_3 = R-Omega_D_3;
    Lmat_FTB_D1_4 = R-Omega_D_4;
    
    Lmat_FTB_D2_3 = Omega_D_3 - kappa*R;
    Lmat_FTB_D2_4 = Omega_D_4 - kappa*R;
    
    % Case 3
    Lmat_3_11 = Gamma_1_1 + Gamma_3'*Gamma_3 + Gamma_2_1'*P1*Gamma_2_1;
    Lmat_3_12 = Gamma_1_1 + Gamma_3'*Gamma_3 + Gamma_2_1'*P2*Gamma_2_1;
    
    Lmat_3_21 = Gamma_1_2 + Gamma_3'*Gamma_3 + Gamma_2_2'*P1*Gamma_2_2;
    Lmat_3_22 = Gamma_1_2 + Gamma_3'*Gamma_3 + Gamma_2_2'*P2*Gamma_2_2;
    
    % Case 4
    Lmat_4_13 = Gamma_1_1 + Gamma_3'*Gamma_3 + Gamma_2_1'*(P3 + omi_s3*M_C'*arr_Q_3*M_C + omi_s4*M_C'*arr_Q_4*M_C)*Gamma_2_1;
    Lmat_4_14 = Gamma_1_1 + Gamma_3'*Gamma_3 + Gamma_2_1'*(P4 + omi_s3*M_C'*arr_Q_3*M_C + omi_s4*M_C'*arr_Q_4*M_C)*Gamma_2_1;
    
    Lmat_4_23 = Gamma_1_2 + Gamma_3'*Gamma_3 + Gamma_2_2'*(P3 + omi_s3*M_C'*arr_Q_3*M_C + omi_s4*M_C'*arr_Q_4*M_C)*Gamma_2_2;
    Lmat_4_24 = Gamma_1_2 + Gamma_3'*Gamma_3 + Gamma_2_2'*(P4 + omi_s3*M_C'*arr_Q_3*M_C + omi_s4*M_C'*arr_Q_4*M_C)*Gamma_2_2;
    
    % Static sagment FTB
    Lmat_FTB_S1_1 = R-Omega_S_1;
    Lmat_FTB_S1_2 = R-Omega_S_2;
    
    Lmat_FTB_S2_1 = Omega_S_1 - kappa*R;
    Lmat_FTB_S2_2 = Omega_S_2 - kappa*R;
    
    L_FTB = (1+mu)^T*kappa*c1 + gamma*dd^2 - c2;
    
    flag = all(eig(L_FTB)<0)*all(eig(Lmat_FTB_D1_3))*all(eig(Lmat_FTB_D1_4))*all(eig(Lmat_FTB_D2_3))*all(eig(Lmat_FTB_D2_4))*all(eig(Lmat_FTB_S1_1)<0)*all(eig(Lmat_FTB_S1_2)<0)*all(eig(Lmat_FTB_S2_1)<0)*all(eig(Lmat_FTB_S2_2)<0);
    flag = flag * all(eig(Lmat_1_33)<0)*all(eig(Lmat_1_34)<0)*all(eig(Lmat_1_43)<0)*all(eig(Lmat_1_44)<0);
    flag = flag * all(eig(Lmat_2_31)<0)*all(eig(Lmat_2_32)<0)*all(eig(Lmat_2_41)<0)*all(eig(Lmat_2_42)<0);
    flag = flag * all(eig(Lmat_3_11)<0)*all(eig(Lmat_3_12)<0)*all(eig(Lmat_3_21)<0)*all(eig(Lmat_3_22)<0);
    flag = flag * all(eig(Lmat_4_13)<0)*all(eig(Lmat_4_14)<0)*all(eig(Lmat_4_23)<0)*all(eig(Lmat_4_24)<0);
    flag = flag * all(eig(P1)>0)*all(eig(P2)>0)*all(eig(P3)>0)*all(eig(P4)>0)*all(eig(barQ)>0)*all(eig(R)>0);
end