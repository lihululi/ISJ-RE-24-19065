clc
clear
close all

a = -0.2:0.05:0.2;
b = -0.2:0.05:0.2;

LMIA = [2.0481	2.0865	2.1357	2.1960	2.2690	2.3581	2.4752	2.6644	3.2224];
GA_A = [0.1086 0.1126 0.1170 0.1219	0.1283 0.1381 0.1577 0.1901 0.2976];

LMIB = [1.3585	1.5677	1.7893	2.0230	2.2690	2.5273	2.7982	3.0821	3.3796];
GA_B = [0.0712 0.0837  0.0974 0.1122 0.1283 0.1451 0.1633 0.1829 0.2042];

LMIA = sqrt(LMIA);
GA_A = sqrt(GA_A);

LMIB = sqrt(LMIB);
GA_B = sqrt(GA_B);

figure(1)
subplot(2,1,1)
box on
hold on
plot(a, LMIA, '--+', Linewidth = 2.5, Markersize = 20);
plot(a, GA_A, '--x', Linewidth = 2.5, Markersize = 20);

axis([-0.2, 0.2, 0, 3]);
xlabel('$a$', 'Interpreter', 'latex','FontSize', 28);
ylabel('$\gamma_{\mathrm{min}}$', 'Interpreter', 'latex','FontSize', 28);
legend('The value of $\gamma_{\mathrm{min}}$ obtained by LMI approach (Theorem 2)', 'The value of $\gamma_{\mathrm{min}}$ obtained by solving Algorithm 3', 'interpreter', 'latex', 'FontSize', 28, 'Location', 'Northwest');
set(gca, 'FontSize', 28);
title('$\gamma_{\mathrm{min}}$ with respect to the parameter $a$ with $b=0$', 'Interpreter', 'latex','FontSize', 28);
hold off

figure(1)
subplot(2,1,2)
box on
hold on
plot(b, LMIB, '--+', Linewidth = 2.5, Markersize = 20);
plot(b, GA_B, '--x', Linewidth = 2.5, Markersize = 20);

axis([-0.2, 0.2, 0, 3]);
xlabel('$b$', 'Interpreter', 'latex','FontSize', 28);
ylabel('$\gamma_{\mathrm{min}}$', 'Interpreter', 'latex','FontSize', 28);
legend('The value of $\gamma_{\mathrm{min}}$ obtained by LMI approach (Theorem 2)', 'The value of $\gamma_{\mathrm{min}}$ obtained by solving Algorithm 3', 'interpreter', 'latex', 'FontSize', 28, 'Location', 'Northwest');
set(gca, 'FontSize', 28);
title('$\gamma_{\mathrm{min}}$ with respect to the parameter $b$ with $a=0$', 'Interpreter', 'latex','FontSize', 28);
hold off
