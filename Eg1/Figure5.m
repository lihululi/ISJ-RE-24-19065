clc
clear
close all

T = (0:28);
Tn = 28;

Channel_attacked = [0	1	4	4	4	4	4	4	0	1	1	4	4	4	4	0	0	1	1	4	4	3	3	3	3	0	0	0	0];
Channel_selected = [1	1	2	3	3	1	2	4	4	1	2	4	4	1	2	4	3	1	2	4	4	1	2	4	3	1	2	4	3];
Channel_selected_jump = [1	2	1	3	3	2	1	3	4	2	2	3	3	1	2	4	3	2	2	3	3	1	2	4	4	1	2	4	3];

figure(1)
subplot(2,1,2)
hold on;
box on;
stairs(T,Channel_selected,'LineWidth',1.8); 
stairs(T,Channel_attacked, '--','LineWidth',1.8);
legend('The selected Sensor node', 'The attacked Sensor node', 'Interpreter', 'latex','FontSize', 22, 'Location', 'Northeast');
ylim([-0.1 6]);
xlim([-0.1 Tn]);
xlabel('$k$', 'Interpreter', 'latex','FontSize', 22);
ylabel('$\alpha(k)$', 'Interpreter', 'latex','FontSize', 22);
xticks(0:4:28);
set(gca, 'FontSize', 22);
title('The selected sensor node and the attacked sensor node under traditional FRP','FontSize', 22);

subplot(2,1,1)
hold on;
box on;
stairs(T,Channel_selected_jump,'LineWidth',1.8); 
stairs(T,Channel_attacked, '--','LineWidth',1.8);
legend('The selected Sensor node', 'The attacked Sensor node', 'Interpreter', 'latex','FontSize', 22, 'Location', 'Northeast');
ylim([-0.1 6]);
xlim([-0.1 Tn]);
xlabel('$k$', 'Interpreter', 'latex','FontSize', 22);
ylabel('$\alpha(k)$', 'Interpreter', 'latex','FontSize', 22);
xticks(0:4:28);
set(gca, 'FontSize', 22);
title('The selected sensor node and the attacked sensor node under jump-like FRP','FontSize', 22);
