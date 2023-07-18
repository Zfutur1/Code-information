% ��Ⱥ�����Է���

clc;clear;close all

SearchAgents=50; % ��Ⱥ����

Fun_name='F23';  % ��׼���Ժ���F1��F23
Max_iterations=500; % ����������
[lowerbound,upperbound,dimension,fitness]=fun_info(Fun_name); % ���غ���ϸ��

% ��������
[Score,Best_pos,MENGO_curve,Ic]=GWO(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness);

% ��ͼ
plot(Ic,'Color','r','LineWidth',3)
title(Fun_name)
xlabel('Iteration#','fontsize',15);
ylabel('Population diversity','fontsize',15);
box on
axis tight
grid on
legend('GWO','FontSize',15)


