% 种群多样性分析

clc;clear;close all

SearchAgents=50; % 种群数量

Fun_name='F23';  % 基准测试函数F1到F23
Max_iterations=500; % 最大迭代次数
[lowerbound,upperbound,dimension,fitness]=fun_info(Fun_name); % 加载函数细节

% 函数调用
[Score,Best_pos,MENGO_curve,Ic]=GWO(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness);

% 绘图
plot(Ic,'Color','r','LineWidth',3)
title(Fun_name)
xlabel('Iteration#','fontsize',15);
ylabel('Population diversity','fontsize',15);
box on
axis tight
grid on
legend('GWO','FontSize',15)


