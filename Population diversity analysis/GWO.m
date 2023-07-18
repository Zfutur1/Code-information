function [Alpha_score, Alpha_pos, Convergence_curve,Ic] = GWO(N, Max_iter, lb, ub, dim, fobj)
%% 灰狼优化算法
Positions = initialization(N, dim, ub, lb);
% 初始化alpha, beta和delta_pos
Alpha_pos = zeros(1, dim);
Alpha_score = inf;  % 将此更改为-inf以解决最大化问题
Beta_pos = zeros(1,dim);
Beta_score = inf;    % 将此更改为-inf以解决最大化问题
Delta_pos = zeros(1,dim);
Delta_score = inf;   % 将此更改为-inf以解决最大化问题
% 初始化种群位置
Convergence_curve = zeros(1, Max_iter);  % 收敛曲线
l = 0;    % 循环计数器
% 主要循环
while l < Max_iter
    for i = 1:N
        % 边界处理
        Flag4ub = Positions(i, :)>ub;
        Flag4lb = Positions(i, :)<lb;
        Positions(i, :) = (Positions(i, :).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
%         Positions(i, :) = max(Positions(i, :), lb);
%         Positions(i, :) = min(Positions(i, :), ub);
        % 计算每个种群的目标函数
        fitness = fobj(Positions(i, :));
        % 更新Alpha, Beta和Delta
        if fitness < Alpha_score
            Alpha_score = fitness;        % 更新alpha
            Alpha_pos = Positions(i, :);
        end
        if fitness > Alpha_score && fitness < Beta_score
            Beta_score = fitness;           % 更新beta
            Beta_pos = Positions(i, :);
        end
        if fitness > Alpha_score && fitness > Beta_score && fitness < Delta_score
            Delta_score = fitness;          % 更新delta
            Delta_pos = Positions(i, :);
        end
    end
    a = 2-l*((2)/Max_iter);      % a从2线性减小到0    
    % 更新包括omegas在内的种群的位置
    for i = 1:N
        for j = 1:dim
            r1 = rand();             % r1是[0,1]中的随机数
            r2 = rand();             % r2是[0,1]中的随机数
            
            A1 = 2*a*r1-a;        % 公式（4）
            C1 = 2*r2;               % 公式（5）
            
            D_alpha = abs(C1*Alpha_pos(j)-Positions(i, j));  % 公式（6）-第一部分
            X1 = Alpha_pos(j)-A1*D_alpha;   % 公式 (7)-第一部分
            
            r1 = rand();
            r2 = rand();
            
            A2 = 2*a*r1-a;         % 公式（4）
            C2 = 2*r2;                % 公式（5）
            
            D_beta = abs(C2*Beta_pos(j)-Positions(i, j));   % 公式（6）-第二部分
            X2 = Beta_pos(j)-A2*D_beta;       % 公式 (7)-第二部分
            
            r1 = rand();
            r2 = rand();
            
            A3 = 2*a*r1-a;        % 公式 (4)
            C3 = 2*r2;               % 公式 (5)
            
            D_delta = abs(C3*Delta_pos(j)-Positions(i, j)); % 公式（6）-第三部分
            X3 = Delta_pos(j)-A3*D_delta;      % 公式 (7)-第三部分
            
            Positions(i, j)=(X1+X2+X3)/3;       % 公式 (8)
        end
    end
    l = l+1;
    Cd=(1/dim)*sum(Positions);
    M=(Positions-Cd).^2;
    Ic(l)=(sum(sum(M)))^0.5;
    Convergence_curve(l) = Alpha_score;
    disp(['GWO: At iteration ', num2str(l), ' ,the best fitness is ', num2str(Alpha_score)]);
end



