%----------求解TSP问题----------%
clc;
clear all;
close all;
%给定城市坐标City_Coord，每一列代表一个城市
City_Coord = [0.4000, 0.2439, 0.1707, 0.2293, 0.5171, 0.8732, 0.6878, 0.8488, 0.6683, 0.6195;
              0.4439, 0.1463, 0.2293, 0.7610, 0.9414, 0.6536, 0.5219, 0.3609, 0.2536, 0.2634];
%获取城市数量city_quantity
n = size(City_Coord, 2);     %City_Coord的列数即为城市数量
%计算城市之间的距离city_distance矩阵    %city_quantity行city_quantity列，第i行第j列代表城市i到城市j的距离
for i = 1:n
    for j = 1:n
        if i <= j
            City_Distance(i,j) = norm(City_Coord(:,i)-City_Coord(:,j));
            City_Distance(j,i) = City_Distance(i,j);
        end
    end
end

%----------设置参数----------%
m = 20;    %种群个体数目
pc = 0.5;             %交叉概率
pm = 0.1;           %变异概率
n_loop_max = 100;     %最大迭代次数

%----------一些矩阵说明----------%
History_XBest = [];                   %存储每次循环中最优的路径，每一列为依次循环中的路径
History_X_FunctionBest = [];          %存储每次循环中最优的路径长度，行向量

%----------主循环----------%
for n_loop = 1:n_loop_max
    if n_loop == 1
        X = X_init(n, m);           %种群初始化
    end
    
    %将当前种群最优路径存储进历史记录
    X_Function = GetFunction(n, m, City_Distance, X);
    [min_X_Function, min_X_Function_position] = min(X_Function);
    History_XBest(:, n_loop) = X(:, min_X_Function_position);         %存储当前最优路径
    History_X_FunctionBest(n_loop) = min_X_Function;                  %存储当前最优路径的长度
   
    fprintf('第%d次迭代\n',n_loop)
    fprintf('初始种群：')
    display(X)

    %洗牌（打乱顺序）
    shuffle = randperm(m);
    X = X(:,shuffle);
    fprintf('洗牌后：')
    display(X)
        
    %交叉
    X = Cross(m, n, X, pc);
    fprintf('交叉后：')
    display(X)
    
    %变异
    X = mutate(n, m, X, pm);
    fprintf('变异后：')
    display(X)
    
    %选择
    X = Select(n, m, City_Distance, X);
    fprintf('选择后的新种群：')
    display(X)
end

%将最后一次迭代种群最优路径存储进历史记录
X_Function = GetFunction(n, m, City_Distance, X);
[min_X_Function, min_X_Function_position] = min(X_Function);
History_XBest(:, n_loop) = X(:, min_X_Function_position);         %存储当前最优路径
History_X_FunctionBest(n_loop) = min_X_Function;                  %存储当前最优路径的长度

Xbest = X(:, min_X_Function_position);
XFbest = min_X_Function;

%**********************打印结果**********************%
fprintf('二进制编码混合遗传算法：\n')
fprintf('****************参数设置****************\n')
fprintf('种群个体数目: %d\n', m)
fprintf('交叉概率: %f\n', pc)
fprintf('变异概率: %f\n', pm)
fprintf('最大迭代次数: %d\n', n_loop_max)
fprintf('****************迭代结果****************\n')
fprintf('最优路径为： %d %d %d %d %d %d %d %d %d %d\n', Xbest)
fprintf('最优路径的长度为： %f\n', XFbest)
fprintf('****************************************\n')
fprintf('end\n')

                                                                           
%**********************************************************************%
%****************************种群X初始化函数****************************%
%**********************************************************************%
function X = X_init(n, m)
    for i = 1:m                          %随机生成城市顺序
        X(:, i) = randperm(n);
    end
end


%**********************************************************************%
%*********************计算函数值（路径总长度）函数**********************%
%**********************************************************************%
function X_Function = GetFunction(n, m, City_Distance, X)
    for i = 1:m
        X_Function(i) = 0;
        for ii = 1:(n - 1)
            X_Function(i) = X_Function(i) + City_Distance(X(ii, i), X(ii+1, i));
        end
        X_Function(i) = X_Function(i) + City_Distance(X(n, i), X(1, i));
    end
end

%**********************************************************************%
%****************************计算适应度函数****************************%
%**********************************************************************%
function X_Fitness = GetFitness(n, m, City_Distance, X)
    X_Function = GetFunction(n, m, City_Distance, X);
    for i = 1:m
        X_Fitness(i) = 1000 / X_Function(i);
    end
end

%**********************************************************************%
%*********************相邻两个染色体两点交叉算法函数**********************%
%**********************************************************************%
function X = Cross(m, n, X, pc)
    if mod(m ,2) == 0
        cross_num = m;        %如果个体数为偶数，参与交叉个体数为全部个体数
    else
        cross_num = m - 1;    %如果个体数为奇数，参与交叉个体数为全部个体数-1，,最后一个不参与交叉
    end
    for i = 1:2:cross_num
        hand_of_god = rand;            %上帝之手随机做选择
        if hand_of_god < pc            %如果上帝手小于交叉概率pc，就进行交叉，否则不交叉
            %A、B为相邻两条染色体，减少索引次数
            A = X(:, i);               
            B = X(:, i+1);
            %确定交叉范围
            first = round(rand * (n - 1)) + 1;              %交叉段起点
            last =  round(rand * (n - first)) + first;      %交叉段结束点
            %交叉
            Cross_Temp = A(first:last);         %Cross_Temp为临时存储
            A(first:last) = B(first:last);
            B(first:last) = Cross_Temp;
            %检查并修正冲突
            A = Reslove_Conflict(n, first, last, A, B);
            B = Reslove_Conflict(n, first, last, B, A);
            X(:, i) = A;               
            X(:, i+1) = B;
        end
    end 
end
%**********************************************************************%
%***********************检查并修正冲突(交叉结束后)************************%
%**********************************************************************%
function A = Reslove_Conflict(n, first, last, A, B)
    %寻找A中重复出现的值（冲突）
    while 1
        k = 0;
        for ii = 1:n
            for iii = ii+1 : n
                if A(ii) == A(iii)
                    k = k+1;
                    Save(1, k) = A(ii);
                    Save(2, k) = ii;     %Save为3行重复值的个数列的矩阵，第一行是重复值的数值，第二三行是重复值的位置
                    Save(3, k) = iii;
                end
            end
        end
        if k == 0
            break
        end
        %替换A中的重复值
        for ii = 1:k
            if Save(2, ii) >= first && Save(2, ii) <= last        %找到在交叉范围内的重复值的位置，使其为position1
                position1 = Save(2, ii);
                position2 = Save(3, ii);
            else
                position1 = Save(3, ii);
                position2 = Save(2, ii);
            end
            A(position2) = B(position1);
        end
    end
end

%**********************************************************************%
%****************************变异函数****************************%
%**********************************************************************%
 function X = mutate(n, m, X, pm)
    for i = 1:m
        hand_of_god = rand;            %上帝之手随机做选择
        if hand_of_god < pm            %如果上帝手小于变异概率pc，就进行变异，否则不变异
            Temp = X(:, i);
            %确定变异范围
            first = round(rand * (n - 1)) + 1;              %变异段起点
            last =  round(rand * (n - first)) + first;      %变异段结束点
            if last == first
                if last ~= n
                    last = last + 1;
                else
                    first ~= 1;
                    first = first - 1;
                end
            end
            Temp(first:last) = flipud(Temp(first:last));          %逆序（变异）
            X(:, i) = Temp;
        end
    end
 end

%**********************************************************************%
%********************************选择函数********************************%
%**********************************************************************%
function X = Select(n, m, City_Distance, X)
    X_Fitness = GetFitness(n, m, City_Distance, X);
    fitness_sum = sum(X_Fitness);                          %适应度求和
    Possibility = cumsum(X_Fitness / fitness_sum);         %根据适应度计算，X_Fitness / fitness_Sum是被选择到的概率，Possibility为累积概率
    for i = 1:(m - 1)                           %n_Population选n_Population-1，可重复选
        hand_of_god = rand;                                %“上帝之手”转动轮盘，皮一下好开心
        for ii = 1:m
            if hand_of_god < Possibility(ii)
                Temp_X(:,i) = X(:,ii);                     %Temp_X为临时存储选择到的值的矩阵
                break
            end
        end
    end
    [~, max_x_position] = max(X_Fitness);                 %找出最大适应度索引
    Temp_X(:, m) = X(:, max_x_position);
    X = Temp_X;                                           %得到选择后的种群
end