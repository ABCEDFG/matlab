% 计算斯皮尔曼相关系数（数据维数尽量大于 30）
% 输入：
%   A, B：输入的数值序列
% 输出：
%   r：A 和 B 的斯皮尔曼相关系数
%   t：t 检验，需查t检验（t-test）临界值表-t检验表以验证皮尔逊相关系数的关系，大于表中数据通过
function [r, t] = spearman(A, B)

% 判断两个数据的维数是否相同
A_len = length(A);
B_len = length(B);
if A_len ~= B_len
    error('两数据的维数不相等');
    return;
end

% 计算 A_rank 排列序号
A_rank = zeros(1, A_len);
for i = 1 : A_len
    cont_big = 0;  % 记录大于特定元素的个数
    cont_equal = -1;  % 记录等于特定元素的个数
    
    for j = 1 : A_len
        if A(i) < A(j)
            cont_big = cont_big + 1;
        elseif A(i) == A(j)
            cont_equal = cont_equal + 1;
        end
    end
    A_rank(i) = cont_big + mean([0 : cont_equal]);
end

% 计算 B_rank 排列序号
B_rank = zeros(1, B_len);
for i = 1 : B_len
    cont_big = 0;  % 记录大于特定元素的个数
    cont_equal = -1;  % 记录等于特定元素的个数 
    
    for j = 1 : B_len
        if B(i) < B(j)
            cont_big = cont_big + 1;
        elseif B(i) == B(j)
            cont_equal = cont_equal + 1;
        end
    end
    B_rank(i) = cont_big + mean([0 : cont_equal]);
end

% 计算斯皮尔曼等级相关系数
len = A_len;
r_fenzi = 6 * sum((A_rank - B_rank) .^ 2);  % 计算斯皮尔曼等级相关系数的分子
r_fenmu = len * (len ^ 2 - 1);  % 计算斯皮尔曼等级相关系数的分母
r = 1 - r_fenzi / r_fenmu;  % 计算斯皮尔曼等级相关系数

t = r * (((len - 2) / (1 - r ^ 2)) ^ 0.5);  % t 统计检验
