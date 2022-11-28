% 计算皮尔逊相关系数（数据维数尽量大于 30）
% 输入：
%   A, B：输入的数值序列
% 输出：
%   r：A 和 B 的皮尔逊相关系数
%   t：t 检验，需查t检验（t-test）临界值表-t检验表以验证皮尔逊相关系数的关系，大于表中数据通过
function [r, t] = pearson(A, B)

% 判断两个数据的维数是否相同
A_len = length(A);
if A_len ~= length(B)
    error('两数据的维数不相等');
    return;
end

% 计算皮尔逊相关系数
rr_fenzi = sum((A - sum(A) / length(A)) .* (B - sum(B) / length(B)));  % 计算皮尔逊相关系数的分子
rr_fenmu = sqrt(sum((A - sum(A) / length(A)) .^ 2)) * sqrt(sum((B - sum(B) / length(B)) .^ 2));  % 计算皮尔逊相关系数的分母
rr = rr_fenzi / rr_fenmu;  % 计算皮尔逊相关系数

% 修正
if A_len > 30
    r = rr;
elseif A_len > 4 && A_len < 30
    r = rr * (1 + (1 - rr ^ 2) / (2 * (A_len - 4)));  % 计算无偏相关系数加以矫正
else
    r = rr;
    fprintf('数据长度小于5 %8.5\n\n', r); 
end


t = abs(r / (sqrt((1 - r ^ 2) / A_len - 2)));  % t 检验

