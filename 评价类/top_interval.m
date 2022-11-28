% TOPSIS 中间型转极大型
function R = top_interval(A, minx, manx)

    A_r = size(A, 1);  % 求向量长度
    R = zeros(A_r, 1);  % 生成一个和 A 相同长度的 0 向量

    A_mm = max([minx - min(A), max(A) - manx]);  % 求基准值 M

    % 分类计算，化为及大型
    for i = 1 : A_r
        if A(i) < minx
            R(i) = 1 - (minx - A(i)) / A_mm;
        elseif A(i) > manx
            R(i) = 1 - (A(i) - manx) / A_mm;
        else
            R(i) = 1;
        end
    end
end