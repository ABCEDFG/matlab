%   新信息 GM（1，1）模型
%   输入：
%       data：数据
%       forecast_num：预测的期数
%   输出：
%       r：预测值
function [r] = gm11_new(data, forecast_num)

r = zeros(forecast_num, 1);

for i = 1 : forecast_num
    r(i) = gm11(data, 1);
    data = [data; r(i)];
end
