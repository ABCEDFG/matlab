%   使用传统的 GM(1, 1) 模型对数据进行预测
%   输入：
%       data：数据
%       forecast_num：预测的期数
%   输出：
%       r：预测值
%       data_hat：拟合值
%       relative_residuals：相对残差
%       data_eta：级比偏差
function [r, data_hat, relative_residuals, data_eta] = gm11(data, forecast_num)

data_len = length(data);  % 数据的长度
data_ago = cumsum(data);  % 计算累加值
data_z = (data_ago(1 : end - 1) + data_ago(2 : end)) / 2;  % 计算紧邻均值生成数列（长度为n-1） 


y = data(2 : end); 
x = data_z;

k = ((data_len - 1) * sum(x .* y) - sum(x) * sum(y)) / ((data_len - 1) * sum(x .* x) - sum(x) * sum(x));  % 发展系数
b = (sum(x .* x) * sum(y) - sum(x) * sum(x .* y)) / ((data_len - 1) * sum(x .* x) - sum(x) * sum(x));  % b 就是灰作用量
a = -k;  

fprintf("发展系数为：%f\n", k);
fprintf("灰作用量为：%f\n\n", b);

data_hat = zeros(data_len, 1);  % 用于存放拟合值
data_hat(1) = data(1);

for i = 1 : data_len - 1
    data_hat(i + 1) = (1 - exp(a)) * (data(1) - b / a) * exp(-a * i);  
end

r = zeros(forecast_num, 1);
for i = 1 : forecast_num
    r(i) = (1 - exp(a)) * (data(1) - b / a) * exp(-a * (data_len + i - 1)); 
end

% 计算绝对残差和相对残差
absolute_residuals = data(2 : end) - data_hat(2 : end);  % 绝对残差
relative_residuals = abs(absolute_residuals) ./ data(2 : end);  % 相对残差


% 计算级比和级比偏差
class_ratio = data(2 : end) ./ data(1 : end - 1) ;  % 计算级比 sigma(k) = x0(k)/x0(k-1)
data_eta = abs(1 - (1 - 0.5 * a) / (1 + 0.5 * a) * (1 ./ class_ratio));  % 计算级比偏差

