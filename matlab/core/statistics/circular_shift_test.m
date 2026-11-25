function [r_obs, p_val] = circular_shift_test(metric_ts, task_design, n_perms)
% ANALYZE_TASK_CORRELATION
% 输入:
%   metric_ts   : 你计算得到的指标时间序列 (N x 1)
%   task_design : 原始的 0/1 向量 (N x 1), 1=Task, 0=Fixation
%
% 输出:
%   r_obs         : 观测到的相关系数 (Pearson r)
%   p_val         : 基于循环移位的 P 值


% 计算原始相关性
r_obs = corr(metric_ts, task_design);

% 循环移位检验 (Circular Shift Test)
if ~exist('n_perms','var')
    n_perms = 5000;
end
null_dist = zeros(n_perms, 1);
N = length(task_design);

for i = 1:n_perms
    % 生成一个随机的时间偏移量 (至少偏移 10 个 TR 以打破局部依赖)
    shift_amount = randi([10, N-10]);

    % 对任务回归因子进行循环移位
    reg_shifted = circshift(task_design, shift_amount);

    % 计算零分布的相关性
    null_dist(i) = corr(metric_ts, reg_shifted);
end

% 计算 P 值 (双尾)
p_val = (sum(abs(null_dist) >= abs(r_obs)) + 1) / (n_perms + 1);
