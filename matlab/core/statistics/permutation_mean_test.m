function [p_val, obs_t, perm_dist] = permutation_mean_test(x, y, n_perms, tail)
% PERM_TEST_IND 独立样本的置换检验 (Permutation Test for Independent Samples)
%
% 输入:
%   x, y    : 两个维度的输入向量 (可以 nx != ny)
%   n_perms : 置换次数 (建议 >= 5000)
%   tail    : 'both' (双尾, 默认), 'right' (x>y), 'left' (x<y)
%
% 输出:
%   p_val     : 计算得到的 P 值
%   obs_t     : 原始数据的 t 统计量
%   perm_dist : 置换产生的零分布 (n_perms x 1)
%
% 参考文献: 
%   Good, P. I. (2005). Permutation, Parametric and Bootstrap Tests of Hypotheses.

    if nargin < 3, n_perms = 10000; end
    if nargin < 4, tail = 'both'; end

    % 确保列向量
    x = x(:); 
    y = y(:);
    
    nx = length(x);
    ny = length(y);
    n_total = nx + ny;
    
    % 合并数据
    combined_data = [x; y];
    
    % 1. 计算原始观测的 t 统计量
    obs_t = get_t_stat(x, y);
    
    % 2. 置换循环 (预分配内存加速)
    perm_dist = zeros(n_perms, 1);
    
    % 并行计算提示: 如果 n_perms 很大或计算体素很多，可以将 for 改为 parfor
    for i = 1:n_perms
        % 随机打乱索引
        perm_idx = randperm(n_total);
        
        % 根据原始样本量切分数据
        x_perm = combined_data(perm_idx(1:nx));
        y_perm = combined_data(perm_idx(nx+1:end));
        
        % 计算置换后的 t 统计量
        perm_dist(i) = get_t_stat(x_perm, y_perm);
    end
    
    % 3. 计算 P 值 (注意 +1 修正，防止 p=0)
    switch tail
        case 'both'
            % 双尾: 看绝对值是否更极端
            p_val = (sum(abs(perm_dist) >= abs(obs_t)) + 1) / (n_perms + 1);
        case 'right'
            % 右尾: H1: x > y, t 越大越显著
            p_val = (sum(perm_dist >= obs_t) + 1) / (n_perms + 1);
        case 'left'
            % 左尾: H1: x < y, t 越小越显著
            p_val = (sum(perm_dist <= obs_t) + 1) / (n_perms + 1);
    end
    
    % 可视化 (可选)
    % figure;
    % histogram(perm_dist, 50, 'Normalization', 'pdf', 'FaceColor', [0.7 0.7 0.7]);
    % hold on;
    % xline(obs_t, 'r', 'LineWidth', 2);
    % title(['Permutation Test Distribution (p = ' num2str(p_val, '%.4f') ')']);
    % xlabel('t-statistic'); ylabel('Density');
    % legend('Null Distribution', 'Observed t');
    % box off; set(gca, 'FontSize', 12);
    
end

function t = get_t_stat(g1, g2)
    % 辅助函数：计算双样本 t 统计量 (不等方差假设, Welch's t-test 形式)
    m1 = mean(g1); m2 = mean(g2);
    v1 = var(g1);  v2 = var(g2);
    n1 = length(g1); n2 = length(g2);
    
    % Welch's t formula
    se = sqrt(v1/n1 + v2/n2);
    t = (m1 - m2) / se;
end