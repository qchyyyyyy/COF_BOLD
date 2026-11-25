function [p_val, obs_ks, perm_dist] = permutation_ks_test(x, y, n_perms)
% PERM_KS_TEST_IND 基于 K-S 统计量的置换检验 (高性能优化版)
%
% 输入:
%   x, y    : 两个独立样本 (维度可不一致)
%   n_perms : 置换次数
%
% 输出:
%   p_val     : P 值
%   obs_ks    : 原始数据的 K-S 统计量 (最大垂直距离)
%   perm_dist : 零分布
%
% 数学原理: D = max | F_n(x) - G_m(x) |
% 优化策略: 预先对合并数据排序，循环内只打乱标签计算累积和，避免重复排序。

    if nargin < 3, n_perms = 5000; end

    x = x(:); nx = length(x);
    y = y(:); ny = length(y);
    n_total = nx + ny;
    
    % --- 1. 计算原始 K-S 统计量 ---
    % 使用 MATLAB 内置函数做基准
    [~, ~, obs_ks] = kstest2(x, y);
    
    % --- 2. 预处理 (核心优化步骤) ---
    % 将所有数据合并并固定排序。
    % 在 K-S 检验中，数值本身不重要，重要的是它们的相对秩次(Rank)。
    % 我们只需要知道在排序后的序列中，哪些点属于X，哪些属于Y。
    combined_data = [x; y];
    [~, sort_idx] = sort(combined_data);
    
    % 创建一个逻辑索引模板，对应排序后的数据
    % 真实的标签是: 原数据的第 sort_idx 个位置来源于哪里
    % 但在置换检验中，我们不需要追踪原始位置，只需要随机分配 n_x 个位置给 X
    
    perm_dist = zeros(n_perms, 1);
    
    % 预计算分母，避免循环内除法
    inv_nx = 1 / nx;
    inv_ny = 1 / ny;
    
    % --- 3. 高速置换循环 ---
    % 即使 n_total = 10000, n_perms = 10000 也能在数秒内完成
    for i = 1:n_perms
        % 随机生成一个逻辑向量，表示排序后的序列中，哪些被分配给了 "伪组 X"
        % randsample 或 randperm 生成 1~n_total 的索引
        perm_indices = randperm(n_total, nx);
        
        % 构造逻辑向量 (boolean mask)
        is_x = false(n_total, 1);
        is_x(perm_indices) = true;
        
        % 计算两个经验分布函数 (ECDF)
        % 因为数据已经是排好序的，cumsum 直接就是未归一化的 CDF
        cdf_x = cumsum(is_x) * inv_nx;
        cdf_y = cumsum(~is_x) * inv_ny;
        
        % 计算 K-S 统计量: 两个 CDF 差值的最大绝对值
        perm_dist(i) = max(abs(cdf_x - cdf_y));
    end
    
    % --- 4. 计算 P 值 ---
    % K-S 检验天生是单尾的 (统计量 D >= 0，越大越显著)
    p_val = (sum(perm_dist >= obs_ks) + 1) / (n_perms + 1);
    
    % % 可视化
    % figure;
    % histogram(perm_dist, 50, 'Normalization', 'pdf', 'EdgeColor', 'none');
    % hold on;
    % xline(obs_ks, 'r', 'LineWidth', 2);
    % title(['K-S Permutation Test (p = ' num2str(p_val, '%.4f') ')']);
    % subtitle('Null Distribution of D_{stat}');
    % xlabel('K-S Statistic (D)');
    % box off;
end