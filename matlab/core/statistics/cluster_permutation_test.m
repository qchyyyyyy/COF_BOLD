function out = cluster_permutation_test(y, varargin)
% CLUSTER_PERM_1D  Cluster-based permutation test for a single time series.
% 目标：在单条时间序列 y 中，找出“显著 > 0”的时段与“显著 < 0”的时段（单侧）。
%
% 输入:
%   y           : 1 x T 或 T x 1 的时间序列（等间隔采样）
%
% 可选参数（Name-Value）:
%   'Alpha'     : 家族错误率 FWER 水平（默认 0.05）
%   'StartP'    : 形成集簇的单侧“起始阈值”p0（默认 0.05；越小越少碎片）
%   'B'         : 置换/替代次数（默认 4000；≥2000 更稳）
%   'MinLen'    : 最小集簇长度（默认 3 点）
%   'LocalStd'  : true/false，是否用滑窗标准差做局部标准化（默认 true）
%   'StdWin'    : 滑窗长度（默认 31；需为奇数）
%   'NullMode'  : 'shift' | 'aaft' | 'iaaft'（默认 'shift'）
%   'Seed'      : 随机种子（默认 []）
%
% 输出结构 out:
%   .pos.clusters   : K_pos x 2 的 [start_idx, end_idx]（显著 >0）
%   .pos.pvals      : 每个显著正向集簇的 p 值
%   .neg.clusters   : K_neg x 2 的 [start_idx, end_idx]（显著 <0）
%   .neg.pvals      : 每个显著负向集簇的 p 值
%   .params         : 所有参数与阈值、置换数等
%   .mass           : 真实数据各方向的集簇质量
%
% 参考：
%   - FieldTrip/MNE 的 cluster-based permutation 教程（集簇质量/最大统计校正） 
%   - Nichols & Holmes (2002) 非参数置换综述（FWER 控制）
%   - AAFT/IAAFT 相位随机化替代序列（保持功率谱/幅度分布）
%
% 作者：你可爱且严谨的 AI 助教（中文+英文关键词）

%% ------------- 参数 -------------
p = inputParser;
addParameter(p, 'Alpha',    0.05);
addParameter(p, 'StartP',   0.05);
addParameter(p, 'B',        4000);
addParameter(p, 'MinLen',   3);
addParameter(p, 'LocalStd', true);
addParameter(p, 'StdWin',   31);
addParameter(p, 'NullMode', 'shift'); % 'shift' | 'aaft' | 'iaaft'
addParameter(p, 'Seed',     []);
parse(p, varargin{:});
P   = p.Results;

if ~isempty(P.Seed), rng(P.Seed); end

y = y(:);
T = numel(y);
assert(T >= 10, '时间序列太短，至少需要 ~10 点。');

%% ------------- 第0步：标准化到“近似 t 统计量” -------------
if P.LocalStd
    s = local_std(y, P.StdWin); % 滑窗标准差
    t0 = y ./ s;
else
    s = std(y, 1); % 用总体标准差刻度（ddof=0，偏保守）
    t0 = y / s;
end

% 起始阈值：单侧未校正 p0 -> z/t 临界值（用正态近似）
tau = -norminv(P.StartP);   % 单侧：P(Z>tau)=p0 => tau = norminv(1-p0) = -norminv(p0)

%% ------------- 第1步：真实数据的集簇（正/负分开） -------------
[posClus, posMass] = find_clusters(t0, +tau, +1, P.MinLen); % >0 集簇
[negClus, negMass] = find_clusters(t0, -tau, -1, P.MinLen); % <0 集簇（注意阈值符号）

%% ------------- 第2步：构造零分布（保留自相关） -------------
maxNullPos = zeros(P.B,1);
maxNullNeg = zeros(P.B,1);

for b = 1:P.B
    switch lower(P.NullMode)
        case 'shift'
            yb = circshift(y, randi(T));
        case 'aaft'
            yb = aaft_surrogate(y);
        case 'iaaft'
            yb = iaaft_surrogate(y, 100); % 迭代次数可调
        otherwise
            error('NullMode 仅支持 shift/aaft/iaaft');
    end

    % 用相同的“噪声刻度”规则得到 t 值
    if P.LocalStd
        sb = local_std(yb, P.StdWin); 
        tb = yb ./ sb;
    else
        tb = yb / s;
    end

    [~, massBpos] = find_clusters(tb, +tau, +1, P.MinLen);
    [~, massBneg] = find_clusters(tb, -tau, -1, P.MinLen);

    maxNullPos(b) = (isempty(massBpos)) * 0 + (~isempty(massBpos)) * max(massBpos);
    maxNullNeg(b) = (isempty(massBneg)) * 0 + (~isempty(massBneg)) * max(massBneg);
end

%% ------------- 第3步：集簇 p 值（单侧，max-cluster 校正） -------------
posP = arrayfun(@(m) (1 + sum(maxNullPos >= m)) / (1 + P.B), posMass);
negP = arrayfun(@(m) (1 + sum(maxNullNeg >= m)) / (1 + P.B), negMass);

% 显著集簇筛选
posSig = posP < P.Alpha;
negSig = negP < P.Alpha;

out.pos.clusters = posClus(posSig, :);
out.pos.pvals    = posP(posSig);
out.neg.clusters = negClus(negSig, :);
out.neg.pvals    = negP(negSig);

out.params = P;
out.params.Threshold = tau;
out.mass.pos = posMass;
out.mass.neg = negMass;

end % main


%% ==================== 辅助函数 ====================

function s = local_std(y, win)
% 滑窗标准差（反射填充，稳健对局部方差变化）
if mod(win,2)==0, win = win+1; end
pad = floor(win/2);
y2 = [flipud(y(1:pad)); y; flipud(y(end-pad+1:end))];
s  = zeros(size(y));
for i = 1:numel(y)
    seg = y2(i:i+win-1);
    s(i) = std(seg, 1); % ddof=0（总体标准差），更保守
end
s(s==0) = median(s(s>0)); % 防止除零
end

function [clusters, masses] = find_clusters(t, thr, side, minlen)
    if side > 0
        idx = find(t >= thr);
    else
        idx = find(t <= thr);
    end
    clusters = zeros(0,2);
    masses   = zeros(0,1);
    if isempty(idx), return; end

    breaks = find(diff(idx) > 1);   % breaks 是列向量
    starts = [1; breaks + 1];       % 用分号做“竖着拼接”
    ends   = [breaks; numel(idx)];  % 同理

    for k = 1:numel(starts)
        seg = idx(starts(k):ends(k));
        if numel(seg) >= minlen
            clusters(end+1,:) = [seg(1), seg(end)]; %#ok<AGROW>
            if side > 0
                masses(end+1,1) = sum(t(seg));      %#ok<AGROW>
            else
                masses(end+1,1) = sum(-t(seg));     %#ok<AGROW>
            end
        end
    end
end

function ys = aaft_surrogate(y)
% AAFT：匹配原始幅度分布 + 随机相位（一次性）
% 参考：AAFT/IAAFT 文献；生产级可用更强的 iAAFT（见下）。
y = y(:);
T = numel(y);
% 将标准正态的秩映射到 y 的秩分布
[ysort, ir] = sort(y);
z = randn(T,1); z = sort(z);
[~, r2] = sort(ir); y0 = z(r2);

% 频域随机相位
Y0 = fft(y0);
ph = exp(1i*2*pi*rand(floor(T/2)-1,1));
Y = zeros(size(Y0));
Y(1) = Y0(1);
Y(2:floor(T/2)) = abs(Y0(2:floor(T/2))).*ph;
Y(end-floor(T/2)+2:end) = conj(Y(2:floor(T/2)));
if mod(T,2)==0, Y(T/2+1)=Y0(T/2+1); end

x = real(ifft(Y));

% 将 x 的秩映射回 y 的幅度分布
[xs, ix] = sort(x);
ys = zeros(T,1);
ys(ix) = ysort;
end

function ys = iaaft_surrogate(y, maxIter)
% 简化 iAAFT：迭代匹配功率谱与幅度分布
y = y(:);
T = numel(y);
ys = aaft_surrogate(y);
amp = abs(fft(y));
ys_sorted = sort(y);

for it = 1:maxIter
    % 1) 替换幅度谱
    YS = fft(ys);
    YS = amp .* exp(1i*angle(YS));
    x = real(ifft(YS));

    % 2) 映回幅度分布
    [xs, ix] = sort(x);
    tmp = zeros(T,1);
    tmp(ix) = ys_sorted;
    % 收敛性粗判
    if norm(tmp - ys) < 1e-8*norm(ys)
        ys = tmp; break;
    end
    ys = tmp;
end
end
