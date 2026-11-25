import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import scipy.io as sio
import networkx as nx
import os
from tqdm import tqdm
from scipy.spatial import distance_matrix
import argparse


# ---------------------------------------------------------
# 1. 认知任务生成器 (带上下文线索)
# ---------------------------------------------------------
def generate_cognitive_data(seq_len, task_type='2-back', num_stimuli=5, batch_size=1):
    input_dim = num_stimuli + 2
    Y = np.zeros((batch_size, seq_len), dtype=np.float32)
    X_features = np.zeros((batch_size, seq_len, input_dim), dtype=np.float32)
    stimuli_seq = np.zeros((batch_size, seq_len), dtype=int)

    if task_type == 'rest':
        X_features = np.random.normal(0, 0.05, size=(batch_size, seq_len, input_dim)).astype(np.float32)
        return torch.tensor(X_features), torch.tensor(Y).unsqueeze(-1), stimuli_seq

    TARGET_RATIO = 0.35

    for b in range(batch_size):
        if task_type == '0-back':
            target_stim = 0
            num_targets = int(seq_len * TARGET_RATIO)
            target_indices = np.random.choice(np.arange(seq_len), num_targets, replace=False)

            for t in range(seq_len):
                if t in target_indices:
                    stimuli_seq[b, t] = target_stim
                    Y[b, t] = 1.0
                else:
                    choices = [s for s in range(num_stimuli) if s != target_stim]
                    stimuli_seq[b, t] = np.random.choice(choices)
            X_features[b, :, num_stimuli] = 1.0

        elif task_type == '2-back':
            stimuli_seq[b, 0] = np.random.randint(0, num_stimuli)
            stimuli_seq[b, 1] = np.random.randint(0, num_stimuli)
            num_targets = int((seq_len - 2) * TARGET_RATIO)
            target_indices = np.random.choice(np.arange(2, seq_len), num_targets, replace=False)

            for t in range(2, seq_len):
                if t in target_indices:
                    stimuli_seq[b, t] = stimuli_seq[b, t - 2]
                    Y[b, t] = 1.0
                else:
                    avoid_stim = stimuli_seq[b, t - 2]
                    choices = [s for s in range(num_stimuli) if s != avoid_stim]
                    stimuli_seq[b, t] = np.random.choice(choices)
            X_features[b, :, num_stimuli + 1] = 1.0

        for t in range(seq_len):
            X_features[b, t, stimuli_seq[b, t]] = 1.0

    return torch.tensor(X_features), torch.tensor(Y).unsqueeze(-1), stimuli_seq


# ---------------------------------------------------------
# 2. 生物物理模型：Leaky RNN
# ---------------------------------------------------------
class LeakyRNN(nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim, alpha=0.2):
        super(LeakyRNN, self).__init__()
        self.hidden_dim = hidden_dim
        self.alpha = alpha

        self.W_ih = nn.Linear(input_dim, hidden_dim, bias=False)
        self.W_hh = nn.Linear(hidden_dim, hidden_dim, bias=False)
        self.fc_out = nn.Linear(hidden_dim, output_dim)

        self.b = nn.Parameter(torch.zeros(hidden_dim))

        # 物理正交初始化与信号放大 （输入层用 Xavier 保证方差一致）
        nn.init.orthogonal_(self.W_hh.weight)
        nn.init.xavier_uniform_(self.W_ih.weight, gain=3.0)

        # 贝叶斯偏置初始化
        nn.init.constant_(self.fc_out.bias, -0.619)

    def forward(self, x, h0=None, return_diagnostics=False):
        """Run the leaky RNN.

        Parameters
        ----------
        x : torch.Tensor
            Input sequence with shape [batch, time, input_dim].
        h0 : torch.Tensor or None
            Optional initial hidden state with shape [batch, hidden_dim].
            If None, the same stochastic initialization used in the original
            experiments is used.
        return_diagnostics : bool
            If True, also return h0_used, preactivation history, hidden states
            before each update, and readout logits. These variables are saved
            for later exact intervention / reconstruction analyses.
        """
        batch_size, seq_len, _ = x.size()
        device = x.device

        # 引入基线噪声波动，模拟人脑初始状态的差异性
        if h0 is None:
            h_t = torch.randn(batch_size, self.hidden_dim, device=device, dtype=x.dtype) * 0.05
        else:
            h_t = h0.to(device=device, dtype=x.dtype)

        h0_used = h_t.clone()
        h_history = []
        h_pre_update_history = []
        preactivation_history = []

        for t in range(seq_len):
            h_pre_update_history.append(h_t.unsqueeze(1))
            rnn_input = self.W_ih(x[:, t, :]) + self.W_hh(h_t) + self.b
            preactivation_history.append(rnn_input.unsqueeze(1))
            h_t = (1 - self.alpha) * h_t + self.alpha * torch.tanh(rnn_input)
            h_history.append(h_t.unsqueeze(1))

        out = torch.cat(h_history, dim=1)
        logits = self.fc_out(out)
        predictions = torch.sigmoid(logits)

        if return_diagnostics:
            diagnostics = {
                'h0': h0_used,
                'hidden_pre_update': torch.cat(h_pre_update_history, dim=1),
                'preactivations': torch.cat(preactivation_history, dim=1),
                'logits': logits,
            }
            return predictions, out, diagnostics

        return predictions, out


# ---------------------------------------------------------
# 3. 核心训练引擎
# ---------------------------------------------------------
def train_and_extract(model_type, device, config, model_idx):
    print(f"\n" + "=" * 50)
    print(f"🚀 启动实验分支: {model_type.upper()} RNN ({config['HIDDEN_DIM']} Neurons) - 被试/模型 {model_idx}")
    print("=" * 50)

    model = LeakyRNN(input_dim=config['INPUT_DIM'],
                     hidden_dim=config['HIDDEN_DIM'],
                     output_dim=1, alpha=0.2).to(device)

    optimizer = optim.Adam(model.parameters(), lr=0.001)
    criterion = nn.BCELoss(reduction='none')

    D_tensor = None
    spatial_coords = None
    if model_type == 'spatial':
        grid_rows = int(np.sqrt(config['HIDDEN_DIM']))
        grid_cols = int(np.sqrt(config['HIDDEN_DIM']))
        x_vals = np.linspace(0, 1, grid_cols)
        y_vals = np.linspace(0, 1, grid_rows)
        xx, yy = np.meshgrid(x_vals, y_vals)
        spatial_coords = np.column_stack((xx.ravel(), yy.ravel())).astype(np.float32)
        dist_matrix = distance_matrix(spatial_coords, spatial_coords)
        beta = 3.0
        penalty_matrix = np.exp(beta * dist_matrix)
        D_tensor = torch.tensor(penalty_matrix, dtype=torch.float32, device=device)

    print(f"-> 开始多任务认知训练 ({config['EPOCHS']} Epochs)...")
    model.train()
    with tqdm(total=config['EPOCHS']) as pbar:
        history = {'0-back': {'loss': [], 'acc': []}, '2-back': {'loss': [], 'acc': []}}

        for epoch in range(config['EPOCHS']):
            if epoch < 2000:
                task = '2-back'
            else:
                task = '0-back' if np.random.rand() > 0.5 else '2-back'

            x_train, y_train, _ = generate_cognitive_data(
                config['SEQ_LEN_TRAIN'], task_type=task,
                num_stimuli=config['NUM_STIMULI'], batch_size=128
            )
            x_train, y_train = x_train.to(device), y_train.to(device)

            optimizer.zero_grad()
            predictions, _ = model(x_train)

            # 代价敏感学习退火
            current_weight = 1.0 + 1.5 * min(epoch / 500.0, 1.0)
            base_loss = criterion(predictions, y_train)
            weight_tensor = torch.where(y_train == 1.0, current_weight, 1.0)
            loss = (base_loss * weight_tensor).mean()

            if model_type == 'spatial':
                if epoch < 1000:
                    current_lambda = 0.0
                elif epoch < 2000:
                    current_lambda = 0.5 * ((epoch - 1000) / 1000.0)
                else:
                    current_lambda = 0.5
                wiring_cost = torch.mean(torch.abs(model.W_hh.weight) * D_tensor)
                loss = loss + current_lambda * wiring_cost

            loss.backward()
            dynamic_max_norm = np.sqrt(config['HIDDEN_DIM'] / 256.0)
            nn.utils.clip_grad_norm_(model.parameters(), max_norm=dynamic_max_norm)
            optimizer.step()

            acc = ((predictions > 0.5).float() == y_train).float().mean().item()
            history[task]['loss'].append(loss.item())
            history[task]['acc'].append(acc)

            acc_0 = np.mean(history['0-back']['acc'][-2:]) if len(history['0-back']['acc']) > 0 else 0
            acc_2 = np.mean(history['2-back']['acc'][-2:]) if len(history['2-back']['acc']) > 0 else 0
            pbar.set_postfix({'0-back Acc ': f"{acc_0:.3f}", ' 2-back Acc ': f"{acc_2:.3f}"})
            pbar.update(1)

    print("-> 通过 FD 算法事后提取二维坐标...")
    W_hh = model.W_hh.weight.detach().cpu().numpy()
    W_abs = np.abs(W_hh)
    threshold = np.percentile(W_abs, config['SPARSITY'] * 100)
    A_sparse = np.where(W_abs > threshold, W_abs, 0.0)
    G = nx.from_numpy_array(A_sparse, create_using=nx.DiGraph)
    # 【注意】这里解除了刚性随机种子(seed=42)，利用 model_idx 保证每个模型的投射互相独立！
    pos_dict = nx.spring_layout(G, weight='weight', iterations=150, seed=model_idx)

    xy_coords = np.zeros((config['HIDDEN_DIM'], 2), dtype=np.float32)
    for i in range(config['HIDDEN_DIM']):
        xy_coords[i, :] = pos_dict[i]

    # --- 功能时空矩阵与完整突触信息提取 ---
    print(f"-> 提取 {config['NUM_TEST_TRIALS']} 个 Trial 的脑功能状态时间序列...")
    model.eval()

    data_dict = {
        # Full trained dynamics parameters. The original version saved weights
        # but omitted the learned recurrent bias and readout bias; both are
        # required for exact forward reconstruction and intervention analysis.
        'W_ih': model.W_ih.weight.detach().cpu().numpy(),
        'W_hh': W_hh,
        'W_out': model.fc_out.weight.detach().cpu().numpy(),
        'b': model.b.detach().cpu().numpy(),
        'fc_out_bias': model.fc_out.bias.detach().cpu().numpy(),
        'alpha': np.array([[model.alpha]], dtype=np.float32),

        # Spatial embedding and metadata.
        'coordinates_2d_FD': xy_coords,
        'model_type': model_type,
        'subject_id': model_idx,
        'training_history': history,
        'num_stimuli': np.array([[config['NUM_STIMULI']]], dtype=np.int32),
        'input_dim': np.array([[config['INPUT_DIM']]], dtype=np.int32),
        'hidden_dim': np.array([[config['HIDDEN_DIM']]], dtype=np.int32),
        'seq_len_test': np.array([[config['SEQ_LEN_TEST']]], dtype=np.int32),
        'num_test_trials': np.array([[config['NUM_TEST_TRIALS']]], dtype=np.int32),
    }

    if model_type == 'spatial':
        data_dict['coordinates_2d'] = spatial_coords

    for task in ['rest', '0-back', '2-back']:
        x_test, y_test, stimuli_seq = generate_cognitive_data(
            config['SEQ_LEN_TEST'], task_type=task,
            num_stimuli=config['NUM_STIMULI'], batch_size=config['NUM_TEST_TRIALS']
        )
        x_test = x_test.to(device)

        # Save the RNG state just before the stochastic h0 draw. This is
        # diagnostic only; h0_used below is the authoritative initial state.
        torch_rng_state_before = torch.get_rng_state().cpu().numpy()
        if device.type == 'cuda':
            cuda_rng_state_before = torch.cuda.get_rng_state(device).cpu().numpy()
        else:
            cuda_rng_state_before = np.array([], dtype=np.uint8)

        with torch.no_grad():
            predictions, hidden_states, diagnostics = model(x_test, return_diagnostics=True)

        # Transpose 为 [Neurons, Time, Trials] 方便 MATLAB 后续读取
        activations = hidden_states.cpu().numpy().transpose(2, 1, 0)
        preactivations = diagnostics['preactivations'].cpu().numpy().transpose(2, 1, 0)
        hidden_pre_update = diagnostics['hidden_pre_update'].cpu().numpy().transpose(2, 1, 0)
        h0_used = diagnostics['h0'].cpu().numpy().T  # [Neurons, Trials]
        logits = diagnostics['logits'].cpu().numpy().squeeze(-1)  # [Trials, Time]
        probabilities = predictions.cpu().numpy().squeeze(-1)    # [Trials, Time]

        prefix = task.replace('-', '')
        data_dict[f'activations_{prefix}'] = activations
        data_dict[f'hidden_pre_update_{prefix}'] = hidden_pre_update
        data_dict[f'preactivations_{prefix}'] = preactivations
        data_dict[f'h0_{prefix}'] = h0_used
        data_dict[f'x_features_{prefix}'] = x_test.detach().cpu().numpy()
        data_dict[f'stimuli_{prefix}'] = stimuli_seq
        data_dict[f'true_targets_{prefix}'] = y_test.cpu().numpy().squeeze()
        data_dict[f'predictions_{prefix}'] = probabilities
        data_dict[f'logits_{prefix}'] = logits
        data_dict[f'model_responses_{prefix}'] = (probabilities > 0.5).astype(float)
        data_dict[f'torch_rng_state_before_{prefix}'] = torch_rng_state_before
        data_dict[f'cuda_rng_state_before_{prefix}'] = cuda_rng_state_before

    num_neuron = config['HIDDEN_DIM']
    # 动态命名，防止覆盖 (如: activation_dynamics_baseline_1024_sub01.mat)
    save_path = os.path.join(config['SAVE_DIR'],
                             f'activation_dynamics_{model_type}_{num_neuron}_sub{model_idx:02d}.mat')
    checkpoint_path = save_path.replace('.mat', '.pth')
    data_dict['checkpoint_path'] = checkpoint_path

    # A PyTorch checkpoint is saved in parallel for exact model restoration.
    # The .mat file contains all numeric arrays needed by MATLAB/Python COF and
    # intervention code; the .pth file keeps the complete state_dict.
    torch.save({
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'config': config,
        'history': history,
        'model_type': model_type,
        'subject_id': model_idx,
    }, checkpoint_path)

    sio.savemat(save_path, data_dict, do_compression=True)
    print(f"✅ 模型 {model_idx} 的数据已保存至: {save_path}")
    print(f"✅ PyTorch checkpoint 已保存至: {checkpoint_path}")


# ---------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="activation dynamics of Leaky RNN performing working memory tasks")
    parser.add_argument('--model', type=str, default='baseline', choices=['baseline', 'spatial'],
                        help="Choose the method for spatial location.")
    parser.add_argument('--num_neuron', type=int, default=1024,
                        choices=[144, 256, 400, 576, 784, 1024, 1296, 1600, 1936, 2304, 4096, 16384],
                        help="Choose the dataset.")
    # 【新增参数】：模拟群体分析的样本量
    parser.add_argument('--num_models', type=int, default=10,
                        help="Number of independent models (subjects) to train.")
    parser.add_argument('--output-dir', default='results/rnn_training', help='Fresh directory for synthetic model outputs.')
    parser.add_argument('--seed', type=int, default=None, help='Optional new-run seed; original training seeds were not recorded.')
    args = parser.parse_args()
    if args.seed is not None:
        np.random.seed(args.seed); torch.manual_seed(args.seed)

    DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    config = {
        'NUM_STIMULI': 5,
        'INPUT_DIM': 7,
        'HIDDEN_DIM': args.num_neuron,
        'SEQ_LEN_TRAIN': 60,
        'SEQ_LEN_TEST': 200,
        'NUM_TEST_TRIALS': 3,  # 缩小 Trial 数量，3次足以平滑噪声
        'EPOCHS': 4000,
        'SPARSITY': 0.1,
        'SAVE_DIR': args.output_dir
    }
    os.makedirs(config['SAVE_DIR'], exist_ok=True)

    # 循环训练多个独立的被试模型
    for m in range(1, args.num_models + 1):
        train_and_extract(model_type=args.model, device=DEVICE, config=config, model_idx=m)

    print(f"\n🎉 规模为 {args.num_neuron} 的 {args.num_models} 个被试群体实验全部完成！")


if __name__ == '__main__':
    main()
