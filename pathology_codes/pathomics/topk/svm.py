import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable

class SmoothTop1SVM(nn.Module):
    def __init__(self, n_classes, alpha=1., tau=1.0):
        super(SmoothTop1SVM, self).__init__()
        self.n_classes = n_classes
        self.alpha = alpha
        self.tau = tau

    def forward(self, x, y):
        """
        x: output logits (before softmax)
        y: ground truth labels (indices)
        """
        if y.dim() == 1:
            # 将标签转换为 one-hot 编码
            y = F.one_hot(y, num_classes=self.n_classes).float()
        
        # 计算 Smooth Top-1 Loss
        # 公式参考: "Smooth Loss Functions for Deep Top-k Classification"
        
        # 1. 提取正确类别的分数
        y_true = (x * y).sum(dim=1, keepdim=True)
        
        # 2. 提取非正确类别的分数
        y_false = x - 1e12 * y # 屏蔽掉正确类别的分数
        
        # 3. Smooth Max approximation
        # using LogSumExp trick for numerical stability
        y_topk = self.tau * torch.logsumexp(y_false / self.tau, dim=1, keepdim=True)
        
        # 4. Hinge Loss
        loss = F.relu(y_topk - y_true + self.alpha)
        
        return loss.mean()