import torch
import torch.nn as nn
import torch.optim as optim

class StabilityPredictor(nn.Module):
    def __init__(self):
        super(StabilityPredictor, self).__init__()
        self.fc1 = nn.Linear(5, 64)  # 3 coords + 2 features
        self.fc2 = nn.Linear(64, 32)
        self.fc3 = nn.Linear(32, 1)

    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        x = self.fc3(x)
        return x
