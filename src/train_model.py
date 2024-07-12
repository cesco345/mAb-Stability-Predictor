import pandas as pd
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
from stability_predictor import StabilityPredictor
import torch.nn as nn
import torch.optim as optim

class StabilityDataset(Dataset):
    def __init__(self, feature_file, label_file):
        self.features = pd.read_csv(feature_file).values
        self.labels = pd.read_csv(label_file)['stability_score'].values.reshape(-1, 1)

    def __len__(self):
        return len(self.features)

    def __getitem__(self, idx):
        return torch.tensor(self.features[idx], dtype=torch.float32), torch.tensor(self.labels[idx], dtype=torch.float32)

def train_model(feature_file, label_file, model_path):
    dataset = StabilityDataset(feature_file, label_file)
    dataloader = DataLoader(dataset, batch_size=32, shuffle=True)

    model = StabilityPredictor()
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)

    num_epochs = 1000
    for epoch in range(num_epochs):
        model.train()
        for features, labels in dataloader:
            optimizer.zero_grad()
            outputs = model(features)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()

        if (epoch + 1) % 100 == 0:
            print(f'Epoch [{epoch + 1}/{num_epochs}], Loss: {loss.item():.4f}')

    torch.save(model.state_dict(), model_path)
