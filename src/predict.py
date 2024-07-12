import pandas as pd
import torch
from load_model import load_model
import numpy as np

def prepare_data(file_path):
    df = pd.read_csv(file_path)
    features = df.drop(columns=['residue_id']).values.astype(np.float32)
    return torch.tensor(features, dtype=torch.float32), df['residue_id'].values

def predict(model, features):
    with torch.no_grad():
        predictions = model(features)
    return predictions.numpy()

def main():
    model_path = '../models/stability_model.pth'
    new_data_path = '../data/new_feature_data.csv'  # Replace with your new data file path

    model = load_model(model_path)
    features, residue_ids = prepare_data(new_data_path)

    predictions = predict(model, features)
    result_df = pd.DataFrame({'residue_id': residue_ids, 'predicted_stability_score': predictions.flatten()})
    result_df.to_csv('../data/predictions.csv', index=False)

if __name__ == '__main__':
    main()
