import torch
from stability_predictor import StabilityPredictor

def load_model(model_path):
    model = StabilityPredictor()
    model.load_state_dict(torch.load(model_path))
    model.eval()
    return model
