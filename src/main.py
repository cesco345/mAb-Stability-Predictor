import os
from pdb_parser.py import extract_residues
from feature_extractor import extract_features
import pandas as pd

# Step 1: Extract residue coordinates and IDs
pdb_file = 'data/input_pdbs/antibody.pdb'
residue_coords, residue_ids = extract_residues(pdb_file)

# Step 2: Extract features using PyRosetta
features = extract_features(pdb_file)

# Step 3: Combine features and coordinates
feature_data = []
for coord, feature, residue_id in zip(residue_coords, features, residue_ids):
    feature_data.append(coord + feature + [residue_id])

# Step 4: Save feature data to CSV
feature_df = pd.DataFrame(feature_data, columns=['x', 'y', 'z', 'sasa', 'packstat', 'residue_id'])
feature_df.to_csv('data/feature_data.csv', index=False)

# Step 5: Train the model
from train_model import train_model

train_model('data/feature_data.csv', 'data/stability_scores.csv', 'models/stability_model.pth')
