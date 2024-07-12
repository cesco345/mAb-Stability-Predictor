import os
from pdb_parser import extract_residues
from feature_extractor import extract_features
import pandas as pd

# Step 1: Extract residue coordinates and IDs
pdb_file = '../data/input_pdbs/5jxr.pdb'
residue_coords, residue_ids = extract_residues(pdb_file)

# Step 2: Extract features using PyRosetta
features = extract_features(pdb_file)

# Step 3: Combine features and coordinates
feature_data = []
for coord, feature, residue_id in zip(residue_coords, features, residue_ids):
    feature_data.append(coord.tolist() + feature + [residue_id])

# Ensure the output directory exists
output_dir = '../data'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Step 4: Save feature data to CSV
feature_df = pd.DataFrame(feature_data, columns=['x', 'y', 'z', 'sasa', 'buried_unsat', 'residue_id'])
feature_df.to_csv(os.path.join(output_dir, 'feature_data.csv'), index=False)

# Check number of rows in both files
feature_df = pd.read_csv(os.path.join(output_dir, 'feature_data.csv'))
label_df = pd.read_csv('../data/stability_scores.csv')
print(f"Number of features: {len(feature_df)}")
print(f"Number of labels: {len(label_df)}")

# Step 5: Train the model
from train_model import train_model

train_model(os.path.join(output_dir, 'feature_data.csv'), '../data/stability_scores.csv', '../models/stability_model.pth')



