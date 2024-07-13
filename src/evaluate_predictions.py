import pandas as pd
from sklearn.metrics import mean_squared_error


def evaluate_predictions(predictions_file, true_labels_file):
    predictions = pd.read_csv(predictions_file)
    true_labels = pd.read_csv(true_labels_file)

    merged = pd.merge(predictions, true_labels, on='residue_id', suffixes=('_pred', '_true'))
    print(merged.columns)  # Debugging: Print the columns of the merged DataFrame to verify the suffixes
    mse = mean_squared_error(merged['predicted_stability_score'], merged['stability_score'])

    print(f'Mean Squared Error: {mse:.4f}')


if __name__ == '__main__':
    predictions_file = '../data/predictions.csv'
    true_labels_file = '../data/true_labels.csv'  # Replace with the true labels file path

    evaluate_predictions(predictions_file, true_labels_file)

