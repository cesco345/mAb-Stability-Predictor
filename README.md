# mAb-Stability-Predictor

## Overview

**mAb-Stability-Predictor** is a computational pipeline designed to predict the stability of monoclonal antibodies (mAbs) using deep learning with PyTorch and structural analysis with PyRosetta. Stability is a critical factor in the development of effective and safe therapeutic antibodies. This tool aims to facilitate drug discovery by providing accurate stability predictions, thereby aiding in the design and selection of stable antibody candidates.

## Features

- **PDB Parsing**: Extracts residue coordinates and identifiers from PDB files.
- **Feature Extraction**: Uses PyRosetta to calculate structural features such as solvent-accessible surface area (SASA) and packstat.
- **Deep Learning Model**: Trains a neural network using PyTorch to predict stability scores based on structural features.
- **Data Integration**: Combines experimental stability scores with computational features for robust model training.

## Directory Structure



## Getting Started

### Prerequisites

- Python 3.6+
- PyTorch
- PyRosetta
- Biopython
- Pandas
- NumPy
- Jupyter

### Installation

1. **Clone the repository**:
    ```bash
    git clone https://github.com/cesco345/mAb-Stability-Predictor.git
    cd mAb-Stability-Predictor
    ```

2. **Install dependencies**:
    ```bash
    pip install -r requirements.txt
    ```

3. **Download and install PyRosetta**:
    Follow the instructions on the [PyRosetta website](https://www.pyrosetta.org/dow) to download and install PyRosetta.

### Usage

1. **Prepare Input Data**:
    - Place your antibody PDB files in the `data/input_pdbs/` directory.
    - Create a `stability_scores.csv` file in the `data/` directory with stability scores for each residue.

2. **Run the Pipeline**:
    - Use the main script to extract features, combine data, and train the model:
        ```bash
        python src/main.py
        ```

3. **Analyze Results**:
    - Open the Jupyter notebook in the `notebooks/` directory to analyze and visualize the model's predictions:
        ```bash
        jupyter notebook notebooks/analysis.ipynb
        ```

## Example

Here’s an example of the `stability_scores.csv` file format:

```csv
residue_id,stability_score
TYR_15,0.90
TYR_22,0.85
TRP_10,0.88
TRP_34,0.80
```

## Contributing

We welcome contributions to enhance the functionality and usability of this tool. Please fork the repository, create a new branch, and submit a pull request with your changes.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgements

- [PyRosetta](https://www.pyrosetta.org/)
- [PyTorch](https://pytorch.org/)
- [Biopython](https://biopython.org/)

## Contact

For questions or issues, please open an issue on GitHub or contact the repository maintainer at [fpiscani@stem-apks.com](mailto:fpiscani@stem-apks.com).

