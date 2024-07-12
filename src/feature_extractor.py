import pyrosetta
from pyrosetta import rosetta

pyrosetta.init()

def extract_features(pdb_file):
    pose = pyrosetta.pose_from_pdb(pdb_file)
    features = []

    # Initialize SASA calculation
    sasa_calculator = rosetta.core.scoring.sasa.SasaCalc()
    sasa_calculator.calculate(pose)
    sasa_values = sasa_calculator.get_residue_sasa()

    # Initialize BuriedUnsatHbondFilter
    buried_unsat_filter = rosetta.protocols.simple_filters.BuriedUnsatHbondFilter()

    for i in range(1, pose.size() + 1):
        res = pose.residue(i)
        if res.name1() in ['Y', 'W']:
            # Get SASA for the specific residue
            sasa = sasa_values[i]
            # Calculate Buried Unsatisfied Hydrogen Bonds (as an example feature)
            buried_unsat = buried_unsat_filter.compute(pose)
            features.append([sasa, buried_unsat])

    return features



