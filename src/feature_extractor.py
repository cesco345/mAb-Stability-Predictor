import pyrosetta

pyrosetta.init()

def extract_features(pdb_file):
    pose = pyrosetta.pose_from_pdb(pdb_file)
    features = []

    for i in range(1, pose.size() + 1):
        res = pose.residue(i)
        if res.name1() in ['Y', 'W']:
            sasa = pyrosetta.rosetta.protocols.simple_filters.SasaFilter().report_sm(pose)
            packstat = pyrosetta.rosetta.protocols.simple_filters.PackStatFilter().report_sm(pose)
            features.append([sasa, packstat])

    return features
