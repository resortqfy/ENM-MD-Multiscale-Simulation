import pandas as pd
import MDAnalysis as mda

def write_pdb(coords: np.ndarray, output_file: str, template: mda.Universe):
    """写入PDB文件"""
    with mda.Writer(output_file) as w:
        for pos in coords:
            template.atoms.positions = pos
            w.write(template)

def write_csv(data: dict, output_file: str):
    """写入CSV"""
    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)

