import numpy as np
import matplotlib.pyplot as plt
import argparse
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
import os
import sys

def calculate_and_plot_per_residue_rmsf(ref_file, traj_file, start_res, end_res, atom_selection, output_filename=None):
    """
    Calculates average per-residue RMSF over a trajectory and plots the results as a bar chart.

    Args:
        ref_file (str): Path to the reference structure file (.pdb, .gro, etc.) for topology.
        traj_file (str): Path to the trajectory file (.xtc, .dcd, .trr, etc.).
        start_res (int): The starting residue number (inclusive) for calculation.
        end_res (int): The ending residue number (inclusive) for calculation.
        atom_selection (str): MDAnalysis selection string for atoms *within each residue*
                              to include in the RMSF calculation (e.g., 'backbone', 'all', 'name CA').
        output_filename (str, optional): Path to save the output plot image.
                                         If None, the plot is displayed.
    """
    try:
        # Load the universe
        print(f"Loading universe from topology: {ref_file} and trajectory: {traj_file}")
        u = mda.Universe(ref_file, traj_file)

        # Define selection for alignment (e.g., protein backbone)
        # This selection determines which atoms are used to superimpose frames
        # to remove overall translation/rotation. Choose a stable region.
        # We will align to the first frame's coordinates for this selection.
        alignment_selection = 'protein and backbone'
        print(f"Using '{alignment_selection}' for trajectory alignment.")

        # Perform trajectory alignment
        # This modifies the coordinates of the trajectory in place.
        # Align to the first frame (index 0) of the trajectory.
        try:
            print("Aligning trajectory...")
            # Create a reference atom group from the first frame's selection
            ref_atoms = u.select_atoms(alignment_selection)
            # Ensure the trajectory is at the first frame before aligning
            u.trajectory[0]
            aligner = align.AlignTraj(u, ref_atoms, select=alignment_selection, in_memory=True).run()
            print("Alignment complete.")
        except Exception as e:
             print(f"Error during trajectory alignment: {e}")
             print("Alignment is crucial for meaningful RMSF calculation.")
             return

        # Define the selection for which atoms to calculate RMSF for
        # This selects all atoms in the specified residue range that match atom_selection
        rmsf_selection_string = f"resid {start_res}:{end_res} and {atom_selection}"
        rmsf_atoms = u.select_atoms(rmsf_selection_string)

        if len(rmsf_atoms) == 0:
            print(f"Error: No atoms selected for RMSF calculation with selection '{rmsf_selection_string}'.", file=sys.stderr)
            return

        print(f"Calculating RMSF for selected atoms ({len(rmsf_atoms)} atoms): '{rmsf_selection_string}'...")

        # Calculate atomic RMSF for the selected atoms
        # The calculation happens on the aligned trajectory
        R = rms.RMSF(rmsf_atoms).run()

        # The result R.results.rmsf is a numpy array of RMSF values, one for each atom in rmsf_atoms
        atomic_rmsf_values = R.results.rmsf

        # --- Process atomic RMSF to get per-residue average RMSF ---
        print("Calculating average per-residue RMSF...")

        # Create a dictionary to store summed RMSF and atom counts for each residue
        residue_rmsf_sum = {}
        residue_atom_count = {}

        # Iterate through the atoms that were included in the RMSF calculation
        # and assign their RMSF value to their parent residue
        for i, atom in enumerate(rmsf_atoms):
            res_id = atom.resnum
            # Ensure the residue is within the target range, although the selection string should handle this
            if start_res <= res_id <= end_res:
                 if res_id not in residue_rmsf_sum:
                     residue_rmsf_sum[res_id] = 0.0
                     residue_atom_count[res_id] = 0
                 # Add the atom's RMSF to the residue's sum
                 residue_rmsf_sum[res_id] += atomic_rmsf_values[i]
                 residue_atom_count[res_id] += 1

        # Calculate the average RMSF for each residue
        # Store results in lists for plotting
        residue_numbers = []
        average_per_residue_rmsf = []

        # Iterate through the residues in the specified range to maintain order
        for res_id in range(start_res, end_res + 1):
             if res_id in residue_rmsf_sum and residue_atom_count[res_id] > 0:
                 avg_rmsf = residue_rmsf_sum[res_id] / residue_atom_count[res_id]
                 residue_numbers.append(res_id)
                 average_per_residue_rmsf.append(avg_rmsf)
             else:
                 # This residue might exist in the range but had no atoms matching the selection
                 print(f"Warning: No atoms matching selection '{atom_selection}' found for residue {res_id}.")
                 # Optionally append NaN or 0 if you want to show a gap/zero bar
                 # For a clean bar plot, we only add residues that had valid calculations
                 pass


        if not residue_numbers:
             print("No valid per-residue RMSF data calculated. Cannot plot.", file=sys.stderr)
             return

        # --- Plotting ---
        print("Generating plot...")
        plt.figure(figsize=(15, 7)) # Adjust figure size as needed

        # Create the bar plot
        bars = plt.bar(residue_numbers, average_per_residue_rmsf, width=0.8, color='lightcoral')

        plt.xlabel("Residue Number")
        plt.ylabel(f"Average RMSF ({atom_selection}) ($\AA$)") # Add units, usually Angstroms
        plt.title(f"Average Per-Residue {atom_selection} RMSF (Residues {start_res}-{end_res})")
        plt.grid(axis='y', linestyle='--', alpha=0.7) # Add horizontal grid lines

        # Improve x-axis ticks visibility if many residues
        num_residues_in_range = end_res - start_res + 1
        if num_residues_in_range <= 30:
             plt.xticks(residue_numbers) # Show every residue number that was plotted
        elif num_residues_in_range <= 80:
             # Select ticks from the list of plotted residue numbers
             tick_positions = residue_numbers[::2]
             plt.xticks(tick_positions, tick_positions)
        else:
             tick_positions = residue_numbers[::5]
             plt.xticks(tick_positions, tick_positions)

        # Adjust x-axis limits slightly for better visual
        plt.xlim(start_res - 0.5, end_res + 0.5)

        plt.tight_layout() # Adjust layout to prevent labels overlapping

        # Save or show the plot
        if output_filename:
            try:
                plt.savefig(output_filename, dpi=300) # Save with higher resolution
                print(f"Plot saved successfully to {output_filename}")
            except Exception as e:
                print(f"Error saving plot to {output_filename}: {e}")
        else:
            plt.show()

    except FileNotFoundError as e:
        print(f"Error: File not found - {e}", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc() # Print detailed error info for debugging


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate and plot average per-residue RMSF from trajectory."
    )
    parser.add_argument(
        "-ref", "--reference", required=True, help="Path to the reference structure file (.pdb, .gro, etc.) for topology."
    )
    parser.add_argument(
        "-traj", "--trajectory", required=True, help="Path to the trajectory file (.xtc, .dcd, .trr, etc.)."
    )
    parser.add_argument(
        "-start_res", type=int, default=1, help="Starting residue number (inclusive) for RMSF calculation. Default is 1."
    )
    parser.add_argument(
        "-end_res", type=int, default=76, help="Ending residue number (inclusive) for RMSF calculation. Default is 76."
    )
    parser.add_argument(
        "-atom_selection", type=str, default="backbone",
        help="MDAnalysis atom selection string *within each residue* (e.g., 'backbone', 'all', 'name CA'). Default is 'backbone'."
    )
    parser.add_argument(
        "-o", "--output", help="Output filename to save the plot (e.g., per_residue_rmsf.png). If not specified, the plot will be displayed."
    )

    args = parser.parse_args()

    # Basic validation for residue numbers
    if args.start_res > args.end_res:
        print("Error: start_res cannot be greater than end_res.", file=sys.stderr)
        sys.exit(1)
    if args.start_res < 0 or args.end_res < 0:
         print("Error: Residue numbers must be non-negative.", file=sys.stderr)
         sys.exit(1)


    calculate_and_plot_per_residue_rmsf(
        args.reference,
        args.trajectory,
        args.start_res,
        args.end_res,
        args.atom_selection,
        args.output
    )

