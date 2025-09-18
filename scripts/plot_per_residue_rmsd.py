import numpy as np
import matplotlib.pyplot as plt
import argparse
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
import os
import sys

def calculate_and_plot_per_residue_rmsd(ref_file, traj_file, start_res, end_res, atom_selection, output_filename=None):
    """
    Calculates average per-residue RMSD over a trajectory and plots the results as a bar chart.

    Args:
        ref_file (str): Path to the reference structure file (.pdb, .gro, etc.).
        traj_file (str): Path to the trajectory file (.xtc, .dcd, .trr, etc.).
        start_res (int): The starting residue number (inclusive) for calculation.
        end_res (int): The ending residue number (inclusive) for calculation.
        atom_selection (str): MDAnalysis selection string for atoms within each residue
                              (e.g., 'backbone', 'all', 'name CA').
        output_filename (str, optional): Path to save the output plot image.
                                         If None, the plot is displayed.
    """
    try:
        # Load the universe
        print(f"Loading reference file: {ref_file}")
        print(f"Loading trajectory file: {traj_file}")
        u = mda.Universe(ref_file, traj_file)
        ref = mda.Universe(ref_file) # Load reference separately for stable comparison

        # Define selection for alignment (e.g., protein backbone)
        # This selection determines which atoms are used to superimpose frames
        # to remove overall translation/rotation. Choose a stable region.
        alignment_selection = 'protein and backbone'
        print(f"Using '{alignment_selection}' for trajectory alignment.")

        # Perform trajectory alignment
        # This modifies the coordinates of the trajectory in place
        try:
            print("Aligning trajectory...")
            aligner = align.AlignTraj(u, ref, select=alignment_selection, in_memory=True).run()
            print("Alignment complete.")
        except Exception as e:
             print(f"Error during trajectory alignment: {e}")
             print("Alignment is crucial for meaningful per-residue RMSD.")
             return


        # Prepare to store per-residue RMSD values for each frame
        num_residues = end_res - start_res + 1
        # Initialize a list to store lists of RMSDs for each frame
        all_residue_rmsd_per_frame = []

        print(f"Calculating {atom_selection} RMSD for residues {start_res} to {end_res}...")

        # Iterate through the trajectory frames
        # Note: Alignment was already run, so u.trajectory now iterates over aligned frames
        for ts in u.trajectory:
            current_frame_rmsd = []
            # Iterate through the specified residues
            for res_id in range(start_res, end_res + 1):
                # Select the atoms for the current residue in the current frame (after alignment)
                try:
                    current_res_selection = u.select_atoms(f"resid {res_id} and {atom_selection}")
                    # Select the corresponding atoms in the reference structure
                    ref_res_selection = ref.select_atoms(f"resid {res_id} and {atom_selection}")

                    if len(current_res_selection) == 0 or len(ref_res_selection) == 0:
                         print(f"Warning: No atoms found for residue {res_id} with selection '{atom_selection}'. Skipping.")
                         # Append a placeholder (e.g., NaN or 0) or skip. Let's append NaN.
                         current_frame_rmsd.append(np.nan)
                         continue # Skip to next residue
                    if len(current_res_selection) != len(ref_res_selection):
                         print(f"Warning: Mismatch in atom count for residue {res_id} ({len(current_res_selection)} vs {len(ref_res_selection)}). Skipping.")
                         current_frame_rmsd.append(np.nan)
                         continue # Skip to next residue

                    # Calculate RMSD between the two selections
                    # rmsd function takes coordinates directly
                    res_rmsd = rms.rmsd(current_res_selection.positions, ref_res_selection.positions, center=False, superposition=False)
                    current_frame_rmsd.append(res_rmsd)

                except Exception as e:
                    print(f"Error processing residue {res_id} at frame {ts.frame}: {e}. Skipping.")
                    current_frame_rmsd.append(np.nan) # Append NaN if calculation fails

            all_residue_rmsd_per_frame.append(current_frame_rmsd)

            # Optional: Print progress
            if ts.frame % 100 == 0:
                print(f"Processed frame {ts.frame}/{u.trajectory.n_frames}")

        print("Finished calculating RMSD for all frames.")

        # Convert the list of lists to a numpy array
        # Shape will be (n_frames, n_residues)
        rmsd_data = np.array(all_residue_rmsd_per_frame)

        # Calculate the average RMSD for each residue (average over frames)
        # Ignore NaNs if any residues were skipped
        average_per_residue_rmsd = np.nanmean(rmsd_data, axis=0)

        # Handle cases where a residue had no valid calculations
        residue_numbers = np.arange(start_res, end_res + 1)
        valid_residue_numbers = residue_numbers[~np.isnan(average_per_residue_rmsd)]
        valid_average_rmsd = average_per_residue_rmsd[~np.isnan(average_per_residue_rmsd)]


        # --- Plotting ---
        print("Generating plot...")
        plt.figure(figsize=(15, 7)) # Adjust figure size as needed

        # Create the bar plot
        # Use valid_residue_numbers for x positions and valid_average_rmsd for heights
        bars = plt.bar(valid_residue_numbers, valid_average_rmsd, width=0.8, color='skyblue')

        plt.xlabel("Residue Number")
        plt.ylabel(f"Average RMSD ({atom_selection}) ($\AA$)") # Add units, usually Angstroms
        plt.title(f"Average Per-Residue {atom_selection} RMSD (Residues {start_res}-{end_res})")
        plt.grid(axis='y', linestyle='--', alpha=0.7) # Add horizontal grid lines

        # Improve x-axis ticks visibility if many residues
        if num_residues <= 30:
             plt.xticks(residue_numbers) # Show every residue number
        elif num_residues <= 80:
             plt.xticks(residue_numbers[::2]) # Show every 2nd residue number
        else:
             plt.xticks(residue_numbers[::5]) # Show every 5th residue number
        plt.xlim(start_res - 0.5, end_res + 0.5) # Set x-axis limits to better frame the bars

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
        description="Calculate and plot average per-residue RMSD from trajectory."
    )
    parser.add_argument(
        "-ref", "--reference", required=True, help="Path to the reference structure file (.pdb, .gro, etc.)"
    )
    parser.add_argument(
        "-traj", "--trajectory", required=True, help="Path to the trajectory file (.xtc, .dcd, .trr, etc.)"
    )
    parser.add_argument(
        "-start_res", type=int, default=1, help="Starting residue number (inclusive) for RMSD calculation. Default is 1."
    )
    parser.add_argument(
        "-end_res", type=int, default=76, help="Ending residue number (inclusive) for RMSD calculation. Default is 76."
    )
    parser.add_argument(
        "-atom_selection", type=str, default="backbone",
        help="MDAnalysis atom selection string within each residue (e.g., 'backbone', 'all', 'name CA'). Default is 'backbone'."
    )
    parser.add_argument(
        "-o", "--output", help="Output filename to save the plot (e.g., per_residue_rmsd.png). If not specified, the plot will be displayed."
    )

    args = parser.parse_args()

    # Basic validation for residue numbers
    if args.start_res > args.end_res:
        print("Error: start_res cannot be greater than end_res.", file=sys.stderr)
        sys.exit(1)

    calculate_and_plot_per_residue_rmsd(
        args.reference,
        args.trajectory,
        args.start_res,
        args.end_res,
        args.atom_selection,
        args.output
    )

