import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

def extract_quoted_string(line):
    """
    Helper function to extract text within the first pair of double quotes in a string.
    Returns the extracted string or None if not found.
    """
    try:
        start_idx = line.find('"')
        if start_idx == -1:
            return None
        end_idx = line.find('"', start_idx + 1)
        if end_idx == -1:
            return None
        return line[start_idx + 1:end_idx]
    except Exception:
        return None

def read_xvg(filename):
    """
    Reads data and header information from a GROMACS .xvg file.

    Args:
        filename (str): The path to the .xvg file.

    Returns:
        tuple: A tuple containing:
               - x_data (np.ndarray): The data for the x-axis (first column).
               - y_data (np.ndarray): The data for the y-axis (remaining columns).
               - plot_title (str or None): The extracted plot title.
               - x_label (str or None): The extracted x-axis label.
               - y_label (str or None): The extracted y-axis label.
               Returns (None, None, None, None, None) if reading fails or no data.
    """
    plot_title = None
    x_label = None
    y_label = None
    header_lines_count = 0

    # First pass to find header end and extract labels
    try:
        with open(filename, 'r') as f:
            for line in f:
                stripped_line = line.strip()
                if not stripped_line or stripped_line.startswith('#') or stripped_line.startswith('@'):
                    header_lines_count += 1

                    # Attempt to extract metadata
                    if stripped_line.startswith('@ title'):
                        plot_title = extract_quoted_string(stripped_line)
                    elif stripped_line.startswith('@ xaxis label'):
                        x_label = extract_quoted_string(stripped_line)
                    elif stripped_line.startswith('@ yaxis label'):
                         y_label = extract_quoted_string(stripped_line)
                    # Add parsing for @ legend if needed for multiple curves
                    # elif stripped_line.startswith('@ legend'):
                    #    legend_entries.append(extract_quoted_string(stripped_line))

                else:
                    # Found the first non-header line (start of data)
                    break
    except FileNotFoundError:
        print(f"Error: File not found at {filename}")
        return None, None, None, None, None
    except Exception as e:
        print(f"Error reading header of {filename}: {e}")
        return None, None, None, None, None

    # Second pass to load the actual data using numpy
    try:
        # Use unpack=True to get columns as separate arrays if needed,
        # but keeping it as one array is easier for plotting multiple columns.
        # delimiter=None allows loadtxt to auto-detect whitespace (space, tab)
        data = np.loadtxt(filename, skiprows=header_lines_count, delimiter=None)

    except Exception as e:
        print(f"Error reading data from {filename}: {e}")
        # If loadtxt fails, it might mean the file is empty or has format issues after headers
        return None, None, plot_title, x_label, y_label # Return extracted info even if data fails

    # Check if data was loaded and has enough columns
    if data.ndim < 2 or data.shape[1] < 2:
        print(f"Warning: No data or insufficient columns found in {filename} after skipping {header_lines_count} header lines.")
        # If data is 1D, it might be a single column file, which we can't plot X vs Y
        return None, None, plot_title, x_label, y_label

    x_data = data[:, 0]     # First column is X
    y_data = data[:, 1:]    # Remaining columns are Y

    return x_data, y_data, plot_title, x_label, y_label

def plot_xvg_data(x_data, y_data, plot_title, x_label, y_label, output_filename=None):
    """
    Plots the provided data using matplotlib.

    Args:
        x_data (np.ndarray): X-axis data.
        y_data (np.ndarray): Y-axis data (can be multiple columns).
        plot_title (str or None): Title for the plot.
        x_label (str or None): Label for the x-axis.
        y_label (str or None): Label for the y-axis.
        output_filename (str or None): If provided, saves the plot to this file.
                                       Otherwise, displays the plot.
    """
    if x_data is None or y_data is None:
        print("No data to plot.")
        return

    fig, ax = plt.subplots(figsize=(10, 6)) # Create a figure and axes

    # Plot the data. Matplotlib's plot handles multiple y-columns correctly.
    ax.plot(x_data, y_data)

    # Set title and labels if available
    if plot_title:
        ax.set_title(plot_title)
    else:
         ax.set_title(f"Plot of Data") # Default title if none extracted

    if x_label:
        ax.set_xlabel(x_label)
    else:
        ax.set_xlabel("X-axis") # Default label

    if y_label:
        ax.set_ylabel(y_label)
    else:
        ax.set_ylabel("Y-axis") # Default label

    ax.grid(True) # Add a grid for better readability

    # Add a legend if there are multiple Y columns
    if y_data.shape[1] > 1:
         # If you parsed legend entries from headers, you could use them here
         # ax.legend(legend_entries, loc='best')
         # Otherwise, show a default legend
         ax.legend(loc='best') # 'best' chooses the least obstructive location

    # Save or show the plot
    if output_filename:
        try:
            plt.savefig(output_filename, bbox_inches='tight')
            print(f"Plot saved to {output_filename}")
        except Exception as e:
            print(f"Error saving plot to {output_filename}: {e}")
    else:
        plt.show()

    plt.close(fig) # Close the figure after showing or saving

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Read and plot data from a GROMACS .xvg file.")
    parser.add_argument("-f", "--file", required=True, help="Input .xvg file path.")
    parser.add_argument("-o", "--output", help="Output filename to save the plot (e.g., plot.png). If not specified, the plot will be displayed.")

    args = parser.parse_args()

    print(f"Attempting to read {args.file}...")
    x_data, y_data, plot_title, x_label, y_label = read_xvg(args.file)

    if x_data is not None and y_data is not None:
        print("Data read successfully. Plotting...")
        plot_xvg_data(x_data, y_data, plot_title, x_label, y_label, args.output)
    else:
        print("Could not read data or file is empty/malformed.")

