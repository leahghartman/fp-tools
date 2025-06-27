import os
import sys
import matplotlib.pyplot as plt

def main():
    if len(sys.argv) < 2:
        print("Usage: plot_msd.py <msd_data_file> [format]")
        return
    
    msd_file = sys.argv[1]
    output_format = sys.argv[2] if len(sys.argv) > 2 else "png"
    base_name = os.path.splitext(os.path.basename(msd_file))[0]

    # Load data
    lags, msd, counts = [], [], []
    with open(msd_file) as f:
        for line in f:
            if line.startswith("#") : continue
            parts = line.strip().split()
            lags.append(int(parts[0]))
            msd.append(float(parts[1]))
            counts.append(int(parts[2]))

    plt.figure()
    plt.plot(lags, msd, marker='o')
    plt.xlabel("Lag")
    plt.ylabel("MSD")
    plt.title("Mean Squared Displacement")

    # Create 'output/' directory if it doesn't exist
    output_dir = "output"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Output image path
    plot_path = os.path.join(output_dir, "msd.svg")

    # Save with the correct format
    plt.savefig(plot_path)
    print(f"Plot saved to {output_dir}")

if __name__ == "__main__":
    main()
