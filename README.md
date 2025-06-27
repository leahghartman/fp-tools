# fp-tools

A high-performance toolkit for **melting curve prediction** from molecular dynamics (MD) trajectories, leveraging first-passage time analysis. Accelerated with C++/MPI and optimized for LAMMPS/VASP/XYZ inputs.

*A Python/C++ toolkit for molecular dynamics post-processing, specializing in liquid-phase melting diagnostics and general analysis. Supports CPU and GPU workflows.*

---

## Key Features

- Support for LAMMPS, VASP, and XYZ trajectory formats
- Implements advanced analyses including MSD, RDF, FKT, FSKT, MFPT, and melting diagnostics (D(r))
- Efficient C++ backend with optional MPI parallelism and OpenMP acceleration
- Python integration for easy plotting and post-processing
- Modular, extensible design to add new analysis modules or input formats

---

## Quick Start

### Requirements

- **C++17 compatible compiler** (e.g., GCC 7+, Clang 7+, MSVC 2017+)
- **Python 3** (for plotting and analysis scripts)  
- **YAML-CPP** library for configuration parsing  
  > *Note:* If you do not have YAML-CPP installed, the CMake build system will automatically download and build it for you.
- Python libraries (e.g., `matplotlib`, `numpy`) for plotting  
  > *You can install Python dependencies via:*  
  > - **Recommended:** Create and activate a Conda or Mamba environment, then run:  
  >   ```bash
  >   pip install -r requirements.txt
  >   ```  
  > - Alternatively, install directly using pip:  
  >   ```bash
  >   pip install -r requirements.txt
  >   ```

---

## Build and Install

### Clone the repository

```bash
git clone https://github.com/yourusername/fptools.git $HOME/src/fp-tools
cd $HOME/src/fp-tools
