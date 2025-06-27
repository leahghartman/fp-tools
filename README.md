# Welcome to FP-Tools!

`fp-tools` is an open-source toolkit for analyzing molecular dynamics (MD) simulations using first-passage time diagnostics. It is designed to provide a fast, flexible, and extensible analysis framework capable of working with large-scale trajectory data from formats such as LAMMPS, VASP, and XYZ.

## Quick Start

### Requirements

- **C++17 compatible compiler** (e.g., GCC 7+, Clang 7+, MSVC 2017+)
- **Python 3** (for plotting and analysis scripts)  
- **yaml-cpp** library for YAML parsing 

  > *Note:* If you do not have yaml-cpp installed, the CMake build system will automatically download and build it for you.

- **Python libraries** (e.g., `matplotlib`, `numpy`) for plotting  

  > *You can install Python dependencies via:*  
  > - **Recommended:** Create and activate a Conda or Mamba environment, then run:  
  >   ```bash
  >   pip install -r requirements.txt
  >   ```  
  > - Alternatively, install directly using pip:  
  >   ```bash
  >   pip install -r requirements.txt
  >   ```

### Build and Install

If you have not downloaded ```fp-tools``` yet, please clone it from GitHub via

```bash
git clone https://github.com/leahghartman/fp-tools.git $HOME/src/fp-tools # or whatever path you prefer
```

From the base of the ```fp-tools``` source directory, execute:

```bash
# Configure project and dependencies
cmake -S . -B build

# Build project
cmake --build build
```

### Add the Executable to Your PATH

To run ```fp-tools``` from anywhere in your terminal, add the build output directory (where the executable resides) to your PATH environment variable.

For example, if the executable is in ```$HOME/src/fp-tools/build/bin```, you can add this to your shell profile:

```bash
export PATH="$HOME/src/fp-tools/build/bin:$PATH"
```

Add this line to your ```~/.bashrc```, ```~/.zshrc```, or appropriate shell configuration file and reload your terminal or source the file:

```bash
source ~/.bashrc
```

### Basic Usage

Run the tool with a YAML configuration file:

```bash
fptools config.yaml
```
If you have **not added** the executable to your PATH, run it by specifying the full or relative path to the executable. For example, if your executable is located in the ```build/bin``` directory inside the project root, run:

```bash
./build/bin/fptools config.yaml
```

or with a full absolute path:

```bash
$HOME/src/fp-tools/build/bin/fptools config.yaml
```

---

## Input/Output Formats





