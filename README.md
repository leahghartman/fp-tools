# Welcome to FP-Tools!

`fp-tools` is an open-source toolkit for analyzing molecular dynamics (MD) simulations using first-passage time diagnostics. It is designed to provide a fast, flexible, and extensible analysis framework capable of working with large-scale trajectory data from formats such as LAMMPS, VASP, and XYZ.

## Table of Contents

- [Quick Start](#quick-start)
    - [Build and Install](#build-and-install)
    - [Add Executable to Path](#add-the-executable-to-your-path)
    - [Basic Usage](#basic-usage)
- [Input/Output Formats](#inputoutput-formats)
- [Contributing](#contributing)
- [License](#license)

## Quick Start

### Requirements

#### C++ Build Environment
- **C++20 compatible compiler** (e.g., GCC 10+, Clang 11+, MSVC 2019+)
- CMake 3.15+ for build configuration
- **Python 3.8+** for plotting and analysis scripts

#### Python Dependencies
These are required to run the plotting scripts:
- `numpy`
- `matplotlib`

  > **Recommended setup:** 
  > Create and activate a virtual environment with Conda or Mamba:
  >   ```bash
  >   conda create -n fptools python=3.10
  >   conda activate fptools
  >   pip install -r requirements.txt
  >   ```  

  > **Or install directly with pip:** 
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





