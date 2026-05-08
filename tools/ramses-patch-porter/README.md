# RAMSES-CPP AI Patch Porter ✨

This tool automatically translates legacy RAMSES-2025 Fortran patches (`.f90`) into modern C++17 override classes for RAMSES-CPP.

## Installation
1. Install Python dependencies:
   ```bash
   pip install -r tools/ramses-patch-porter/requirements.txt
   ```
2. Set your preferred AI API key:
   ```bash
   export GEMINI_API_KEY="your-key-here"
   # OR
   export OPENAI_API_KEY="your-key-here"
   ```

## Usage
Run the script providing the input Fortran file and the target solver:
```bash
python tools/ramses-patch-porter/port_patch.py legacy/patch/cooling_fine.f90 --solver Cooling
```

The tool will generate a C++ file in `patch/`.

## Compilation
Simply re-run CMake from your build directory. It will automatically detect the new file and compile it into the simulation:
```bash
cd build
cmake ..
make ramses_3d
```

---
*This tool was architected with the help of Gemini-chan ✨*
