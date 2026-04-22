# Hologram Analysis and Creation App
 
This MATLAB-based application is designed for analyzing and creating holograms using experimentally acquired measurement data, such as photos captured from a CCD camera. The app provides robust tools to manage hologram data and perform essential transformations and analysis using a modular architecture.
 
## Features
 
- **Data Management**:
  - Efficient storage and retrieval of hologram data.
  - Support for various data formats, including CCD camera images.
  
- **Hologram Transformation**:
  - Implementation of key transformation algorithms (e.g., Fourier transforms, filtering).
  - Tools for image preprocessing and reconstruction.
 
- **Visualization**:
  - Interactive visualization of holograms and transformations.
  - Ability to compare input data with processed results.
 
- **Modular Design**:
  - Includes a dedicated `QCI_Model` class to manage data and methods for hologram processing.
  - Flexible and extendable for incorporating new algorithms or data types.
 
## Requirements
 
- MATLAB R2025b or newer
- Signal Processing Toolbox
- Image Processing Toolbox
- Optional: Parallel Computing Toolbox (for enhanced performance)
 
## Installation
 
1. Clone this repository to your local machine:
   ```bash
   git clone https://github.com/QCI-LAB/QCI_HoLLoApp.git
   ```
 
## Important Note for v1.0.0 Users
 
In version **v1.0.0**, the application relies on MATLAB's `pwd` function to resolve internal paths at runtime. `pwd` returns the folder that is **currently open in MATLAB**, not the location of the script itself. As a result, you **must open MATLAB's working directory to the project root folder** before running `main.m`, otherwise the application will fail to locate required classes and functions (e.g., `QCI_Model`) and throw an error such as:
 
```
Error using main (line 5) Invalid default value for property 'QCIModel'
in class 'HoloApp': Unrecognized function or variable 'QCI_Model'.
```
 
**To run the app correctly in v1.0.0:**
1. Open MATLAB.
2. In the MATLAB file browser, navigate to the root folder of the cloned repository (e.g., `QCI_HoLLoApp/`).
3. Make sure this folder is set as the current working directory (it should be visible in the address bar at the top of the MATLAB window).
4. Run `main.m`.
 
This limitation is resolved in **v1.1.0** and later, where `fileparts(mfilename('fullpath'))` is used instead, allowing the app to be launched from any working directory.