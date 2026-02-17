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
