# NeXus azint writer

[`azint-writer`](https://github.com/maxiv-science/azint_writer) is a Python package for writing HDF5 files in the [NXazint1d](https://nxazint-hdf5-nexus-3229ecbd09ba8a773fbbd8beb72cace6216dfd5063e1.gitlab-pages.esrf.fr/classes/contributed_definitions/NXazint1d.html) or [NXazint2d](https://nxazint-hdf5-nexus-3229ecbd09ba8a773fbbd8beb72cace6216dfd5063e1.gitlab-pages.esrf.fr/classes/contributed_definitions/NXazint2d.html) format, a new extension to the [NeXus](https://www.nexusformat.org/) standard. The NXazint1d (NXazint2d) format is specifically designed for storing data related to azimuthal integration in diffraction experiments.

> **Note:** The **azint-writer** package is currently coupled to the [`azint`](https://github.com/maxiv-science/azint) library for performing azimuthal integration.


## Features

- **HDF5 File Writing**: Easily generate NeXus-compliant HDF5 files.
- **azint-writer Format**: Supports the new azint-writer format for azimuthal integration.
- **Customizable Configuration**: Flexible options for writing metadata and experimental data.
- **Seamless Integration**: Compatible with existing NeXus tools and libraries.

## Installation

Install **azint-writer** via `conda`:

```bash
conda install -c maxiv azint-writer
```

## azint-writer Format

- **Experiment Metadata**: Beamline details, sample information, and experimental setup.
- **Azimuthal Integration Data**: One-dimensional and two-dimensional (cake) integration results.
- **Calibration Information**: Detector geometry and calibration parameters.

For more information about the NeXus standard, see the [official documentation](https://www.nexusformat.org/).

## Example Usage

Here is a sample Python snippet using **azint-writer**:
([link](https://zenodo.org/records/15744977?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjQxMTA4MjgzLWQ5ODUtNGE3MS04MGU4LTI4MzgwYzAwNDNlYiIsImRhdGEiOnt9LCJyYW5kb20iOiIzOTI4ZmM1YzRhODgwODI3ZDU0ZGVjYTYxNmViNTg0NyJ9.lt_pXKDYR6t29tGKqcm6huHzvgeqlwc5U5I9TXJ5-LQlL865aGLQE7B6-h6ZS7PLQ7yEGf6M3jV1HacXVpSPiA) to file examples and output)

```python
import azint, azint_writer, h5py
import numpy as np

# ------------------------------
# Load detector image and setup
# ------------------------------

# Path to the input HDF5 file containing detector data
h5name = "scan-1737_pilatus.h5"
h = h5py.File(h5name, 'r')

# Extract a single image frame from the dataset
img = h['/entry/instrument/pilatus/data'][10]

# Path to the detector calibration file
poni = 'Si_135mm.poni'

# Path to the mask file for bad pixels or beamstop
mask = 'hot_px_bs_mask.npy'

# ------------------------------
# Azimuthal integration settings
# ------------------------------

config = {
    'poni': poni,                      # Detector geometry calibration file
    'mask': mask,                      # Mask to ignore hot/dead pixels
    'radial_bins': 3000,              # Number of radial bins
    'azimuth_bins': 180,              # Number of azimuthal bins (set to None for 1D only)
    'n_splitting': 21,                # Pixel subdivision for integration precision
    'error_model': 'poisson',         # Error model for propagation
    'solid_angle': True,              # Apply solid angle correction
    'polarization_factor': 0.965,     # Correction for polarization effects
    'normalized': True,               # Normalize intensity values
    'unit': '2th',                    # Output units (e.g., '2th' for 2θ)
}

# ------------------------------
# Output file and metadata config
# ------------------------------

# Create the azimuthal integrator instance
ai = azint.AzimuthalIntegrator(**config)

# Configuration for writing NeXus-compliant HDF5 output
init_writer_config = {
    'ai': ai,
    'output_file': 'nx_azint1d_azint2d.h5',
    'write_1d': True,
    'write_2d': True,  # Set to False if azimuth_bins is None
    'instrument_name': 'DanMAX',
    'source_name': 'MAX IV',
    'source_type': 'Synchrotron X-ray Source',
    'source_probe': 'X-ray',
}

# Create writer instance
nx = azint_writer.NXWriter(**init_writer_config)

# ------------------------------
# Perform integration and save
# ------------------------------

# Integrate the selected image
data = ai.integrate(img)

# Add the integrated data to the NeXus HDF5 file
nx.add_data(data)

# ------------------------------
# Optionally add monitor data
# ------------------------------
from azint_writer import add_monitor

# A 1D array of monitor counts per image (must match number of images)
monitor = [1, ..., nImg]  # Replace with actual monitor data
add_monitor("nx_azint1d_azint2d.h5", monitor)
```

## Contributing

We welcome and encourage contributions to **azint-writer**! Here's how you can get involved:

- Fork the repository.
- Create a new branch for your feature or bug fix.
- Submit a pull request with a clear and concise description of your changes.

Please make sure your contributions follow the [NeXus standard](https://www.nexusformat.org/) guidelines, especially when adding features related to `azint-writer`.

## License

This project is licensed under the MIT License. See the [LICENSE file on GitHub](https://github.com/maxiv-science/azint_writer/blob/main/LICENSE) for details.

## Support

If you encounter any issues, have questions, or would like to suggest improvements, feel free to open an issue on the project’s [GitHub issue tracker](https://github.com/maxiv-science/azint_writer/issues).
