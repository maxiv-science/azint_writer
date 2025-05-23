# NeXus azint writer

**azint-writer** is a Python package for writing HDF5 files in the **NXazint1d** or **NXazint2d** format, a new extension to the [NeXus](https://www.nexusformat.org/) standard. The NXazint1d (NXazint2d) format is specifically designed for storing data related to azimuthal integration in diffraction experiments.

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

```python
import azint, azint_writer, h5py
import numpy as np

h5name = "lab6-lab6_6_master.h5"
h = h5py.File(h5name, 'r')
img = h['/entry/data/data_000001'][0]
poni = '78.poni'
msk = np.zeros(img.shape, dtype=np.uint8)  # a NumPy array

config = {
    'poni': poni,
    'mask': msk,
    'radial_bins': 500,
    'azimuth_bins': 180,
    'n_splitting': 2,
    'error_model': 'poisson',
    'solid_angle': True,
    'polarization_factor': 0.99997,
    'normalized': True,
    'unit': '2th',
}

ai = azint.AzimuthalIntegrator(**config)
data = ai.integrate(img)

init_writer_config = {
    'ai': ai,
    'output_file': 'nxazin_2D.h5',
    'write_1d': False,
    'write_2d': True,
    'instrument_name': 'DanMAX',
    'source_name': 'MAX IV',
    'source_type': 'Synchrotron X-ray Source',
    'source_probe': 'x-ray',
}

nx = azint_writer.NXWriter(**init_writer_config)
nx.add_data(data)
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
