# azint-writer

**azint-writer** is a Python package for writing HDF5 files in the **azint-writer** format, a new extension to the [NeXus](https://www.nexusformat.org/) standard. The azint-writer format is specifically designed for storing data related to azimuthal integration in diffraction experiments.

---

## Features

- **HDF5 File Writing**: Easily generate NeXus-compliant HDF5 files.
- **azint-writer Format**: Supports the new azint-writer format for azimuthal integration.
- **Customizable Configuration**: Flexible options for writing metadata and experimental data.
- **Seamless Integration**: Compatible with existing NeXus tools and libraries.

---

## Installation

Install Azint Writer via `conda`:

```bash
conda install -c maxiv azint-writer
```
---

## azint-writer Format

- **Experiment Metadata**: Beamline details, sample information, and experimental setup.
- **Azimuthal Integration Data**: One-dimensional and two-dimensional (cake) integration results.
- **Calibration Information**: Detector geometry and calibration parameters.

For more details on the NeXus standard, see the NeXus documentation.

---

## Contributing
We welcome contributions to Azint Writer! Here's how you can help:
- Fork the repository.
- Create a new branch for your feature or bug fix.
- Submit a pull request with a clear description of your changes.

Please follow the NeXus standard guidelines when adding features related to azint-writer.

---

## License

This project is licensed under the MIT License. See the LICENSE file for details.

---

## Support
If you encounter any issues or have questions, feel free to open an issue on GitHub.

