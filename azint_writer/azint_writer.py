from collections.abc import Iterable
from datetime import datetime
import h5py
import numpy as np
from enum import Enum
import azint
import os
import logging
from . import __version__

DEFAULT_LOG_LEVEL = "INFO"

log_level = os.getenv("LOG_LEVEL", DEFAULT_LOG_LEVEL).upper()
LOG_LEVELS = {
    "CRITICAL": logging.CRITICAL,
    "ERROR": logging.ERROR,
    "WARNING": logging.WARNING,
    "INFO": logging.INFO,
    "DEBUG": logging.DEBUG,
    "NOTSET": logging.NOTSET,
}
numeric_level = LOG_LEVELS.get(log_level, logging.INFO)
logging.basicConfig(
    level=numeric_level,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


class NXWriter:
    """
    NXWriter class for writing azimuthal integration data to NeXus HDF5 files.
    This class handles the creation of the NeXus hierarchy, including metadata
    such as instrument configuration, source details, monochromator properties,
    and integration parameters. It supports both 1D and 2D data formats, allowing
    for flexible data storage and retrieval.

    Attributes:
        - ai: Azimuthal integrator object (should contain integration results and parameters)
        - output_file (str): Path to output HDF5 file
        - write_1d (bool): Whether to include 1D data in the file
        - write_2d (bool): Whether to include 2D data in the file
        - instrument_name (str): Name of the instrument
        - source_name (str): Name of the source
        - source_type (str): Type of the source (e.g., 'Synchrotron')
        - source_probe (str): Type of probe (e.g., 'x-ray')

    Methods:
        - write_header: Creates and writes the NeXus hierarchy for the dataset,
          including metadata such as instrument configuration, source details,
          monochromator properties, and integration parameters.
        - add_data: Adds azimuthal integration data to the HDF5 file under the
          proper NXdata group.
        - write_radial_axis: Writes radial axis information to the specified group.
    """

    def __init__(
        self, 
        ai, 
        output_file,
        write_1d=True, 
        write_2d=True, 
        instrument_name=None, 
        source_name=None, 
        source_type=None, 
        source_probe=None
    ):
        self.ai = ai
        self.output_file = output_file
        self.write_1d = write_1d
        self.write_2d = write_2d
        self.instrument_name = instrument_name
        self.source_name = source_name
        self.source_type = source_type
        self.source_probe = source_probe

        with h5py.File(self.output_file, "w") as fh_w:
            self.fh = fh_w
            if "entry" not in fh_w: # this condition can be omitted, check that
                self.write_header()
        
    def write_header(self):
        """
        Creates and writes the NeXus hierarchy for the dataset,
        including metadata such as instrument configuration, 
        source details, monochromator properties, and integration parameters.

        This method will automatically determine what type(s) of data 
        (1D/2D) are to be written and populate corresponding subentries.
        """
        logging.info(f"writing header started, azint version is {azint.__version__}")

        entry = self.fh.create_group("entry", track_order=True)
        entry.attrs["NX_class"] = "NXentry"
        entry.attrs["default"] = "data" 
        logging.debug("entry is created in the file")


        if not (self.write_1d and self.write_2d):
            logging.info(f"Creating {'NXazint1d' if self.write_1d else 'NXazint2d'}")
            try:
                self.fh.attrs.setdefault("file_name", np.string_(os.path.basename(self.output_file)))
                self.fh.attrs.setdefault("file_time", np.string_(datetime.now().isoformat()))
                self.fh.attrs.setdefault("HDF5_Version", np.string_(h5py.version.hdf5_version))
            except Exception as e:
                logging.debug(f"Could not set global attributes: {e}")
            
            definition = entry.create_dataset("definition", data='NXazint1d' if self.write_1d else 'NXazint2d')       
        
            solid_angle = entry.create_dataset("solid_angle_applied", data=True if self.ai.solid_angle else False)

            polarization = self.ai.polarization_factor if self.ai.polarization_factor is not None else 0
            polarization_applied = entry.create_dataset("polarization_applied", data=True if self.ai.polarization_factor  is not None else False)

            normalization = entry.create_dataset("normalization_applied", data=True if self.ai.normalized else False)
            monitor = entry.create_dataset("monitor_applied", data=False)

            logging.info("solid_angle_applied and polarization_applied data sets are created")
            logging.info(f"solid_angle: {self.ai.solid_angle}")
            logging.info(f"polarization_factor: {self.ai.polarization_factor}")
        

        # Add instrument
        instrument = entry.create_group("instrument", track_order=True)
        instrument.attrs["NX_class"] = "NXinstrument"
        instrument.attrs["default"] = "name" 
        logging.info(f"Instrument: {self.instrument_name}")
        
        instrument["name"] = np.string_(self.instrument_name)

        # Add monochromator
        mono = instrument.create_group("monochromator", track_order=True)
        mono.attrs["NX_class"] = "NXmonochromator"
        mono.attrs["default"] = "energy"  

        # Add source
        source = instrument.create_group("source", track_order=True)
        source.attrs["NX_class"] = "NXsource"
        source.attrs["default"] = "name" 
        source['name'] = self.source_name
        source['type'] = self.source_type
        source['probe'] = self.source_probe

        poni_file = self.ai.poni
        if isinstance(self.ai.poni, str):
            with open(poni_file, "r") as pf:
                try:
                    logging.info(f"Reading poni file ...")
                    ponif = pf.read()
                    wavelength_found = False
                    for line in ponif.splitlines():
                        if line.startswith("Wavelength:"):
                            wlength_str = line.split(":")[1].strip()
                            try:
                                wlength = float(wlength_str)
                                logging.info(f"From poni file: wavelength: {wlength * 1e10} Å")
                                wavelength_found = True
                            except ValueError as e:
                                logging.error(f"Error converting wavelength to float: {e}")
                                wlength = None
                            break
                    if not wavelength_found:
                        logging.error("Wavelength not found in poni file.")
                        wlength = None
                except Exception as e:
                    logging.error(f"Cannot open poni file: {e}")
        elif isinstance(self.ai.poni, dict):
            ponif = "\n".join(f"{key}: {value}" for key, value in self.ai.poni.items())
            wavelength_found = False
            try:
                wlength = float(self.ai.poni['wavelength'])
                logging.info(f"From poni file: wavelength: {wlength * 1e10} Å")
                wavelength_found = True
            except ValueError as e:
                logging.error(f"Error converting wavelength to float: {e}")
                wlength = None
            if not wavelength_found:
                logging.error("Wavelength not found in poni dict.")
                wlength = None
        else:
            logging.error("Provided format for poni is wrong.")

        # Now handle data splitting
        if (self.write_1d and self.write_2d):
            logging.info(f"Creating 1D and 2D data ...")
            if self.ai.azimuth_axis is None:
                logging.error("2d data is not available.")
            # 1D data subentry
            azint1dSE = entry.create_group("azint1d", track_order=True)
            azint1dSE.attrs["NX_class"] = "NXsubentry"

            definition = azint1dSE.create_dataset("definition", data="NXazint1d")
                 
        
            solid_angle = azint1dSE.create_dataset("solid_angle_applied", data=True if self.ai.solid_angle else False)

            polarization = self.ai.polarization_factor if self.ai.polarization_factor is not None else 0
            polarization_applied = azint1dSE.create_dataset("polarization_applied", data=True if self.ai.polarization_factor  is not None else False)

            normalization = azint1dSE.create_dataset("normalization_applied", data=True if self.ai.normalized else False)
            monitor = azint1dSE.create_dataset("monitor_applied", data=False)

            azint1dSE["instrument"] = h5py.SoftLink('/entry/instrument')

            reduction = azint1dSE.create_group("reduction", track_order=True)
            reduction.attrs["NX_class"] = "NXprocess"
            prog = reduction.create_dataset("program", data="azint-pipeline")
            ver = reduction.create_dataset("version", data=f"azint {azint.__version__}\nazint-writer {__version__}")
            reduction.create_dataset("date", data=np.string_(datetime.now().isoformat()))
            ref = reduction.create_dataset("reference", data="Jensen, A. B., et al., (2022). J. Synchrotron Rad. 29, 1420-1428.\nhttps://doi.org/10.1107/S1600577522008232,\nhttps://maxiv-science.github.io/azint/")
            note = reduction.create_dataset("note", data="Geometry convention:\nAzimuthal origin in the horizontal plane to the right of the beam position, i.e., at 3 o’clock,\non the detector face. Positive azimuthal direction: clockwise.")

            input = reduction.create_group("input", track_order=True)
            input.attrs["NX_class"] = "NXparameters"
            dset = input.create_dataset("poni", data=ponif, track_order=True)
            dset.attrs["filename"] = poni_file if isinstance(self.ai.poni, str) else "Poni dict."

            wavelength = mono.create_dataset("wavelength", data=wlength * 1e10, track_order=True)
            wavelength.attrs["units"] = "angstrom"
            energy = mono.create_dataset("energy", data=(1.2398 * 1e-9) / wlength, track_order=True)
            energy.attrs["units"] = "keV"

            wavelength2 = input.create_dataset("wavelength", data=wlength * 1e10)
            wavelength2.attrs["units"] = "angstrom"

            input.create_dataset("n_splitting", data=self.ai.n_splitting)
            input.create_dataset("radial_bins", data=self.ai.radial_bins)
            input.create_dataset("azimuth_bins", data=1)
            input.create_dataset("unit", data=self.ai.unit)    
            input.create_dataset("mask", data=self.ai.mask_path if self.ai.mask_path else ("A numpy array was provided" if self.ai.mask is not None else "None"))
            
            input.create_dataset("solid_angle", data=True if self.ai.solid_angle else False)

            polarization = self.ai.polarization_factor if self.ai.polarization_factor is not None else 0
            input.create_dataset("polarization_factor", data=polarization)
            error_model = self.ai.error_model if self.ai.error_model else "None"
            input.create_dataset("error_model", data=error_model)


            azint1d = azint1dSE.create_group("data")
            azint1d.attrs["NX_class"] = "NXdata"
            azint1d.attrs["signal"] = "I"
            azint1d.attrs["axes"] = [".", "radial_axis"]
            azint1d.attrs["interpretation"] = "spectrum"
            self.write_radial_axis(azint1d, self.ai.unit, self.ai.radial_axis, self.ai.radial_bins)

            azint2dSE = entry.create_group("azint2d", track_order=True)
            azint2dSE.attrs["NX_class"] = "NXsubentry"

            definition = azint2dSE.create_dataset("definition", data="NXazint2d")
                 
        
            solid_angle = azint2dSE.create_dataset("solid_angle_applied", data=True if self.ai.solid_angle else False)

            polarization = self.ai.polarization_factor if self.ai.polarization_factor is not None else 0
            polarization_applied = azint2dSE.create_dataset("polarization_applied", data=True if self.ai.polarization_factor  is not None else False)

            normalization = azint2dSE.create_dataset("normalization_applied", data=True if self.ai.normalized else False)
            monitor = azint2dSE.create_dataset("monitor_applied", data=False)

            azint2dSE["instrument"] = h5py.SoftLink('/entry/instrument')

            reduction = azint2dSE.create_group("reduction", track_order=True)
            reduction.attrs["NX_class"] = "NXprocess"
            prog = reduction.create_dataset("program", data="azint-pipeline")
            ver = reduction.create_dataset("version", data=f"azint {azint.__version__}\nazint-writer {__version__}")
            reduction.create_dataset("date", data=np.string_(datetime.now().isoformat()))
            ref = reduction.create_dataset("reference", data="Jensen, A. B., et al., (2022). J. Synchrotron Rad. 29, 1420-1428.\nhttps://doi.org/10.1107/S1600577522008232,\nhttps://maxiv-science.github.io/azint/")
            note = reduction.create_dataset("note", data="Geometry convention:\nAzimuthal origin in the horizontal plane to the right of the beam position, i.e., at 3 o’clock,\non the detector face. Positive azimuthal direction: clockwise.")

            input = reduction.create_group("input", track_order=True)
            input.attrs["NX_class"] = "NXparameters"

            wavelength2 = input.create_dataset("wavelength", data=wlength * 1e10)
            wavelength2.attrs["units"] = "angstrom"

            input.create_dataset("n_splitting", data=self.ai.n_splitting)
            input.create_dataset("radial_bins", data=self.ai.radial_bins)
            input.create_dataset("azimuth_bins", data=self.ai.azimuth_bins)
            input.create_dataset("unit", data=self.ai.unit)
            input.create_dataset("mask", data=self.ai.mask_path if self.ai.mask_path else ("A numpy array was provided" if self.ai.mask is not None else "None"))
            
            input.create_dataset("solid_angle", data=True if self.ai.solid_angle else False)

            polarization = self.ai.polarization_factor if self.ai.polarization_factor is not None else 0
            input.create_dataset("polarization_factor", data=polarization)
            error_model = self.ai.error_model if self.ai.error_model else "None"
            input.create_dataset("error_model", data=error_model)

            azint2d = azint2dSE.create_group("data")
            azint2d.attrs["NX_class"] = "NXdata"
            azint2d.attrs["signal"] = "I"
            azint2d.attrs["axes"] = [".", "azimuthal_axis", "radial_axis"]
            azint2d.attrs["interpretation"] = "image"
            entry["data"] = h5py.SoftLink('/entry/azint1d/data')
            self.write_radial_axis(azint2d, self.ai.unit, self.ai.radial_axis, self.ai.radial_bins)
            dset = azint2d.create_dataset("azimuthal_axis", data=self.ai.azimuth_axis)
            dset.attrs["units"] = "degrees"
            dset.attrs["long_name"] = "azimuthal bin center"

            if isinstance(self.ai.azimuth_bins, Iterable):
                aedges = self.ai.azimuth_bins
            else:
                acentres = self.ai.azimuth_axis
                awidth = acentres[1] - acentres[0]
                aedges = (acentres - 0.5 * awidth)
                aedges = np.append(aedges, aedges[-1] + awidth)

            dset2 = azint2d.create_dataset("azimuthal_axis_edges", data=aedges)
            dset2.attrs["units"] = "degrees"
            dset2.attrs["long_name"] = "azimuthal bin edges"

            azint1dSE.attrs["default"] = "data"
            azint2dSE.attrs["default"] = "data"
            entry.attrs["default"] = "data"
            
        elif self.write_1d:
            logging.info(f"Creating just 1D data ...")
            # ONLY 1D DATA section
            reduction = entry.create_group("reduction", track_order=True)
            reduction.attrs["NX_class"] = "NXprocess"
            prog = reduction.create_dataset("program", data="azint-pipeline")
            ver = reduction.create_dataset("version", data=f"azint {azint.__version__}\nazint-writer {__version__}")
            reduction.create_dataset("date", data=np.string_(datetime.now().isoformat()))
            ref = reduction.create_dataset("reference", data="Jensen, A. B., et al., (2022). J. Synchrotron Rad. 29, 1420-1428.\nhttps://doi.org/10.1107/S1600577522008232,\nhttps://maxiv-science.github.io/azint/")
            note = reduction.create_dataset("note", data="Geometry convention:\nAzimuthal origin in the horizontal plane to the right of the beam position, i.e., at 3 o’clock,\non the detector face. Positive azimuthal direction: clockwise.")

            input = reduction.create_group("input", track_order=True)
            input.attrs["NX_class"] = "NXparameters"
            dset = input.create_dataset("poni", data=ponif, track_order=True)
            dset.attrs["filename"] = poni_file
            # Add wavelength and energy to mono
            wavelength = mono.create_dataset("wavelength", data=wlength * 1e10, track_order=True)
            wavelength.attrs["units"] = "angstrom"
            energy = mono.create_dataset("energy", data=(1.2398 * 1e-9) / wlength, track_order=True)
            energy.attrs["units"] = "keV"
        
            # Add other info from poni to input
            wavelength2 = input.create_dataset("wavelength", data=wlength * 1e10)     
            wavelength2.attrs["units"] = "angstrom"
            
            input.create_dataset("n_splitting", data=self.ai.n_splitting)
            input.create_dataset("radial_bins", data=self.ai.radial_bins)
            azimuth_bins = self.ai.azimuth_bins if self.ai.azimuth_bins else 1
            input.create_dataset("azimuth_bins", data=azimuth_bins)
            input.create_dataset("unit", data=self.ai.unit)
            input.create_dataset("mask", data=self.ai.mask_path if self.ai.mask_path else ("A numpy array was provided" if self.ai.mask is not None else "None"))
            
            input.create_dataset("solid_angle", data=True if self.ai.solid_angle else False)

            polarization = self.ai.polarization_factor if self.ai.polarization_factor is not None else 0
            input.create_dataset("polarization_factor", data=polarization)
            error_model = self.ai.error_model if self.ai.error_model else "None"
            input.create_dataset("error_model", data=error_model)
            
            azint1d = entry.create_group("data", track_order=True)
            azint1d.attrs["NX_class"] = "NXdata"
            azint1d.attrs["signal"] = "I"
            azint1d.attrs["axes"] = [".", "radial_axis"]
            azint1d.attrs["interpretation"] = "spectrum"
            self.write_radial_axis(azint1d, self.ai.unit, self.ai.radial_axis, self.ai.radial_bins)

            entry.attrs["default"] = "data"
        
        elif self.write_2d:
            logging.info(f"Creating just 2D data ...")
            # ONLY 2D DATA section
            reduction = entry.create_group("reduction", track_order=True)
            reduction.attrs["NX_class"] = "NXprocess"
            prog = reduction.create_dataset("program", data="azint-pipeline")
            ver = reduction.create_dataset("version", data=f"azint {azint.__version__}\nazint-writer {__version__}")
            reduction.create_dataset("date", data=np.string_(datetime.now().isoformat()))
            ref = reduction.create_dataset("reference", data="Jensen, A. B., et al., (2022). J. Synchrotron Rad. 29, 1420-1428.\nhttps://doi.org/10.1107/S1600577522008232,\nhttps://maxiv-science.github.io/azint/")
            note = reduction.create_dataset("note", data="Geometry convention:\nAzimuthal origin in the horizontal plane to the right of the beam position, i.e., at 3 o’clock,\non the detector face. Positive azimuthal direction: clockwise.")

            input = reduction.create_group("input", track_order=True)
            input.attrs["NX_class"] = "NXparameters"
            dset = input.create_dataset("poni", data=ponif, track_order=True)
            dset.attrs["filename"] = poni_file
            # Add wavelength and energy to mono
            wavelength = mono.create_dataset("wavelength", data=wlength * 1e10, track_order=True)
            wavelength.attrs["units"] = "angstrom"
            energy = mono.create_dataset("energy", data=(1.2398 * 1e-9) / wlength, track_order=True)
            energy.attrs["units"] = "keV"
        
            # Add other info from poni to input
            wavelength2 = input.create_dataset("wavelength", data=wlength * 1e10)     
            wavelength2.attrs["units"] = "angstrom"
            
            input.create_dataset("n_splitting", data=self.ai.n_splitting)
            input.create_dataset("radial_bins", data=self.ai.radial_bins)
            azimuth_bins = self.ai.azimuth_bins
            input.create_dataset("azimuth_bins", data=azimuth_bins)
            input.create_dataset("unit", data=self.ai.unit)
            input.create_dataset("mask", data=self.ai.mask_path if self.ai.mask_path else ("A numpy array was provided" if self.ai.mask is not None else "None"))
            
            input.create_dataset("solid_angle", data=True if self.ai.solid_angle else False)

            polarization = self.ai.polarization_factor if self.ai.polarization_factor is not None else 0
            input.create_dataset("polarization_factor", data=polarization)
            error_model = self.ai.error_model if self.ai.error_model else "None"
            input.create_dataset("error_model", data=error_model)
            
            azint2d = entry.create_group("data", track_order=True)
            azint2d.attrs["NX_class"] = "NXdata"
            azint2d.attrs["signal"] = "I"
            azint2d.attrs["axes"] = [".", "azimuthal_axis", "radial_axis"]
            azint2d.attrs["interpretation"] = "image"
            self.write_radial_axis(azint2d, self.ai.unit, self.ai.radial_axis, self.ai.radial_bins)
            dset = azint2d.create_dataset("azimuthal_axis", data=self.ai.azimuth_axis)
            dset.attrs["units"] = "degrees"
            dset.attrs["long_name"] = "azimuthal bin center"

            if isinstance(self.ai.azimuth_bins, Iterable):
                aedges = self.ai.azimuth_bins
            else:
                acentres = self.ai.azimuth_axis
                awidth = acentres[1] - acentres[0]
                aedges = (acentres - 0.5 * awidth)
                aedges = np.append(aedges, aedges[-1] + awidth)

            dset2 = azint2d.create_dataset("azimuthal_axis_edges", data=aedges)
            dset2.attrs["units"] = "degrees"
            dset2.attrs["long_name"] = "azimuthal bin edges"

            entry.attrs["default"] = "data"

    def add_data(self, integrated_data):
        """
        Add azimuthal integration data to the HDF5 file under the proper NXdata group.

        Parameters:
        - integrated_data (tuple): Tuple of form (I, errors_1d, cake, errors_2d), where:
            - I: 1D intensity array
            - errors_1d: 1D error array
            - cake: 2D intensity array
            - errors_2d: 2D error array
        """

        I, errors_1d, cake, errors_2d = integrated_data
        data = {}
        with h5py.File(self.output_file, "r+") as fh_u:
            if (self.write_1d and self.write_2d):
            # if cake is not None: # will have eta bins
                data["/entry/azint1d/data/I"] = I
                data["/entry/azint2d/data/I"] = cake
                if self.ai.normalized:
                    if "/entry/azint2d/data/norm" not in fh_u:
                        data["/entry/azint2d/data/norm"] = self.ai.norm_2d
                    if "/entry/azint1d/data/norm" not in fh_u:
                        data["/entry/azint1d/data/norm"] = self.ai.norm_1d
                if errors_2d is not None:
                    data["/entry/azint2d/data/I_errors"] = errors_2d
                if errors_1d is not None:
                    data["/entry/azint1d/data/I_errors"] = errors_1d
            elif self.write_1d:  # must be radial bins only, no eta, ie 1d.
                data["/entry/data/I"] = I
                if self.ai.normalized:
                    if "/entry/data/norm" not in fh_u:
                        data["/entry/data/norm"] = self.ai.norm_1d
                if errors_1d is not None:
                    data["/entry/data/I_errors"] = errors_1d
            elif self.write_2d:
                data["/entry/data/I"] = cake
                if self.ai.normalized:
                    if "/entry/data/norm" not in fh_u:
                        data["/entry/data/norm"] = self.ai.norm_2d
                if errors_2d is not None:
                    data["/entry/data/I_errors"] = errors_2d
            else:
                logging.error(f"At least one of 1D or 2D should be written")

            for key, value in data.items():
                new_dset = fh_u.get(key)
                if not new_dset:
                    if "norm" in key:
                        new_dset = fh_u.create_dataset(key, data=value, track_order=True)
                    else:
                        new_dset = fh_u.create_dataset(key, dtype=value.dtype,
                                                   shape=(0, *value.shape),
                                                   maxshape=(None, *value.shape),
                                                   chunks=(1, *value.shape))
                    # I and I_error created here
                    new_dset.attrs["units"] = "arbitrary units"
                    new_dset.attrs["long_name"] = "intensity"
                    if "I_error" in key:
                        new_dset.attrs.modify("long_name", "estimated errors on intensity")
                    if "norm" in key:
                        new_dset.attrs.modify("long_name", "effective number of pixels contributing to the corresponding bin")

                if "norm" not in key:
                    n = new_dset.shape[0]
                    new_dset.resize(n + 1, axis=0)
                    new_dset[n] = value

    def write_radial_axis(self, group, unit, radial_axis, radial_bins):
        # real dataset for radial axis is always "radial axis"
        dset = group.create_dataset("radial_axis", data=radial_axis, track_order=True)
        dset.attrs["long_name"] = "q" if unit == "q" else "2theta"
        dset.attrs["units"] = "1/angstrom" if unit == "q" else "degrees"
        
        # Calculate edges
        if isinstance(radial_bins, Iterable):
            edges = radial_bins
        else:
            centres = radial_axis
            width = centres[1]-centres[0]
            edges = (centres-0.5*width)
            edges = np.append(edges,edges[-1]+width)

        dsete = group.create_dataset("radial_axis_edges", data=edges, track_order=True)
        dsete.attrs["long_name"] = "q bin edges" if unit == "q" else "2theta bin edges"
        dsete.attrs["units"] = "1/angstrom" if unit == "q" else "degrees"

def add_monitor(h5file, monitor_data):
    """
    Add or update a 'monitor' group with a 'data' dataset in a NeXus HDF5 file.

    This function checks for the existence of specific groups in the NeXus file
    and inserts monitor data into `/entry/monitor` if `/entry/definition` exists,
    or into `/entry/azint1d/monitor` and `/entry/azint2d/monitor` otherwise.

    If the 'monitor' group or 'data' dataset already exists, they are reused or
    updated. If the shape of the existing dataset differs from the new data, it
    is deleted and recreated.

    Parameters
    ----------
    h5file : str
        Path to the HDF5 (.nxs) file to be modified.
    monitor_data : array-like
        NumPy array or array-like object containing monitor signal data to write.
    """

    if not isinstance(monitor_data, np.ndarray):
        try:
            monitor_data = np.asarray(monitor_data)
        except Exception as e:
            raise TypeError(f"Could not convert monitor_data to a NumPy array: {e}")

    with h5py.File(h5file, "r+") as fh_u:
        if "/entry/definition" in fh_u:
            entry_paths = ["entry"]
        else:
            entry_paths = ["entry/azint1d", "entry/azint2d"]
        n_image = fh_u[f"entry/data/I"].shape[0]
        if monitor_data.shape[0] != n_image:
            raise ValueError(
                f"Monitor data length ({monitor_data.shape[0]}) does not match "
                f"number of images ({n_image})."
            )
        for entry_path in entry_paths:
            try:
                # Get or create the monitor group
                monitor = fh_u.require_group(f"{entry_path}/monitor")
                monitor.attrs["NX_class"] = np.string_("NXmonitor")

                # Check if 'data' exists and update if possible
                if "data" in monitor:
                    dset = monitor["data"]
                    if dset.dtype != monitor_data.dtype or dset.shape != monitor_data.shape:
                        del monitor["data"]
                        monitor.create_dataset("data", data=monitor_data, track_order=True)
                    else:
                        dset[...] = monitor_data
                else:
                    monitor.create_dataset("data", data=monitor_data, track_order=True)
                fh_u[entry_path]["monitor_applied"][...] = True

                logging.info(f"Monitor data added to {entry_path}/monitor")
            except Exception as e:
                logging.error(f"Error handling monitor data at {entry_path}: {e}")
