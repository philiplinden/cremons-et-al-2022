import logging
from pathlib import Path
import sys
import pandas as pd

# where to look for the data if no path is given
DEFAULT_DATA_DIR = Path('../data')
# Temperature of step-heating experiment in degrees Celsius
EXPERIMENT_TEMPERATURE = [650, 700, 750, 800]
# Total water measured from step-heating experiments in ppm
EXPERIMENT_WATER_PPM = [1522, 762, 176, 22]

# ------------------------------------------------------------------------------
#   Logging
# ------------------------------------------------------------------------------
class Colors:
    GREY = "\x1b[37m"
    GREEN = "\x1b[32m"
    YELLOW = "\x1b[33m"
    RED = "\x1b[31m"
    PURPLE = "\x1b[35m"
    BLUE = "\x1b[34m"
    LIGHT_BLUE = "\x1b[36m"
    BLINK_RED = "\x1b[5m\x1b[31m"
    RESET = "\x1b[0m"
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


class CustomFormatter(logging.Formatter):
    """Logging Formatter to add colors and count warning / errors"""
    def __init__(self):
        super(CustomFormatter, self).__init__()
        module = f"{Colors.BLUE}%(module)-6s{Colors.RESET} "
        lvl = f"{Colors.BOLD}%(levelname)+8s:{Colors.END} "
        msg = "%(message)s"
        self.FORMATS = {
            logging.DEBUG:
            module + Colors.PURPLE + lvl + Colors.PURPLE + msg + Colors.RESET,
            logging.INFO:
            module + Colors.GREY + lvl + Colors.GREY + msg + Colors.RESET,
            logging.WARNING:
            module + Colors.YELLOW + lvl + Colors.YELLOW + msg + Colors.RESET,
            logging.ERROR:
            module + Colors.RED + lvl + Colors.RED + msg + Colors.RESET,
            logging.CRITICAL:
            module + Colors.BLINK_RED + lvl + Colors.BLINK_RED + msg +
            Colors.RESET
        }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


console_handler = logging.StreamHandler(sys.stdout)
console_handler.setFormatter(CustomFormatter())
# logfile_handler = logging.handlers.RotatingFileHandler(
#     etl.resolve_path_as(LOG_PATH, 'log'), maxBytes=512)
# logfile_handler.setFormatter(CustomFormatter())
logging.basicConfig(datefmt="%Y-%m-%dT%H:%M:%S",
                    handlers=[
                        # logfile_handler,
                        console_handler,
                    ])
log = logging.getLogger()
log.setLevel(logging.INFO)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#   Importing data from files
# ------------------------------------------------------------------------------
def spectrum_from_file(
    path: Path,
    name='reflectance',
    header=None,
    delimiter=',',
) -> pd.Series:
    spectrum = pd.read_csv(
        path,
        index_col=0,
        header=header,
        delimiter=delimiter,
        names=['wavelength', name],
    )
    return spectrum.squeeze()



def import_spectral_data(data_dir: Path | str = DEFAULT_DATA_DIR):
    if isinstance(data_dir, str):
        data_dir = Path(data_dir)

    # Import reflectance spectra (processed to remove hydration and organics)
    samples = {
        'Mature Mare': data_dir / Path('Mare_70181_Spectra.txt'),
        'Mature Highlands': data_dir / Path('Highlands_62231_Spectra.txt'),
        'Pyroxene':
        data_dir / Path('Apollo15Sample15555ReddishBrownPyroxeneB.txt'),
        'Immature Mare': data_dir / Path('Mare_71061_Spectra.txt'),
        'Immature Highlands': data_dir / Path('Highlands_61221_Spectra.txt'),
    }

    lab_spectra = []
    for name, path in samples.items():
        lab_spectra.append(
            spectrum_from_file(path, name, header=2, delimiter='\t'))

    endmember_spectra = pd.concat(lab_spectra, axis=1)

    experiments = zip(EXPERIMENT_TEMPERATURE, EXPERIMENT_WATER_PPM)

    # Import observations of MORB step-wise heating reflectance spectra
    imported_reflectances = []
    for deg_c, ppm in experiments:
        MORB_spectrum = spectrum_from_file(data_dir /
                                                   Path(f"{deg_c}^oC.csv"),
                                                   name=f'{ppm} ppm')
        MORB_spectrum.index = 1e4 / MORB_spectrum.index  # 1/cm to microns
        MORB_spectrum = MORB_spectrum / 100  # percent to decimal

        # since we converted frequency to wavelength, need to resort data
        MORB_spectrum = MORB_spectrum.sort_index()
        imported_reflectances.append(MORB_spectrum)

    heated_MORB_spectra = pd.concat(imported_reflectances, axis=1)

    # Low wavelength portion of MORB spectrum (<1.5 microns)
    MORB_D38A_LowLam = spectrum_from_file(
        data_dir / Path("Morb_D38A_Low_wavelength.txt"),
        name='MORB D38A low wavelengths',
    )

    return endmember_spectra, heated_MORB_spectra, MORB_D38A_LowLam
# ------------------------------------------------------------------------------
