# compSensPack/__init__.py

# Import all the functions from the package
from .utils import printFormatted
from .sparseSignalGenerator import sparseSignal
from .samplingPhase import compressSignal
from .sparseDictionaries import generate_DCT_dictionary, generate_DWT_basis
from .sparseDictionaries import  compute_independent_columns, check_normalization, compute_coherence
from .reconstructionMethods import SL0
from .measurementMatrix import generate_DBDD_matrix, generate_random_matrix
from .recoveryPhase import non_kron_recovery, kron_recovery
from .evaluation import calculate_snr, plot_signals

# Define __all__ for wildcard imports
__all__ = [
    'printFormatted',
    'sparseSignal',
    'compressSignal',
    'generate_DCT_dictionary',
    'generate_DWT_basis',
    'non_kron_recovery',
    'kron_recovery',
    'generate_DBDD_matrix',
    'generate_random_matrix',
    'compute_independent_columns',
    'check_normalization',
    'compute_coherence',
    'calculate_snr',
    'plot_signals'
]

# Optional metadata for the package
__version__ = '1.0.0'
__author__ = 'RosNaviGator'
