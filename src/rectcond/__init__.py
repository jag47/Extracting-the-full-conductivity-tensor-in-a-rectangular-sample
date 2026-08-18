"""rectcond: extract the full 2D conductivity tensor in a rectangular sample.

Implements J. Gelfond and O. Vafek, arXiv:2608.06643 (2026).
"""

from rectcond.core import ConductivityTensor, extract_tensor

__version__ = "0.1.0"
__all__ = ["ConductivityTensor", "extract_tensor", "__version__"]
