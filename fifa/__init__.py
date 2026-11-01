"""FIFA: read and quality-filter ZTF and ATLAS forced-photometry light curves."""

from .reader import SN_reader_ztf, SN_reader_atlas

__version__ = "0.1.0"
__all__ = ["SN_reader_ztf", "SN_reader_atlas"]
