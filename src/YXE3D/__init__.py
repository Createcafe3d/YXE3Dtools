# -*- mode: python; basic-offset: 4 -*-
import os
import sys
import logging

logger = logging.getLogger('peachy')

from YXE3D.api.peachy_printer_api import PrinterAPI
from YXE3D.infrastructure.communicator import MissingPrinterException

try:
    from VERSION import version
except:
    version = "DEV"

try:
    from infrastructure.peachyusb import lib_version
except:
    lib_version = "Unknown"