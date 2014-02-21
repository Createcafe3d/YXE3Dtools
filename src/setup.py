from cx_Freeze import setup, Executable
import VERSION

# Dependencies are automatically detected, but it might need
# fine tuning.
buildOptions = dict(packages = [], excludes = [])

import sys
base = 'Win32GUI' if sys.platform=='win32' else None

executables = [
    Executable('peachyprintertools.py', base=base, targetName = 'PeachyPrinterTools')
]

setup(name='Peachy Printer Tools',
      version = version,
      description = 'Tool Set for calibrating the Peachy Printer and printing models',
      options = dict(build_exe = buildOptions),
      executables = executables)