from setuptools import setup, find_packages
from setuptools.command.install import install as _Install
from VERSION import version
import os
import glob

# data_files = [('YXE3D/resources/dll', ["YXE3D/dll/libusb-1.0.dll"])]


setup(
    name='YXE3DtoolsAPI',
    version=version,
    description='Tool Set for calibrating the Peachy Printer and printing models',
    options={},
    url="http://www.YXE3D.com",
    author="Peachy Printer",
    author_email="software+YXE3Dtools@YXE3D.com",
    package_data={'': ['*.dll', 'YXE3D/dependancies/win/amd64/*'],
                  '': ['*.dll', 'YXE3D/dependancies/win/x86/*'],
                  '': ['*.dylib', 'YXE3D/dependancies/mac/amd64/*'],
                  '': ['*.so', 'YXE3D/dependancies/linux/amd64/*'],
                  '': ['*.bin', 'YXE3D/dependancies/firmware/*'],
                  '': ['*.dfu', 'YXE3D/dependancies/firmware/*'],
                  },
    install_requires=[
      'protobuf>=2.6.1',
      'pyserial>=2.7',
      'numpy>=1.9.2',
      'libusb1>=1.3.1',
      'YXE3DFirmwareAPI==0.0.1.63'
    ],
    packages=find_packages(),
    py_modules=['VERSION'],
    include_package_data=True
      )

class install(_Install):
    def run(self):
        super(install, self).run(self)
        print "BADA-BADA-KABONG"