YXE3D Tools
==================

Note
---------------------------
This is an API only. To run the full suite for the YXE3D printer, use YXE3D at https://github.com/Createcafe3d/YXE3D


Known Issues
--------------------------

Calibration can be poor in some circumstances
* New code is being tried on the calibration branch

Low level errors are not raised properly on usb disconnect


Contributing 
--------------------------

Yes please. 

The YXE3D is continuing the philosophy of it's parent - The Peachy Printer. Its software is community driven, just send us a pull request :)

In order to maintain chaos, please try and ensure:
+ Test Driven Design (TDD) write your tests first then write the code to make them work.
+ Respect the Single Responsibility Principal
+ Follow Onion Architecture Principals
+ PEP8 formatting [Excpeting line length(E501) we are not programming on terminals anymore]


Licence
---------------------------

Please see the LICENSE file


Development 
--------------------------
#### Dependancies

+ python 2.7
+ numpy
+ mock (development only)
+ pyserial
+ python protobuf
+ protobuf
+ libusb1

##### Windows
c++ compiler for python http://www.microsoft.com/en-us/download/details.aspx?id=44266


#### Runing the tests

Run Suite Once

**python test/test-all.py**

Run Suite on Every Change (linux like OS only )

**./runsuite test/test-all.py**



Software Contributers
--------------------------
**CreateCafe Team**
+ Will Topping
+ 

**Original Peachy Team**
+ James Townley (Lead)
+ https://github.com/Pete5746218
