from YXE3D.api.configuration_api import ConfigurationAPI
from YXE3D.api.calibration_api import CalibrationAPI
from YXE3D.api.print_api import PrintAPI, PrintQueueAPI
from YXE3D.api.test_print_api import TestPrintAPI
from YXE3D.api.firmware_api import FirmwareAPI
from YXE3D.infrastructure.configuration_manager import CircutSourcedConfigurationManager


class PrinterAPI(object):
    def __init__(self, ):
        self._configuration_manager = CircutSourcedConfigurationManager()
        self._configuration_api = ConfigurationAPI(self._configuration_manager)
        self._test_print_api = None
        self._firmware_api = None

    def load_printer(self):
        '''Loads a connected printer'''

        self._configuration_api.load_printer()

    def current_printer(self):
        return self._configuration_api.current_printer()

    def get_print_api(self, start_height=0.0):
        return PrintAPI(self._configuration_api.get_current_config(), start_height=start_height)

    def get_print_queue_api(self):
        return PrintQueueAPI(self._configuration_api.get_current_config())

    def get_calibration_api(self, ):
        return CalibrationAPI(self._configuration_manager)

    def get_firmware_api(self):
        if not self._firmware_api:
            self._firmware_api = FirmwareAPI()
        return self._firmware_api

    def get_configuration_api(self):
        return self._configuration_api

    def get_current_config(self):
        return self._configuration_api.get_current_config()

    def get_test_print_api(self, ):
        if not self._test_print_api:
            self._test_print_api = TestPrintAPI()
        return self._test_print_api
