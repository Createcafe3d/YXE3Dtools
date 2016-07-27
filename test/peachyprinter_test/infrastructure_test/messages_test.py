import unittest
import sys
import os
import logging

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'src'))

from YXE3D.infrastructure.messages import MoveMessage, DripRecordedMessage, SetDripCountMessage, MoveToDripCountMessage, IAmMessage, EnterBootloaderMessage, GetAdcValMessage, ReturnAdcValMessage, PrinterStatusMessage


class MoveMesssageTests(unittest.TestCase):

    def test_move_message_encodes_and_decodes(self):
        inital_message = MoveMessage(77, 88, 55)
        proto_bytes = inital_message.get_bytes()
        self.assertTrue(len(proto_bytes) > 0)
        decoded_message = MoveMessage.from_bytes(proto_bytes)
        self.assertEqual(inital_message, decoded_message)


class DripRecordedMesssageTests(unittest.TestCase):

    def test_move_message_encodes_and_decodes(self):
        inital_message = DripRecordedMessage(77)
        proto_bytes = inital_message.get_bytes()
        self.assertTrue(len(proto_bytes) > 0)
        decoded_message = DripRecordedMessage.from_bytes(proto_bytes)
        self.assertEqual(inital_message, decoded_message)


class SetDripCountMesssageTests(unittest.TestCase):

    def test_move_message_encodes_and_decodes(self):
        inital_message = SetDripCountMessage(77)
        proto_bytes = inital_message.get_bytes()
        self.assertTrue(len(proto_bytes) > 0)
        decoded_message = SetDripCountMessage.from_bytes(proto_bytes)
        self.assertEqual(inital_message, decoded_message)


class MoveToDripCountMesssageTests(unittest.TestCase):

    def test_move_message_encodes_and_decodes(self):
        inital_message = MoveToDripCountMessage(77)
        proto_bytes = inital_message.get_bytes()
        self.assertTrue(len(proto_bytes) > 0)
        decoded_message = MoveToDripCountMessage.from_bytes(proto_bytes)
        self.assertEqual(inital_message, decoded_message)


class IAmMesssageTests(unittest.TestCase):

    def test_move_message_encodes_and_decodes(self):
        inital_message = IAmMessage("77", "88", "99", 9600)
        proto_bytes = inital_message.get_bytes()
        self.assertTrue(len(proto_bytes) > 0)
        decoded_message = IAmMessage.from_bytes(proto_bytes)
        self.assertEqual(type("a"), type(decoded_message.sn))
        self.assertEqual(type("a"), type(decoded_message.hwrev))
        self.assertEqual(type("a"), type(decoded_message.swrev))
        self.assertEqual(type(9600), type(decoded_message.dataRate))
        self.assertEqual(inital_message, decoded_message)


class ReturnAdcValMessageTests(unittest.TestCase):

    def test_returnadcval_message_encodes_and_decodes(self):
        initial_message = ReturnAdcValMessage(3)
        proto_bytes = initial_message.get_bytes()
        self.assertTrue(len(proto_bytes) > 0)
        decoded_message = ReturnAdcValMessage.from_bytes(proto_bytes)
        self.assertEqual(initial_message, decoded_message)


class GetAdcValMessageTests(unittest.TestCase):

    def test_getadcval_message_encodes_and_decodes(self):
        initial_message = GetAdcValMessage(3)
        proto_bytes = initial_message.get_bytes()
        self.assertTrue(len(proto_bytes) > 0)
        decoded_message = GetAdcValMessage.from_bytes(proto_bytes)
        self.assertEqual(initial_message, decoded_message)


class EnterBootloaderTests(unittest.TestCase):
    def test_enter_boot_loader_encodes_and_decodes(self):
        initial_message = EnterBootloaderMessage()
        proto_bytes = initial_message.get_bytes()
        decoded_message = EnterBootloaderMessage.from_bytes(proto_bytes)
        self.assertEqual(initial_message, decoded_message)


class PrinterStatusMesssageTests(unittest.TestCase):

    def test_move_message_encodes_and_decodes(self):
        inital_message = PrinterStatusMessage(True, True, False, False, 32768)
        proto_bytes = inital_message.get_bytes()
        self.assertTrue(len(proto_bytes) > 0)
        decoded_message = PrinterStatusMessage.from_bytes(proto_bytes)
        self.assertEqual(inital_message, decoded_message)


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level='DEBUG')
    unittest.main()
