import unittest
import numpy
import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

from test_helpers import TestHelpers
from infrastructure.audio_disseminator import AudioDisseminator
from domain.laser_control import LaserControl


class AudioDisseminatorTests(unittest.TestCase, TestHelpers):
    def setUp(self):
        self.sample_rate = 1000
        self.on_frequency = self.sample_rate / 4
        self.off_frequency = self.sample_rate / 8
        self.offset = [0, 0]
        self._MODULATION_AMPLITUDE_RATIO = 0.25
        self._SOURCE_AMPLITUDE_RATIO = 1.0 - self._MODULATION_AMPLITUDE_RATIO
        self.laser_control = LaserControl()

    def test_when_laser_off_modulate_it_at_off_frequency(self):
        deseminator = AudioDisseminator(self.laser_control, self.sample_rate, self.on_frequency, self.off_frequency, self.offset)
        self.laser_control.set_laser_off()
        sample_data_chunk = numpy.array([(0, 0)])
        po1 = math.cos(0.0 / 8.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po2 = math.cos(1.0 / 8.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po3 = math.cos(2.0 / 8.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po4 = math.cos(3.0 / 8.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po5 = math.cos(4.0 / 8.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po6 = math.cos(5.0 / 8.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po7 = math.cos(6.0 / 8.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po8 = math.cos(7.0 / 8.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        expected_data = numpy.array([[po1, po1], [po2, po2], [po3, po3], [po4, po4], [po5, po5], [po6, po6], [po7, po7], [po8, po8]])

        actual_data = deseminator.modulate(sample_data_chunk).next()

        self.assertNumpyArrayEquals(expected_data, actual_data)

    def test_when_laser_off_modulate_it_at_off_frequency_with_offset(self):
        offset = [0.1, 0.1]
        deseminator = AudioDisseminator(self.laser_control, self.sample_rate, self.on_frequency, self.off_frequency, offset)
        self.laser_control.set_laser_off()
        sample_data_chunk = numpy.array([(0, 0)])
        po1 = math.cos(0.0 / 8.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + (0.1 * self._SOURCE_AMPLITUDE_RATIO))
        po2 = math.cos(1.0 / 8.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + (0.1 * self._SOURCE_AMPLITUDE_RATIO))
        po3 = math.cos(2.0 / 8.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + (0.1 * self._SOURCE_AMPLITUDE_RATIO))
        po4 = math.cos(3.0 / 8.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + (0.1 * self._SOURCE_AMPLITUDE_RATIO))
        po5 = math.cos(4.0 / 8.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + (0.1 * self._SOURCE_AMPLITUDE_RATIO))
        po6 = math.cos(5.0 / 8.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + (0.1 * self._SOURCE_AMPLITUDE_RATIO))
        po7 = math.cos(6.0 / 8.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + (0.1 * self._SOURCE_AMPLITUDE_RATIO))
        po8 = math.cos(7.0 / 8.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + (0.1 * self._SOURCE_AMPLITUDE_RATIO))
        expected_data = numpy.array([[po1, po1], [po2, po2], [po3, po3], [po4, po4], [po5, po5], [po6, po6], [po7, po7], [po8, po8]])

        actual_data = deseminator.modulate(sample_data_chunk).next()

        self.assertNumpyArrayEquals(expected_data, actual_data)

    def test_offset_can_be_changed(self):
        offset = [0.1,  0.1]
        deseminator = AudioDisseminator(self.laser_control,  self.sample_rate, self.on_frequency, self.off_frequency,  self.offset)
        self.laser_control.set_laser_off()
        deseminator.set_offset(offset)
        sample_data_chunk = numpy.array([(0, 0)])
        po1 = math.cos(0.0 / 8.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + (0.1 * self._SOURCE_AMPLITUDE_RATIO))
        po2 = math.cos(1.0 / 8.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + (0.1 * self._SOURCE_AMPLITUDE_RATIO))
        po3 = math.cos(2.0 / 8.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + (0.1 * self._SOURCE_AMPLITUDE_RATIO))
        po4 = math.cos(3.0 / 8.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + (0.1 * self._SOURCE_AMPLITUDE_RATIO))
        po5 = math.cos(4.0 / 8.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + (0.1 * self._SOURCE_AMPLITUDE_RATIO))
        po6 = math.cos(5.0 / 8.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + (0.1 * self._SOURCE_AMPLITUDE_RATIO))
        po7 = math.cos(6.0 / 8.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + (0.1 * self._SOURCE_AMPLITUDE_RATIO))
        po8 = math.cos(7.0 / 8.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + (0.1 * self._SOURCE_AMPLITUDE_RATIO))
        expected_data = numpy.array([[po1, po1], [po2, po2], [po3, po3], [po4, po4], [po5, po5], [po6, po6], [po7, po7], [po8, po8]])

        actual_data = deseminator.modulate(sample_data_chunk).next()

        self.assertNumpyArrayEquals(expected_data, actual_data)

    def test_when_laser_on_modulate_it_at_on_frequency(self):
        deseminator = AudioDisseminator(self.laser_control,  self.sample_rate, self.on_frequency, self.off_frequency,  self.offset)
        self.laser_control.set_laser_on()
        sample_data_chunk = numpy.array([(0, 0)])
        po1 = math.cos(0.0 / 4.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po2 = math.cos(1.0 / 4.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po3 = math.cos(2.0 / 4.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po4 = math.cos(3.0 / 4.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po5 = math.cos(0.0 / 4.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po6 = math.cos(1.0 / 4.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po7 = math.cos(2.0 / 4.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po8 = math.cos(3.0 / 4.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        expected_data = numpy.array([[po1, po1], [po2, po2], [po3, po3], [po4, po4], [po5, po5], [po6, po6], [po7, po7], [po8, po8]])

        actual_data = deseminator.modulate(sample_data_chunk).next()

        self.assertNumpyArrayEquals(expected_data, actual_data)

    def test_when_laser_on_modulate_it_at_on_frequency_without_applying_offset(self):
        offset = [0.1,  0.1]
        deseminator = AudioDisseminator(self.laser_control,  self.sample_rate, self.on_frequency, self.off_frequency,  offset)
        self.laser_control.set_laser_on()
        sample_data_chunk = numpy.array([(0, 0)])
        po1 = math.cos(0.0 / 4.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po2 = math.cos(1.0 / 4.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po3 = math.cos(2.0 / 4.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po4 = math.cos(3.0 / 4.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po5 = math.cos(0.0 / 4.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po6 = math.cos(1.0 / 4.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po7 = math.cos(2.0 / 4.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        po8 = math.cos(3.0 / 4.0 * 2.0 * math.pi) * self._MODULATION_AMPLITUDE_RATIO
        expected_data = numpy.array([[po1, po1], [po2, po2], [po3, po3], [po4, po4], [po5, po5], [po6, po6], [po7, po7], [po8, po8]])

        actual_data = deseminator.modulate(sample_data_chunk).next()

        self.assertNumpyArrayEquals(expected_data, actual_data)

    def test_on_frequency_must_be_an_even_divisor_of_sample_rate(self):
        sample_rate = 1000
        bad_on_frequency = 71
        off_frequency = 125
        with self.assertRaises(Exception):
            AudioDisseminator(self.laser_control,  sample_rate, bad_on_frequency, off_frequency,  self.offset)

    def test_off_frequency_must_be_an_even_divisor_of_sample_rate(self):
        sample_rate = 1000
        on_frequency = 500
        bad_off_frequency = 99
        with self.assertRaises(Exception):
            AudioDisseminator(self.laser_control,  sample_rate, on_frequency, bad_off_frequency,  self.offset)

    def test_number_of_sample_generated_for_on_and_off_should_be_consistant(self):
        sample_rate = 44100
        on_frequency = 11025
        off_frequency = 7350
        sample_data_chunk = numpy.array([(0, 0)])

        deseminator = AudioDisseminator(self.laser_control,  sample_rate, on_frequency, off_frequency,  self.offset)
        self.laser_control.set_laser_on()
        laser_on = len(list(deseminator.modulate(sample_data_chunk)))
        self.laser_control.set_laser_off()
        laser_off = len(list(deseminator.modulate(sample_data_chunk)))

        self.assertEqual(laser_on, laser_off)

    def test_modualtion_should_be_25_percent_of_amplitude(self):
        deseminator = AudioDisseminator(self.laser_control,  self.sample_rate, self.on_frequency, self.off_frequency,  self.offset)
        self.laser_control.set_laser_on()
        sample_data_chunk = numpy.array([(1.0, 1.0)])
        po1 = math.cos(0.0 / 4.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + 0.75)
        po2 = math.cos(1.0 / 4.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + 0.75)
        po3 = math.cos(2.0 / 4.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + 0.75)
        po4 = math.cos(3.0 / 4.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + 0.75)
        po5 = math.cos(0.0 / 4.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + 0.75)
        po6 = math.cos(1.0 / 4.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + 0.75)
        po7 = math.cos(2.0 / 4.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + 0.75)
        po8 = math.cos(3.0 / 4.0 * 2.0 * math.pi) * (self._MODULATION_AMPLITUDE_RATIO + 0.75)
        expected_data = numpy.array([[po1, po1], [po2, po2], [po3, po3], [po4, po4], [po5, po5], [po6, po6], [po7, po7], [po8, po8]])

        actual_data = deseminator.modulate(sample_data_chunk).next()

        self.assertNumpyArrayEquals(expected_data, actual_data)


if __name__ == '__main__':
    unittest.main()
