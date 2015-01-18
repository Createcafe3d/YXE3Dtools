import math
import numpy
import logging
from fractions import gcd
from domain.disseminator import Disseminator


class AudioDisseminator(Disseminator):
    _MODULATION_AMPLITUDE_RATIO = 0.25
    _SOURCE_AMPLITUDE_RATIO = 1.0 - _MODULATION_AMPLITUDE_RATIO

    def __init__(self, laser_control, audio_data_writer, sampling_rate, on_frequency, off_frequency, offset):
        self._laser_control = laser_control
        self._audio_data_writer = audio_data_writer
        self._x_offset, self._y_offset = offset
        logging.info("Laser Control: Modulation On: %s" % on_frequency)
        logging.info("Laser Control: Modulation Off: %s" % off_frequency)
        logging.info("Laser Offset: %2.3f, %2.3f" % (self._x_offset, self._y_offset))
        if sampling_rate % on_frequency != 0:
            raise Exception("The on_frequency must divide evenly into sampling_rate")
        if sampling_rate % off_frequency != 0:
            raise Exception("The off_frequency must divide evenly into sampling_rate")

        off_laser_steps = sampling_rate / off_frequency
        on_laser_steps = sampling_rate / on_frequency
        lcm = self._lcm([off_laser_steps, on_laser_steps])
        self._actual_samples_per_second = sampling_rate / lcm

        self.off_laser_wave = numpy.array(self._get_cos_wave(off_laser_steps, lcm / off_laser_steps))
        self.on_laser_wave = numpy.array(self._get_cos_wave(on_laser_steps, lcm / on_laser_steps))
        logging.info("Started audio modulation with On Frequency: %s and Off Frequency: %s" % (on_frequency, off_frequency))

    @property
    def samples_per_second(self):
        return self._actual_samples_per_second

    def _lcm(self, numbers):
        return reduce(lambda x, y: (x * y)/gcd(x, y), numbers, 1)

    def _get_cos_wave(self, steps, cycles):
        wave = []
        scale = 2.0 * math.pi
        for _ in range(0, cycles):
            for i in range(0, int(steps)):
                cos_wave = math.cos(i * 1.0 / steps * 1.0 * scale)
                wave.append(cos_wave)
        return wave

    def _modulate(self, data):
        if self._laser_control.laser_is_on():
            pattern = self.on_laser_wave
            for (left, right) in data:
                l = numpy.multiply([self._MODULATION_AMPLITUDE_RATIO + (left * self._SOURCE_AMPLITUDE_RATIO)], pattern)
                r = numpy.multiply([self._MODULATION_AMPLITUDE_RATIO + (right * self._SOURCE_AMPLITUDE_RATIO)], pattern)
                yield numpy.column_stack((l, r))
        else:
            pattern = self.off_laser_wave
            for (left, right) in data:
                r = numpy.multiply([self._MODULATION_AMPLITUDE_RATIO + ((right + self._y_offset) * self._SOURCE_AMPLITUDE_RATIO)], pattern)
                l = numpy.multiply([self._MODULATION_AMPLITUDE_RATIO + ((left + self._x_offset) * self._SOURCE_AMPLITUDE_RATIO)], pattern)
                yield numpy.column_stack((l, r))

    def set_offset(self, offset):
        self._x_offset, self._y_offset = offset

    def process(self, data):
        modulated = self._modulate(data)
        self._audio_data_writer.write_chunk(modulated)

    def next_layer(self, height):
        self._audio_data_writer.next_layer(height)

    def close(self):
        self._audio_data_writer.close()
