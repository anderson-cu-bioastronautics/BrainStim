"""
Example 1: Buzzing all tactors.

Buzzing all tactors can be done by broadcasting the command. If the commands
are only being broadcast, very little initialization is required.

IMPORTANT NOTE: This requires Python 2.7! See `requirements.txt` for additional
dependencies.
"""
import random
from time import sleep

from corbus import Bus, Tactor
import tools

# This attempts to automatically find the first CorBus USB interface. You can
# skip this if you already know the serial port's name, but be aware that the
# enumeration of USB serial devices can change between reboots.
try:
    interfaces = tools.getSerialInterfaces()
    port = interfaces[0]
except IndexError:
    print("Could not find any CorBus USB interfaces!")
    exit(1)

# Create the `Bus`, the high-level representation of a CorBus array
belt = Bus(port)

# Do five one-second buzzes, with one second in between
mylist = []
for i in range(1):
    print("Buzz %d" % (i+1))
    # The buzz command. Parameters are:
    #   amp: Buzz amplitude, a normalized percentage (1.0 = max).
    #   freq: Buzz frequency (Hz).
    #   dur: Duration of buzz, in seconds.
    #   waveform: ID of the waveform shape to use. One of:
    #       * Tactor.WAVE_SINE (0)
    #       * Tactor.WAVE_SQUARE (1)
    #       * Tactor.WAVE_TRIANGLE (2)
    #belt.broadcast.buzz(amp=0.15, freq=3, dur=0.4, waveform=Tactor.WAVE_SINE)
    belt.broadcast.buzz(amp=float(0.0003), freq=400, dur=0.2, waveform=Tactor.WAVE_SINE)
    # Buzzes are asynchronous (the `buzz()` method returns immediately), so
    # any delay needs to be done explicitly.
    sleep(2)
print (mylist)