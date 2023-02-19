"""
BrainStim Tactile Hardware 

STIM LEVELS (in PEST fucntion): 


"""

def tactileHardware(amplitude, frequency)

# Import classes to set up tactor
from time import sleep

from corbus import Bus, Tactor, tools

# This attempts to automatically find the first CorBus USB interface. You can
# skip this if you already know the serial port's name, but be aware that the
# enumeration of USB serial devices can change between reboots.
try:
    interfaces = tools.getSerialInterfaces()
    port = interfaces[0]
    print("Found interface on %s" % port.port)
except IndexError:
    print("Could not find any CorBus USB interfaces!")
    exit(1)

# Create the `Bus`, the high-level representation of a CorBus array
belt = Bus(port)

# Present stim level

belt.broadcast.buzz(amp=amplitude, freq=frequency, dur=0.15, waveform=Tactor.WAVE_TRIANGLE)

