#!/usr/bin/python

import sys

def main():
    # print command line arguments
    amp = float(sys.argv[1])
    freq = float(sys.argv[2])
    duration=float(sys.argv[3])
    tactileHardware(amp, freq, duration)


def tactileHardware(amplitude, frequency,duration):

# Import classes to set up tactor
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
        mags=[.1,.2,.3,.4,.5,.6,.7,.8,.9,1]

# Present stim level
        for x in range(10):
            amplitude=mags[x]
            belt.broadcast.buzz(amp=amplitude, freq=frequency, dur=duration, waveform=Tactor.WAVE_SINE)
            sleep(.05)
if __name__ == "__main__":
    main()