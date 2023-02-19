"""
Example 2: Buzzing specific tactors.

To be accessed individually, the tactors on the bus must be initialized. The
enumeration process identifies tactors and assigns them IDs. The previous
example broadcast to all tactors, so it did not need this step.

Once the Bus has been initialized, individual tactors can be
accessed in two ways:
    * belt.tactors: a list of Tactor instances.
    * belt.tactorsById: a dictionary of Tactors, keyed by their unique ID

IMPORTANT NOTE: This requires Python 2.7! See `requirements.txt` for additional
dependencies.
"""

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

# Initialize the bus. Finds all attached tactors. This may take several
# seconds. Optionally, you can speed it up by providing the method a list of
# known tactor IDs via its keyword argument `knownIds`.
belt.hardwareInit()
print("Found %d tactor(s)" % len(belt.tactors))

# Individual tactors can be accessed in two ways:
#   belt.tactors: a list of Tactor instances.
#   belt.tactorsById: a dictionary of Tactors, keyed by their unique ID
# The order of tactors in belt.tactors is consistent between runs.

# Buzz each tactor:
for t in belt.tactors:
    print("Buzzing tactor %#x" % t.id)

    # See notes on the buzz() command in previous example
    t.buzz(amp=1.0, freq=350, dur=0.5, waveform=Tactor.WAVE_SINE)
    sleep(0.5)
