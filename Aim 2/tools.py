'''
CorBus Tools: some specific-purpose utility functions.
'''
from __future__ import absolute_import, print_function

__author__ = "dstokes"
__copyright__ = "Copyright 2019 Mide Technology Corporation"

import serial.tools.list_ports

from corbus.hardware.interfaces import SerialInterface


def getSerialInterfaces(knownPorts=None, exclusive=True, ftdi=False, all=False,
                        **kwargs):
    """ Find known belt interfaces (i.e. their serial ports). Since the
        enumeration process may assign different paths on different machines
        and/or at different times (based on what's currently attached, etc.),
        the high-level system uses user-specified descriptive names which map
        to specific ports/devices.

        Note that this does not (currently) check whether there are any CorBus
        devices on the port, just that the port exists.

        In the future, this might return other types of stream (TCP, etc.).

        :param knownPorts: a dictionary mapping platform-specific
            USB-to-serial chip serial numbers or explicit serial port names
            (\*NIX ``/dev/tty*`` paths or Windows COM port numbers) to
            user-defined names (e.g. ``display1``, ``bus0``, etc.).
        :keyword exclusive: If `True` and `knownPorts` is not `None`, only
            interfaces in `knownPorts` will be returned. If `False`, all serial
            devices identified as CorBus interfaces will be returned. Will
            always be `False` if `knownPorts` is `None`.
        :keyword ftdi: If `True`, return all found FTDI USB-to-Serial devices
            (e.g. older CorBus USB interfaces), even if `exclusive` and their
            serial numbers do not appear in `knownPorts`. Has no effect if
            `exclusive` is `True` and `knownPorts` is not `None`.
        :keyword all: If `True`, return all found serial ports. Overrides
            `exclusive`.
        :returns: A list of :class:`corbus.hardware.interfaces.SerialInterface` instances. The objects
            will be instantiated in a closed state, and will need to be
            explicitly opened with calls to their `open()` methods before use.
        
        Additional keyword arguments are passed to :class:`SerialInterface`.
    """
    exclusive = False if (all or knownPorts is None) else exclusive
    knownPorts = knownPorts if knownPorts is not None else {}
    ports = []
    
    for p in serial.tools.list_ports.comports():
        portId = None
        sn = str(p.serial_number) if p.serial_number is not None else None

        if all:
            portId = sn or p.device
        elif sn in knownPorts and sn is not None:
            portId = sn
        elif p.device in knownPorts:
            portId = p.device 
        elif exclusive:
            continue
        elif p.vid == 0x0403 and (ftdi or sn.startswith("CORBUS")):
            portId = sn
        else:
            continue
        
        label = knownPorts.get(portId)
        ports.append(SerialInterface(port=p.device, label=label, uid=portId,
                                     info=p, open=False, **kwargs))
        
    return ports