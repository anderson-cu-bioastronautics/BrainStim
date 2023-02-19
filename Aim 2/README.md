# corbus-python: High-level API for Midé's CorBus&trade; haptics system

> *Note:* Before this goes public, the commit history should be cleared in order to make sure all the SHIVR-TSAS material is inaccessible. See [this thread](https://stackoverflow.com/questions/9683279/make-the-current-commit-the-only-initial-commit-in-a-git-repository) for more information.

## Overview
<!-- TODO: Update URLs -->
[CorBus&trade;](https://piezo.com) is a protocol and low-voltage, two-wire network technology developed by Midé Technology Corporation, intended primarily for use in haptics and other wearable electronics. Specifically, CorBus is used to control Midé's 'smart' piezoelectric tactors.

A CorBus network consists of a single host (e.g. the computer running the software), and one or more client devices (collectively referred to as Tactors, although they may perform other functions). Each Tactor has a unique 64-bit ID, and is individually addressable.

### Structure
Physically, a bussed haptic array is a 'flat' collection of Tactors, with no intrinsic indication of the physical configuration or topology. In addition to providing the control architecture, the high-level API creates a layer of logical organization, describing the array's spatial arrangement and grouping tactors by function and/or purpose.

![The logical structure of a haptic array](documentation/source/_static/structure.png "The logical structure of a haptic array")

  - **Bus:** The 'root' object, representing all tactors attached to the same physical interface. Typically, all connected haptic arrays on a single person comprise a single Bus.

  - **Site:** A logical grouping of tactors. Typically, a Site
    represents a collection of tactors in one garment or other physical
    assembly, but this is only a matter of convention. A Site can
    contain tactors spread across multiple garments; likewise, tactors
    in one garment can be organized into multiple Sites (e.g. 'left arm'
    and 'right arm' of a shirt).

  - **Tactor:** The representation of a single haptic actuator. Tactors are typically members or one or more Sites. Tactor objects handle metadata such as the actuator's position in space. The Tactor's coordinate system, and its position therein, is relative to each Site; a Tactor located on a wristband can have Cartesian coordinates relative to a 'body' Site, and polar coordinates in a 'wrist' site, with no inherent mapping between the two systems. A tactor may also be part of the Bus without being in any specific Site.

The relevant high-level metadata is stored in nonvolatile flash memory
on each Tactor. Tactors can be added, removed, or replaced without
affecting the other tactors in the array. Only new or repositioned
tactors require reconfiguration.

Tactors can be selected for activation in several different ways:
individual tactors can be addressed by UUID; multiple tactors can be
addressed together by Site; and one or more can be selected by spatial
position (e.g. those falling within an arc or cone in a direction
indicated by a vector). Tactors must belong to a Site in order to have spatial coordinates, but any tactor can be addressed directly by its UUID.

## Installation
There are two installation variants: basic, which installs only the components required for core API function, and advanced, which installs the packages required to run the GUI tools (editor, desktop server, etc.).

### Basic Installation (API only)
- Download or clone the repo.
- From the command line (\*NIX shell, PowerShell, etc.), change to the downloaded `corbus-python` directory.
- Build/install: `python setup.py install`

### Advanced Installation (GUI tools, etc.)
- Perform the basic installation.
- Install the additional packages: `pip install -r requirements.txt`
