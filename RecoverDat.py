#!/usr/bin/env python2.7

"""
This is a modified version of ExtractDat that does not use the scan
index to find the scans. Instead, is searches through the file looking
for scan headers. The heuristic for identifying a scan header is 
based on looking at a few examples and is not definitive. This code
should only be used if ExtractDat cannot process a file because of 
problems with the scan index.
"""

"""
Copyright (c) 2014 Dr. Philip Wenig
Copyright (c) 2015-2020 John H. Hartman

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version
2.1, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License version 2.1 along with this program.
If not, see <http://www.gnu.org/licenses/>.
"""
import struct
from pprint import *
import sys
import os
import optparse
import glob
import datetime
import math
from collections import defaultdict

VERSION = '2.4'

HDR_NUM_FIELDS = 85
HDR_INDEX_OFFSET = 33
HDR_INDEX_LEN = 39
HDR_TIMESTAMP = 40

SCAN_NUMBER = 9
SCAN_DELTA = 7
SCAN_ACF = 12
SCAN_PREV_TIME = 18
SCAN_TIME = 19
SCAN_EDAC = 31
SCAN_FCF = 34

def Debug(msg):
    if options.debug:
        print msg

KEY_EOS = 0xF # end of scan/acquisition
KEY_EOM = 0x8 # end of mass
KEY_BSCAN = 0xC # B-scan
KEY_B = 0xB # ??
KEY_VOLT = 0x4 # accelerating voltage
KEY_TIME = 0x3 # channel time
KEY_MASS = 0x2 # magnet mass
KEY_DATA = 0x1 # data

DATA_ANALOG = 0x0
DATA_PULSE = 0x1
DATA_FARADAY = 0x8

class DATException(Exception):
    """Base exception class for this module."""
    pass

class EOS(DATException):
    """End of scan."""
    pass

class NotOpen(DATException):
    """Dat file is not open."""
    pass

class UnknownKey(DATException):
    """Unknown tag ("key") in record."""
    pass

class UnknownDataType(DATException):
    """Unknown data type in record."""
    pass

class InvalidScanHeader(DATException):
    """ Scan header corrupted."""
    pass

def _CheckOpen(fd):
    if fd is None:
        raise NotOpen()

def _Seek(fd, offset, whence = 0):
    fd.seek(offset, whence)
    Debug("Offset is now %d" % fd.tell())

class Mass(object):
    """Class for one mass in a scan."""
    def __init__(self, scan, fd):
        """Decodes the mass information from the scan at the specified file and offset"""
        _CheckOpen(fd)
        self.fd = fd
        self.offset = self.fd.tell()
        Debug("Getting mass at offset %d" % self.offset)
        self.scan = scan
        self.magnetMass = None
        self.acceleratingVoltage = None
        self.channelTime = None
        self.duration = None
        self.measurements = defaultdict(list)
        while True:
            Debug("Item at offset %d" % fd.tell())
            tmp = struct.unpack("<I", fd.read(4))[0]
            key = (tmp & 0xF0000000) >> 28
            Debug("Key 0x%x" % key)
            value = tmp & 0x0FFFFFFF
            if key == KEY_EOS: # end of scan
                Debug("EOS")
                raise EOS()
            elif key == KEY_EOM: # end of mass
                Debug("EOM")
                self._SetAttr('duration', value)
                self.size = fd.tell() - self.offset
                break
            elif key == KEY_BSCAN:
                pass # not sure what to do with this
            elif key == KEY_B:
                pass # not sure what to do with this either
            elif key == KEY_VOLT:
                value = scan.edac * 1000.0 / value / 2**18 # TODO: verify this formula
                self.acceleratingVoltage = value  
            elif key == KEY_TIME:
                self._SetAttr('channelTime', value)
            elif key == KEY_MASS:
                self.magnetMass = value * 1.0 / 2**18
            elif key == KEY_DATA:
                flag = (value & 0x0F000000) >> 24
                dataType = (value & 0x00F00000) >> 20
                exp = (value & 0x000F0000) >> 16
                value = value & 0x0000FFFF
                if dataType == DATA_ANALOG:
                    value = value << exp
                    if flag != 0:
                        value = -value
                    self.measurements['analog'].append(value)
                elif dataType == DATA_PULSE:
                    value = value << exp
                    if flag != 0:
                        value = -value
                    self.measurements['pulse'].append(value)
                elif dataType == DATA_FARADAY:
                    value = value << exp
                    if flag != 0:
                        value = -value
                    self.measurements['faraday'].append(value)
                else:
                    raise UnknownDataType(str(dataType))
            else:
                raise UnknownKey(str(key))

    def _SetAttr(self, name, value):
        """Sets an attribute if it is not already set."""
        if getattr(self, name) is not None:
            raise Exception(name + " is already set")
        else:
            setattr(self, name, value)

class Scan(object):
    """Class for one scan in a dat file."""
    def __init__(self, dat, number):
        """Decodes the scan information at the current offset."""
        Debug("Looking for scan header at offset %d" % dat.fd.tell())
        self.headerSize = 47 * 4
        self.offset = dat.fd.tell()
        vals = struct.unpack("<%dI" % (self.headerSize / 4), dat.fd.read(self.headerSize))
        if list(vals[3:6]) != [0xd, 0xe, 0xf]:
            Debug("InvalidScanHeader")
            raise InvalidScanHeader
        Debug("Scan Header")
        for i, val in enumerate(vals):
            Debug("%d %d: 0x%x %d" % (i, i*4, val, val))
        self.number = vals[SCAN_NUMBER]
        if self.number != number:
            Debug("scan number mismatch %d != %d" % (self.number, number))
            raise InvalidScanHeader
        self.delta = vals[SCAN_DELTA]
        self.acf = vals[SCAN_ACF]
        self.time = vals[SCAN_TIME]
        self.fcf = vals[SCAN_FCF]
        self.edac = vals[SCAN_EDAC]
        self._vals = vals
        self.fd = dat.fd
        self.dat = dat

    def __iter__(self):
        """Enables iterating over the masses in a scan."""
        return ScanIterator(self)

    def GetMass(self):
        """Creates a Mass object from the data at the current offset."""
        _CheckOpen(self.fd)
        try:
            mass = Mass(self, self.fd)
        except EOS:
            mass = None
        return mass

class ScanIterator(object):
    """Iterates over the masses in a scan."""
    def __init__(self, scan):
        _CheckOpen(scan.fd)
        self._scan = scan

    def next(self):
        mass = self._scan.GetMass()
        if mass is None:
            raise StopIteration
        return mass

class DatFile(object):
    """Class for one dat file."""
    def __init__(self, path):
        self.path = path
        self.fd = None
        # Read the header.
        Debug("Opening " + self.path)
        with open(self.path, 'rb') as fd:
            _Seek(fd, 0x10)
            fields = HDR_NUM_FIELDS
            tmp = fd.read(fields * 4)
            vals = struct.unpack('<%dI' % fields, tmp)
            Debug("Header");
            for i, val in enumerate(vals):
                Debug("%d: 0x%x %d" % (i, val, val))
            self.timestamp = vals[HDR_TIMESTAMP]
            Debug("timestamp %d" % self.timestamp)
            self._indexLen = vals[HDR_INDEX_LEN]
            Debug("index length %d" % self._indexLen)
            self._indexOffset = vals[HDR_INDEX_OFFSET]
            Debug("index offset %d" % self._indexOffset)
            self._vals = vals
            self.endOfHeader = fd.tell()
            Debug("end of header at offset %d" % self.endOfHeader)
            # Read the offsets although we don't use them, helpful for debugging.
            _Seek(fd, self._indexOffset + 4)
            bytes = self._indexLen * 4
            self._offsets = struct.unpack("<%dI" % self._indexLen, fd.read(bytes))
            Debug("Offsets")
            for i, offset in enumerate(self._offsets):
                Debug("%d: %d" % (i, offset))

    def __iter__(self):
        """Enables iteration over the scans in a dat file."""
        return DatFileIterator(self)

    def Open(self):
        """Opens the dat file."""
        self.fd = open(self.path, 'rb')
        #skip over header
        self.fd.seek(self.endOfHeader)

    def Close(self):
        """Closes the dat file."""
        self.fd.close()
        self.fd = None

    def GetScan(self, index):
        """Creates a Scan object from the data at the current offset."""
        _CheckOpen(self.fd)
        found = False
        while found is False:
            try:
                offset = self.fd.tell()
                scan = Scan(self, index)
                found = True
            except InvalidScanHeader:
                _Seek(self.fd, offset+1, 0) # move forward 1 byte
            except:
                return None
        return scan
    
class DatFileIterator(object):
    """Iterates over the scans in a dat file."""
    def __init__(self, dat):
        self._dat = dat
        self._index = 1

    def next(self):
        scan = self._dat.GetScan(self._index)
        if scan is None:
            raise StopIteration
        self._index += 1
        return scan

description = \
"""\
Decodes the specified dat files and produces a CSV file of their contents. If a single file is
specified the output file has the same base name with a ".csv" suffix. If multiple files are
specified then in addition to producing an output file for each input file an aggregate output
file is produced containing the output from all of the input files. This aggregate output file
has the same base name as the first file with "combinedXX.csv" appended, where 'XX' is a sequence 
number to avoid overwriting existing output files. If a directory is specified then all dat files 
in the directory are processed.
"""
def main(args):
    """Stand-alone application"""
    global options
    global magic

    usage = "Usage %prog [options] file [file ...]"

    parser = optparse.OptionParser(version="%prog " + VERSION, usage=usage, description=description)
    parser.add_option("-c", "--comments",
                      action="store_true", dest="comments",
                      default=False,
                      help="add diagnostic comments to the output file")
    parser.add_option("-d", "--debug",
                      action="store_true", dest="debug",
                      default=False,
                      help="print debugging messages")
    parser.add_option("-D", "--dump",
                      action="store_true", dest="dump",
                      default=False,
                      help="dump values")
    (options, args) = parser.parse_args(args[1:])

    parser.print_version()
    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    # Expand any directories into the dat files they contain. Also keep a list of directories.

    dirs = []
    files = []
    for arg in args:
        if os.path.isdir(arg):
            dirs.append(arg)
            files.extend(glob.glob(os.path.join(arg, '*.dat')))
        else:
            files.append(arg)

    # Sort the dat files by their creation time. 

    dats = [DatFile(f) for f in files]
    dats = sorted(dats, key=lambda x: getattr(x, 'timestamp'))

    # Determine the output directory.

    outputdir = os.path.expanduser('~/Desktop')
    if len(dirs) == 1:
        outputdir = dirs[0]
    elif len(dirs) == 0:
        # No directory specified, try to infer it from the input files.
        dirs = set()
        for dat in dats:
            d = os.path.split(dat.path)[0]
            dirs.add(d)
        dirs = list(dirs)
        if len(dirs) == 1:
            outputdir = dirs[0]
    combinedOutput = None
    if len(dats) > 1:
        # create combined output file name from first dat file name
        base = os.path.splitext(os.path.split(dats[0].path)[1])[0] + 'combined'
        i = 0
        while True:
            path = os.path.join(outputdir, base + '%02d' % i + '.csv')
            if not os.path.exists(path):
                break
            i += 1
        try:
            combinedOutput = open(path, "w")
            print "Writing to", path
        except:
            combinedOutput = None

    first = True
    for dat in dats:
        headers = ["Scan", "Time", "ACF"]
        outputfile = os.path.splitext(dat.path)[0] + '.csv'
        with open(outputfile, "w") as output:
            print "Writing to", outputfile
            dat.Open()   # TODO: context
            # Read elements from FIN2 file if it exists.
            try:
                name = os.path.splitext(dat.path)[0] + '.FIN2'
                with open(name, 'rb') as fin2:
                    for i in xrange(0, 8):
                        line = fin2.readline().strip()
                    elements = line.split(',')[1:]
            except:
                elements = None
            if options.comments:
                print >> output, dat.path, dat.timestamp, datetime.datetime.fromtimestamp(dat.timestamp)
                if combinedOutput != None:
                    print >> combinedOutput, dat.path, dat.timestamp, datetime.datetime.fromtimestamp(dat.timestamp)

            for i, scan in enumerate(dat):
                timestamp = dat.timestamp + scan.time / 1000.0
                results = [str(i+1), '%f' % timestamp, '%f' % scan.acf]
                values = []
                valueHeaders = []
                faraday = False
                
                try:
                    for j, mass in enumerate(scan):
                        modes = ['pulse', 'analog']
                        Debug("scan %d mass %d pulse %d analog %d faraday %d" % (i, j, 
                              len(mass.measurements['pulse']), 
                              len(mass.measurements['analog']), 
                              len(mass.measurements['faraday'])))
                        if len(mass.measurements['faraday']) > 0:
                            modes.append('faraday')
                            faraday = True
                        if headers is not None:
                            element = "Mass%02d" % (j+1) if elements is None else elements[j]
                            if headers is not None:
                                for t in modes:
                                    valueHeaders += ["%s%s" % (element, t[0])] * len(mass.measurements[t])
                                valueHeaders.append('')
                        for t in modes:
                            values += map(lambda x: str(x) if not str(x).startswith('-') else str(-x)+'*', mass.measurements[t])
                        values.append('')
                    if faraday:
                        if headers is not None:
                            headers.append('FCF')
                        results.append('%f' % scan.fcf)

                except UnknownDataType, e:
                    print >> sys.stderr, "Warning: unknown data type 0x%x" % int(e.message)
                    continue
                except UnknownKey, e:
                    print >> sys.stderr, "Warning: unknown key 0x%x" % int(e.message)
                    continue
                if headers is not None:
                    msg = ",".join(headers + valueHeaders)
                    print >> output, msg
                    if combinedOutput != None and first:
                        print >> combinedOutput, msg
                    headers = None
                msg = ",".join(results + values)
                print >> output, msg
                if combinedOutput != None:
                    print >> combinedOutput, msg
            dat.Close()
            first = False

if __name__ == '__main__':
    main(sys.argv)

