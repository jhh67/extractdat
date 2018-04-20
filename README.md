# extractdat
Decoder for Thermo Element ICP Mass Spectrometer dat files.

ExtractDat.py is a Python program that can either be used as a library or a stand-alone command. It
requires Python 2.7, so make sure that is installed on your machine. There are two ways to use
ExtractDat.py as a stand-alone command -- drag-and-drop or on the command line. The former is
probably easiest. Put ExtractDat.py on your desktop, then drag-and-drop your DAT file(s) onto it. If
you drag-and-drop a single DAT file ExtractDat.py will produce an output file with the same base
name as the input file and the suffix ".csv". If you drag-and-drop multiple DAT files then in
addition to an output file for each input file ExtractDat.py will produce a single output file with
the same base name as the first input file and ending with "combined.csv". If the input files are
from different folders then the individual output files will appear in the same folder as the
corresponding input files, and the combined output file will appear in the same directory as the
first input file. Finally, if you drag-and-drop a folder ExtractDat.py will process all the DAT
files in that folder as if they all had been drag-and-dropped. If the combined output file already
exists  then ExtractDat.py will append a sequence number to the output file name to avoid
overwriting existing combined output files.

The syntax for running ExtractDat.py on the command line is: ExtractDat.py [options] file [file
...]. Run ExtractDat.py --help for more details.
