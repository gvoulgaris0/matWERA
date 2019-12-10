# WERA-HF-Radar MATLAB scripts
MATLAB(c) Functions andfor reading files created by the WERA HF Radars manufactured by Helzel Messtechnik GmbH (https://helzel-messtechnik.de). The data files the systems generate are created using the Fortran science software package the WERA's use that was created by Dr. Klaus-Werner Gurgel at the University of Hamburg.

The functions in the package are:

Bragg.m                  - for estimating the Bragg frequency

geog2utm.m               - for converting from geographical to UTM / metric coordinates

read_WERA_asc_cur.m      - for reading the 2-D (u,v)current data from a cur_asc file

read_WERA_crad.m         - for reading the current radial binary files

read_WERA_header.m       - It reads and parses the header of any data file

read_WERA_MTfromSORT.m   - for reading the value MT from the WERA sorted (SORT/RFI) binary file

read_WERA_sort.m         - for reading the SORT/RFI files

read_WERA_spec.m         - for reading the SPEC files containing the Doppler spectra estimates

read_WERA_raw.m          - for reading the RAW / CAL files; used if you want to do your own analysis

time2werafile.m          - convert a date / time string to the format used for naming the WERA files

WGS84v.m                 - compute distance and angles between points on the WGS-84 ellipsoidal Earth to within a few millimeters of                                  accuracy using Vincenty's algorithm.
