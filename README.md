# matWERA: A MATLAB package for reading binary data files recorded by a WERA HF Radar

MATLAB(r) functions for reading the binary files created by the WERA HF Radars manufactured by Helzel Messtechnik GmbH (https://helzel-messtechnik.de). The data files the WERA systems generate are created using the Fortran science software package developed by Dr. Klaus-Werner Gurgel at the University of Hamburg and distributed with the WERA systems.

The functions in the package are:

  - Bragg.m                  - for estimating the Bragg frequency

  - geog2utm.m               - for converting from geographical to UTM / metric coordinates

  - read_WERA_asc_cur.m      - for reading the 2-D (u,v)current data from a cur_asc file

  - read_WERA_crad.m         - for reading the current radial binary files

  - read_WERA_header.m       - It reads and parses the header of any data file

  - read_WERA_MTfromSORT.m   - for reading the value MT from the WERA sorted (SORT/RFI) binary file

  - read_WERA_sort.m         - for reading the SORT/RFI files

  - read_WERA_spec.m         - for reading the SPEC files containing the Doppler spectra estimates

  - read_WERA_raw.m          - for reading the RAW / CAL files; used if you want to do your own analysis

  - time2werafile.m          - convert a date / time string to the format used for naming the WERA files

  - WGS84v.m                 - compute distance and angles between points on the WGS-84 ellipsoidal Earth to within a few millimeters of accuracy using Vincenty's algorithm. This function is a vectorized version of the code of Michael Kleder that was downloaded in 2013 and it is  still available at http://www.mathworks.com/matlabcentral/fileexchange/5379.
  
 # Michael Kleder (2013). Geodetic distance on WGS84 earth ellipsoid (https://www.mathworks.com/matlabcentral/fileexchange/5379-geodetic-distance-on-wgs84-earth-ellipsoid), MATLAB Central File Exchange. Retrieved, June 2013.

A few example files are included and the read_examples.m script show how to use them. Not all files are found in here as some files are too large to be uploaded to GitHub.

If you use that code please cite as follows: 

# George Voulgaris, 2019. matWERA (V1.0.0): A collection of MATLAB(r) functions to read Wellen-Radar (WERA) binary files. (Version V1.0.0). Zenodo. http://doi.org/10.5281/zenodo.3570892

George Voulgaris, University of South Carolina, USA.
Email: gvoulgaris@geol.sc.edu
