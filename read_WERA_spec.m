%% read_WERA_spec.m
function [Time,LON,LAT,X,Y,freq,fbragg,PXY]=read_WERA_spec(fname,GEO)
% [Time,LON,LAT,X,Y,freq,fbragg,PXY]=read_WERA_spec(fname,[GEO])
% 
% Read WERA binary spec files created with the Fortran science package used
% in the WERA HF radar systems developed by Dr. Klaus-Werner Gurgel at the 
% University of Hamburg and commercially available by Helzel-Messtechnic GmbH.
% 
% use examples:
%
% [time,lon,lat,x,y,freq,fbragg,PXY]=read_WERA_spec(filename,'UTM');
% [time,lon,lat,x,y,freq,fbragg,PXY]=read_WERA_spec(filename);
%
%% Inputs
%  filename = the name of the spec file you want to read (including extension SPEC)   
%  GEO      = Optional string. if ='UTM', then the grid is also converted 
%             to UTM coordinate system.
%     
%% Output
%  Time:               Time in MATLAB datenum format at the middle of burst
%  freq:               Spectral frequency (in Hz)
%  fbragg:             Bragg Frequency for this system (in Hz)
%  LON,LAT:            if SPEC3: Grid coordinates in Long and Lat
%                      if SPEC2: Rx array location coordinates in Long and Lat
%  X,Y:                if SPEC3: Orthogonal coordinate system with 0,0 being
%                                at the WERA station location (in meters) 
%                      if SPEC2: Empty, as SPEC2 does not support this
%  Pxy{PX,PY}(1:freq): Cell Array of Spectral Energy (in db with an approximate 
%                      scaling applied to correct for FFT size, ADC etc. 
%                      (for more details see inside program).
%
%% Intermediate Variables:
%  PX, PY:             Grid positions in integer numbers
%        
%% Uses    
%  read_WERA_header.m
%  read_WERA_MTfromSORT.m
%  geog2utm.m
%        
%% Copyright 2019, George Voulgaris, University of South Carolina
%
% This file is part of matWERA.
%
% matWERA is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% If you find an error please notify G. Voulgaris at gvoulgaris@geol.sc.edu
%%
echo off
%
%       01/22/2011 - Time at the center of the burst is estimated & output.
%
% ------- Parameters Used for Scaling the Spectrum ------------------------
%            scale = 2*(Po*MT*MAX_GAIN*NSER)^2
%
Po = 2^15;     % Reference level based on ADC bits 
MT = 1200;     % This parameter is approximate. The exact valye can be found only
               % on the header of the SORT file but not the SPEC file
MAX_GAIN=7.9;  % Approximate Value the correct value is the max of the 
               % wera_gain.asc file
%NSER  value is set from the header.
% -------------------- Open/ Start Reading the data file  -----------------
if strcmpi(fname(end-3:end),'spec')==1  
    filename  = fname;
    s         = dir(filename);
    if isempty(s)==1
        error([fname, ' is not found, check path and/or syntax, Aborting'])
    end 
    filebytes = s.bytes;
else
    error([fname, ' is not a SPEC type file file, Aborting'])
end
%
if filebytes>512
    fid  = fopen(filename,'rb','ieee-le.l64');  % Open the file as a binary
    head = fread(fid,512,'char');
else
    error([filename, ' size is less than 512 bytes, Aborting'])
end
% -------------------  Read the Header Information ------------------------
[WERA] = read_WERA_header(head);
%-------------------------- frequency calcs -------------------------------
 fbragg     = WERA.fbragg;   % Bragg Frequency (from header)
 dt         = WERA.T_chirp;  % Sampling time interval (in s)
 fs         = 1./dt;         % Sampling Frequency (in Hz)
 nftl       = WERA.NFTL;     % Number of points in FFT (from header)
 df         =  fs/nftl;      % Doppler spectrum freq. resolution
 freq_start = -fs/2; 
 freq_end   = +fs/2;
 freq       = (freq_start+df):df:freq_end;  % Doppler spectrum freq. array
%-------------------------------------------------------------------------
 NSER  = WERA.NSER;
 NX    = WERA.NX;
 NY    = WERA.NY;
 DGT   = WERA.DGT;
 LONo  = WERA.LONo;
 LATo  = WERA.LATo;
 GRAD  = 180.0/pi;
 NCOV  = WERA.N_COVERED;
 MT    = read_WERA_MTfromSORT(fname);
 SCALE = 2*(Po*NSER*MT*MAX_GAIN)^2;       %New
 %SCALE = 1;
%------ Estimate time stamp at the middle of the data collection period ---
dt           = 0.5*NSER*WERA.T_chirp;     % Middle of data in secs
WERA.TIME(3) = WERA.TIME(3)+dt;           % Add the extra secs
Time         = datenum([WERA.DATE WERA.TIME]);
%--------------------------------------------------------------------------
%                         Read Data from the File
%--------------------------------------------------------------------------
switch upper(WERA.HDCODE)
%--------------------------------------------------------------------------
%                         Read pointers SPEC2
%--------------------------------------------------------------------------
    case{'FMSPC2'}
        blocks =floor(NCOV/64)+1;
        PX    = zeros(blocks*64);        % Pre-allocation for speed
        PY    = PX;
        for icnt = 1:blocks  
            list_locs = fread(fid,128,'int32');
            PX((icnt-1)*64+(1:64)) = list_locs(1:2:128);
            PY((icnt-1)*64+(1:64)) = list_locs(2:2:128);
        end
        data = fread(fid,NCOV*nftl,'float32');
        fclose(fid);
        Pxy  = reshape(data,nftl,NCOV);
        clear data
        PXY  = cell(NX,NY);
        for i=1:NCOV                      %blocks*64
            disp(i)
            PXY{PX(i),PY(i)}=10*log10(Pxy(:,i)/SCALE);
        end
        LAT = LATo;
        LON = LONo;
        X   = 0;
        Y   = 0;
%--------------------------------------------------------------------------
%                         Read pointers SPEC3
%--------------------------------------------------------------------------
    case{'FMSPC3'}
        blocks = (floor(NCOV/32)+1);
        PX    = zeros(blocks*32);        % Pre-allocation for speed
        PY    = PX;
        LAT   = PX;
        LON   = PX;
        for icnt = 1:blocks   
            list_locs =fread(fid,128,'int32');
            PX((icnt-1)*32+(1:32))  = list_locs(1:4:128);
            PY((icnt-1)*32+(1:32))  = list_locs(2:4:128);
            LAT((icnt-1)*32+(1:32)) = list_locs(3:4:128);
            LON((icnt-1)*32+(1:32)) = list_locs(4:4:128);
        end
        PX   = PX(1:NCOV);
        PY   = PY(1:NCOV);
        LAT  = LAT(1:NCOV)*1e-7;
        LON  = LON(1:NCOV)*1e-7;
        data = fread(fid,NCOV*nftl,'float32'); 
        fclose(fid);
        Pxy  = reshape(data,nftl,NCOV);
        PXY  = cell(NX,NY);
        for i = 1:NCOV
            PXY{PX(i),PY(i)} = 10*log10(Pxy(:,i)/SCALE);
        end
        %---- Geographical to orthogonal UTM coordinate transformation ----
        if strcmp(GEO,'UTM')==1
            [X,Y]      = geog2utm(LON,LAT,LONo,LATo);  % Geographic 2 Orthogonal
            X          = X*1000;                       % Convert from km to m
            Y          = Y*1000;                       % Convert from km to m
        else
            X          = 0;
            Y          = 0;
        end
%--------------------------------------------------------------------------
        case {'FUTURE'}
        % To add future cases/formats
        fclose(fid);
%--------------------------------------------------------------------------        
    otherwise
    disp('Header FMSPC2 or FMSPC3 was expected but not found');
    fclose(fid);
    return
end
end
