%% read_WERA_raw.m
function [WERA,t,tc,I,Q,Fs]=read_WERA_raw(fname)
% [WERA,t,tc,I,Q,Fs]=read_WERA_raw(fname)
%
% Read into Matlab a WERA RAW (RAW/CAL) binary file created by the WERA
% supplied Fortran science package. The WERA system was developed 
% by Dr. Klaus-Werner Gurgel at the University of Hamburg and is commercially 
% available by Helzel-Messtechnic GmbH.
%
% Due to the nature of some WERA systems sometimes cutting off data midway 
% into a RAW file (depending on the specific WERA system and collection 
% time), this program will exit if the file is incomplete and will display 
% the last chirp number that had data
%
%% Input
%  fname = filename of the RAW or CAL file (including extension)
%          examples: fname='20193221553_gtn.RAW' or
%                    fname='20193251714_csw.CAL'
%    
%% Output: 
%  WERA             = All header information describing data collection
%  t(SAMPZ)         = Time base (in s) for the start of each chirp 
%  tc(MT)           = Time base (in s) for each chirp
%  Fs               = sampling frequency (in Hz)
%  I(ANT,SAMPZ,MT)  = I signal for antenna #ANT, chirp number #SAMPZ, of length MT 
%  Q(ANT,SAMPZ,MT)  = Q signal for antenna #ANT, chirp number #SAMPZ, of length MT
% 
%
%% Uses    
%  read_WERA_header.m
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
% 
% 
%%
echo off
if nargin<1
    return
end
%
% -------------------- Open/ Start Reading the data file  -----------------
if ((strcmpi(fname(end-2:end),'RAW')==1) || (strcmpi(fname(end-2:end),'CAL')==1))
    filename  = fname;
    s         = dir(filename);
    if isempty(s)==1
        error([fname, ' is not found, check path and/or syntax, Aborting'])
    end 
    filebytes = s.bytes;
else
    error([fname, ' is not a RAW type file file, Aborting'])
end
%
if filebytes>512
    fid = fopen(filename,'rb','ieee-le.l64');  % Open the file as a binary
    head = fread(fid,512,'char');
else
    error([filename, ' size is less than 512 bytes, Aborting'])
end
%--------------------------------------------------------------------------
%                          Read the Header Information
%--------------------------------------------------------------------------
[WERA] = read_WERA_header(head);
%--------------------------------------------------------------------------
%                         Read Data from the File
%--------------------------------------------------------------------------
switch upper(WERA.HDCODE)
%--------------------CAL---------------------------
     case {'FMCASW','FMCALI'}           % RAW or CAL file
        NANT     = WERA.NAnt_SORT;
        SAMPZ    = WERA.SAMPZ;
        MT = WERA.MT;
        WERA.RATE = WERA.T_chirp_or;
        t        = (0:SAMPZ-1)*WERA.RATE;
        tc = (1:MT)*WERA.T_chirp_or/MT;
        Fs = 1/(WERA.T_chirp_or/MT);
        I = nan(NANT,SAMPZ,MT);
        Q = nan(NANT,SAMPZ,MT);
        t2 = nan(SAMPZ,MT);
        for chirp = 1:SAMPZ
            try
                for ANT = 1:NANT
                    T = fread(fid,2*MT,'int16');
                    I(ANT,chirp,:)  = T(1:2:2*MT);
                    Q(ANT,chirp,:)  = T(2:2:2*MT);
                    t2(chirp,:) = (1:MT)*WERA.T_chirp_or/MT + WERA.T_chirp_or*(chirp-1);
                end
            catch
                disp(['no more data after chirp #' num2str(chirp)])
                return
            end
        end
       fclose(fid);
%-------------------Future Use--------------------
    case {'FUTURE'}
     % To add future cases/formats
     fclose(fid);
%-------------------------------------------------    
    otherwise
        disp('Header FMSOCO(SORT) or FMSOCN(RFI) or FMCALI(CAL) was expected but not found');
        fclose(fid);
        return
end
end














