%% read_WERA_sort.m
function [WERA,t,R,I,Q]=read_WERA_sort(fname)
% [WERA,t,R,I,Q]=read_WERA_sort(fname)
% 
% Read into Matlab a WERA sorted (SORT/RFI) binary file created by the WERA
% supplied Fortran science package. The WERA system was developed 
% by Dr. Klaus-Werner Gurgel at the University of Hamburg and commercially 
% available by Helzel-Messtechnic GmbH.
%
%% Input
%  fname = filename of the sort or rfi file (including extension)
%          examples: fname='20100360051_hat.SORT' or
%                    fname='20100360051_hat.RFI'
%          
%% Output
%  WERA        = All header information describing data collection
%  R(RC)       = Range in km
%  t(N)        = Time base (in s) for I and Q
%  I(ANT,RC,N) = I signal for antenna #ANT, range cell #RC, of length N
%  Q(ANT,RC,N) = Q signal for antenna #ANT, range cell #RC, of length N
% 
%  Cell size is WERA.RHF (in km)
%  dt of time series is WERA.RATE (in s)
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
% (Last checked on  1/21/2011)
% 
%%
echo off
if nargin<1
fname='20100360051_hat.SORT';  % For testing with a single file (change the name)
end
%
% -------------------- Open/ Start Reading the data file  -----------------
if ((strcmpi(fname(end-3:end),'SORT')==1) || (strcmpi(fname(end-2:end),'RFI')==1))
    filename  = fname;
    s         = dir(filename);
    if isempty(s)==1
        error([fname, ' is not found, check path and/or syntax, Aborting'])
    end 
    filebytes = s.bytes;
else
    error([fname, ' is not a SORT or RFI type file file, Aborting'])
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
LONo   = WERA.LON;
LATo   = WERA.LAT;
GRAD   = 180.0/pi;
% --- Estimate time stamp at the middle of the data collection period ---
dt           = 0.5*WERA.NSER*WERA.T_chirp;    % Middle of data in secs
WERA.TIME(3) = WERA.TIME(3);  %+dt;           % Add the extra secs
Time         = datenum([WERA.DATE WERA.TIME]);
%--------------------------------------------------------------------------
%                         Read Data from the File
%--------------------------------------------------------------------------
switch upper(WERA.HDCODE)
%--------------------SORT & RFI files -------------
    case {'FMSOCO','FMSOCN'}    % SORT & RFI files
        NANT     = WERA.NAnt_SORT;
        NRRANGES = WERA.NRRANGES;
        SAMPZ    = WERA.SAMPZ;
        t        = (0:SAMPZ-1)*WERA.RATE;
        R        = (0:NRRANGES-1)*WERA.RHF+0.5*WERA.RHF-(0.001*WERA.RAN_OFF_m);
        I        = zeros(NANT,NRRANGES,SAMPZ);   % Pre-allocation of memory
        Q        = I;
        for RC = 1:NRRANGES
            for ANT = 1:NANT
                T = fread(fid,2*SAMPZ,'real*4');
                I(ANT,RC,:)  = T(1:2:2*SAMPZ);
                Q(ANT,RC,:)  = T(2:2:2*SAMPZ);
            end
        end
       fclose(fid);
%-------------------Future Use--------------------
    case {'FUTURE'}
     % To add future cases/formats
     fclose(fid);
%-----------------------------------__------------    
    otherwise
        disp('Header FMSOCO(SORT) or FMSOCN(RFI) was expected but not found');
        fclose(fid);
        return
end
end
