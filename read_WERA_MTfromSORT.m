%% read_WERA_MTfromSORT.m
function [MT]=read_WERA_MTfromSORT(fname)
% [MT]=read_WERA_MTfromSORT(fname)
%
% Read the value MT from the WERA sorted (SORT/RFI) binary file created
% by the WERA systems developed by Dr. Klaus-Werner Gurgel at the University 
% of Hamburg and commercially available by Helzel-Messtechnic GmbH. This
% value is used in read_WERA_spec.m to scale the spectrum.
%
%% Input
%  fname = filename of the sort or rfi file (including extention)
%          
%% Output
%  MT = MT value to be used to scale the spectrum back to original units
%       if no SORT file is found, then a default value of MT=1200 is returned.
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
%  First version developed on September 10, 2011
%
echo off
% -------------------- Open/ Start Reading the data file  -----------------
filename     = [fname,'.SORT'];
s            = dir(filename);
if isempty(s)==1
   disp([filename, ' SORT file not found, using default (approximate) value of MT=1200'])
   MT = 1200;
   return
end 
   filebytes = s.bytes;
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
MT     = WERA.MT;
end
