%% time2werafile
function [fpath, ftime]=time2werafile(year,mo,day,hr,mi)
%% [fpath, ftime]=time2werafile(year,mo,day,hr,mi)
%
% Function to generate character strings with the directory structure 
% and name of files for a given time. Usefull when you want to find a file
% corresponding to a particular time. It utilizes the naming conversion
% used by the WERA systems developed by Dr. Klaus-Werner Gurgel at the 
% University of Hamburg and commercially available by Helzel-Messtechnic GmbH.
%
%% Input
%  year: year of data collection
%  mo  : month of data collection
%  day : day of data collection
%  hr  : hour of data collection
%
%% Output
%  fpath : path for the wera data 
%  ftime : Pro-name of file with the data for the particular type
%           eg. Name=[ftime,'_hat.SORT'] etc
%
%% Example:      
%  [fpath, ftime]=time2werafile(2010,2,20,10,16)
%  will produce: fpath = /2010/2010051  - the directory path the file is
%                ftime = 20100511016    - the main name of the file
%
%  It can be modified to be used in accordance to file structure in
%  particular applications.
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

% George Voulgaris
% University of South Carolina, USA
% gvoulgaris@geol.sc.edu
% first version: August 27, 2008
%
ddd  = floor(datenum(year,mo,day,hr,mi,0))-datenum(year,1,1,0,0,0)+1;
YYYY = num2str(year);
DDD  = '000'; 
HH   = '00'; 
MM   = '00';
Ddd  = num2str(ddd);
Ndd  = length(Ddd);
DDD(end-Ndd+1:end) = Ddd;
hh   = num2str(floor(hr));
Nhh  = length(hh);
HH(end-Nhh+1:end)=hh;
mm   = num2str(mi);
Nmm  = length(mm);
MM(end-Nmm+1:end)=mm;
ftime = [YYYY,DDD,HH,MM];
fpath = ['/',YYYY,'/',YYYY,DDD];
end