%% read_WERA_asc_cur.m
function [IX,IY,U,V,Uer,Ver,KL]=read_WERA_asc_cur(ftime)
%  [IX,IY,U,V,Uer,Ver, KL]=read_WERA_asc_cur(ftime)
%
%  Function to read the uv WERA data from a cur_asc file and write
%  vectors with locations and velocity components. The cur_asc file has
%  been created by the Fortran science package used in the WERA systems developed 
%  by Dr. Klaus-Werner Gurgel at the University of Hamburg and commercially 
%  available by Helzel-Messtechnic GmbH.
%  
%% Input 
%  ftime: (string) where ftime is the pro_name of the data file
%         this has the format: yyyydddhhmm, where
%  yyyy = year (eg 2008)
%   ddd = yearday of yyyy (day 001, Jan 1 of yyyy)
%    hh = hour of data collection
%    mm = minute of data collection
%
%  You can generate the file name ftime using function time2werafile.m as
%  shown below
%
%  ddd  = floor(datenum(year,mo,day,hr,mi,ss))-datenum(year,1,1,1,0,0,0)+1;
%  YYYY = num2str(year);
%  DDD  = '000'; 
%  HH   = '00';
%  Ddd  = num2str(ddd);
%  Ndd  = length(Ddd);
%  DDD(end-Ndd+1:end) = Ddd;
%  hh   = num2str(floor(hr));
%  Nhh  = length(hh);
%  HH(end-Nhh+1:end) = hh;
%  MM   ='00';
%  ftime=[YYYY,DDD,HH,MM];
%  
%% Output
%  IX  = latitude index
%  IY  = longitude index
%  U   = east current component (m/s)
%  V   = north current component (m/s)
%  Uer = accuracy of east current component (m/s)
%  Ver = accuracy of north current component (m/s)
%
%  Read but not provided as output: K Class information
%  0  : radial information of at least 2 stations are present;
%  1  : radial information only from site 1 is available;
%  2  : radial information only from site 2 is available;
%  ..   ...
%  n  : radial information only from site n is available
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
echo off
% Construct the cur_asc filename fin = [ftime,'_',StID,'.cur_asc'];
fin = ftime;
[fid, message]=fopen(fin,'r');
if isempty(message)==0
    IX = nan;
    IY = nan;
    KL = nan;
    U  = nan;
    V  = nan;
    Uer= nan;
    Ver= nan;
    return
else  
a          = fgetl(fid);
NStations  = str2double(a);
Lat        = zeros(1,NStations);
Lon        = zeros(1,NStations);
for Ist = 1:NStations
    b                = fgetl(fid);
    time(Ist)        = datenum(b(2:18));
    St_name(Ist,1:10)= b(25:34);
    Lon(Ist)         = str2double(b(38:45));
    Lat(Ist)         = str2double(b(54:61));
    ew               = b(63);
    if ew=='W' || ew =='w'
        Lat(Ist)=-Lat(Ist);
    end
end
junk = fgetl(fid);
junk = fgetl(fid);
junk = fgetl(fid);
d    = fgetl(fid);
LAT0 = str2double(d(3:10));
LON0 = str2double(d(13:21));
DGT  = str2double(d(25:30));
NX   = str2double(d(33:35));
NY   = str2double(d(38:40));
%
junk = fgetl(fid);
nos  = fgetl(fid);
N    = str2double(nos);
junk = fgetl(fid);
hunk = fgetl(fid);
  MATR = fscanf(fid,'%3i %3i %8f %8f %3i %8f %8f \n',[7,N]);
  IX = MATR(1,:)';
  IY = MATR(2,:)';
  U  = MATR(3,:)';
  V  = MATR(4,:)';
  KL = MATR(5,:)';
  Uer= MATR(6,:)';
  Ver= MATR(7,:)';
fclose(fid);
end