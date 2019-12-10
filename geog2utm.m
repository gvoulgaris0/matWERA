%% geog2utm.m
function [X,Y]=geog2utm(lon,lat,LON,LAT)
% [X,Y]=geog2utm(lon,lat,[LON],[LAT])
%
%  Function to convert geographic coordinates (Longitude & %  latitude) to 
%  coordinates on the Universal Transverse Mercator (Elipsoid) system. 
%  The program finds automatically the UTM zone from the center of latitudes.
%
%  The algorithm is based on: J.P. Snyder, 1987, Map Projections - A Working 
%                            Manual. USGS Professional Paper 1395.
%
%% Input	
%  lon: longitude in degrees and decimals of a degree 
%  lat: latitude in degrees and decimals of a degree
%  LON: [OPTIONAL] Longitude of the origin point for X and Y (0.0,0.0)
%  LAT: [OPTIONAL] Latitude of the origin point for X and Y 
%
%  Note: South & West must be negative values (i.e. 5.2N 3.3W is +5.2 -3.3).
%                          
%% Output:
%  X: Eastings (in km) from the origin point 
%  Y: Northings (in km) from the origin point 
%   
%  The origin point is that of the UTM zone (if no LON, LAT are given) or
%  (0,0) is the location of the Tx array (if the optional values LON, LAT
%  are given).
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
%  Revision 1, Jan. 2011: Took care of the UTM 500km offset automatically.
%
echo off
%
lat     = lat(:);
lon     = lon(:);
UTMzone = 0;
%
if (nargin > 4 || nargin==1 || nargin==3)
   error('geog2utm.m: Only 2 or 4 arguments are allowed')
end
% Find the UTM Zone and corresponding Meridian using the mean longitude
if nargin == 2 
    LL          = min(lon)+(max(lon)-min(lon))/2;
    UTMrange    = (-180:6:180)';
    UTMmeridian = (-177:6:178)';
    for i = 1:60
      if (UTMrange(i)<LL && UTMrange(i+1)>=LL) 
          UTMzone=i; 
      end
    end
    XOFF = 500.0;   % The 500Km offset required for UTM coordinates
    LAT  = 0;
    LON  = UTMmeridian(UTMzone);
else
    XOFF = 0.0;     % In case of local coordinate system with given origin
end
%
% Constants
%
disp([' > UTM Zone:',num2str(UTMzone),' ORIGIN: ',num2str(LON),'E, ',num2str(LAT),'N <'])
%
alpha   = 6378206.4;
e2      = 0.00676866;
DEG2RAD = (2*pi/360);
ko      = 0.9996;
%
% Converting to rads
%
lat = lat*DEG2RAD;
lon = lon*DEG2RAD;
LAT = LAT*DEG2RAD;
LON = LON*DEG2RAD;
%
ee2 = e2/(1-e2);
N   = alpha./sqrt((1-e2*sin(lat).^2));
T   = tan(lat).^2;
C   = ee2*cos(lat).^2;
A   = (lon-LON).*cos(lat);
M   = 111132.0894*lat/DEG2RAD-16216.94*sin(2*lat)+17.21*sin(4*lat)-0.02*sin(6*lat);
Mo  = 111132.0894*LAT/DEG2RAD-16216.94*sin(2*LAT)+17.21*sin(4*LAT)-0.02*sin(6*LAT);
%
X   = ko*N.*(A+(1-T+C).*(A.^3/6)+(5-18*T+T.^2+72*C-58*ee2).*(A.^5/120));
Y   = ko*(M-Mo+N.*tan(lat).*((A.^2/2)+(5-T+9*C+4*C.^2).*(A.^4/24)+(61-58*T+T.^2+600*C-330*ee2).*(A.^6/720)));
%
X   = X/1000+XOFF; % in Km
Y   = Y/1000;      % in Km
end