%% read_WERA_crad.m
function [Time,lat,lon,x,y,u,uvar,uacc,pwr,ang,Range]=read_WERA_crad(fname,GEO,Ix,Jy)
% [Time,lat,lon,x,y,u,u_var,u_acc,pwr,ang,Range]=read_WERA_crad(fname,[GEO],[Ix],[Jy])
% 
% Function that reads the information from a WERA crad file created by the 
% WERA systems developed by Dr. Klaus-Werner Gurgel at the University of 
% Hamburg and commercially available by Helzel-Messtechnic GmbH.
%
%  [Time,lat,lon,x,y,u,u_var,u_acc,pwr,ang,Range]=read_WERA_crad(fname,~,(10:130),(20:120))
%  [Time,lat,lon,x,y,u,u_var,u_acc,pwr,ang,Range]=read_WERA_crad(fname)
%  [Time,lat,lon,x,y,u,u_var,u_acc,pwr,ang,Range]=read_WERA_crad(fname,'UTM')
%
%% Input
%  fname = filename of the crad file (including extention)
%  GEO   = [Optional] if ='UTM', then the grid is also converted to 
%          UTM coordinate system (x,y). Deafault is Geographical coord.
%  Ix    = [Optional] Array of indices to be loaded (i.e., Ix=10:130)
%  Jy    = [Optional] Array of indices to be loaded (i.e., Jy=20:40)
%          If Ix, Jy are ommitted ALL values are loaded 
%
%% Output
%  Time     = Time in MATLAB datenum format at the middle of burst
%  lat, lon = Latitude and longitude in degrees (negative for W and S)
%  x, y     = Orthogonal coordinate system with 0,0 being at the
%             WERA station location (in meters) [ if GEO=='UTM'].
%  u        = radial velocity (m/s)
%  u_var    = standard deviation of radial velocity (m/s)
%  u_acc    = accuracy of radial velocity estimate  (m/s)
%  pwr      = power of signal
%  ang      = angle of radial (calculated with respect to the othog. grid, 
%             0 degs, radial along the east (x) axis going eastward
%             90 degs, radial along the north (y) axis going northward). 
%             To convert to beam angle use: (HDG-pi+ang).
%  Range    = Range of grid point in Km
%
%% Notes
%  (1) To plot the results use:   quiver(x,y,u.*cosd(ang),u.*sind(ang))
%  (2) If you want to select part of the grid only utilize the parameters Ix, Jy
%
%% Uses
%  read_WERA_header.m
%  geog2utm.m  [only if GEO=='UTM']
%  WGS84v.m 
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
%
%       01/22/2011 - Time at the center of the burst is estimated & output.
%
if nargin<4 
    Jy=[]; Ix=[];
end
if nargin<2 || isempty(GEO)
    GEO = 'WGS';
end
% -------------------- Open/ Start Reading the data file  -----------------
if strcmpi(fname(end-2:end),'RAD')==1  
    filename  = fname;
    s         = dir(filename);
    if isempty(s)==1
        error([fname, ' is not found, check path and/or syntax, Aborting'])
    end 
    filebytes = s.bytes;
else
    error([fname, ' is not a RAD type file file, Aborting'])
end
%
if filebytes>512 
    fid = fopen(filename,'rb','ieee-le.l64');  % Open the file as a binary
    fseek(fid, 4, 'cof'); 
    head = fread(fid,512,'char');
    fseek(fid, 4, 'cof');
else
    error([filename, ' size is less than 512 bytes, Aborting'])
end
%--------------------------------------------------------------------------
%                          Read the Header Information
%--------------------------------------------------------------------------
[WERA] = read_WERA_header(head);
NX     = WERA.NX;
NY     = WERA.NY;
DGT    = WERA.DGT;
LONo   = WERA.LONo;
LATo   = WERA.LATo;
GRAD   = 180.0/pi;
NCOV   = WERA.N_COVERED;
% --- Estimate time stamp at the middle of the data collection period ---
dt           = 0.5*WERA.NSER*WERA.T_chirp;    % Middle of data in secs
WERA.TIME(3) = WERA.TIME(3)+dt;               % Add the extra secs
Time         = datenum([WERA.DATE WERA.TIME]);
%
% --- Read all the matrix or selected part as defined in function inputs
if isempty(Jy) || nargin<4
    SELy   = 1:NY;
else
    SELy   = Jy;
end
if isempty(Ix) || nargin<3
    SELx   = 1:NX;
else
    SELx   = Ix;
end
%--------------------------------------------------------------------------
%                         Read Data from the File
%--------------------------------------------------------------------------
switch upper(WERA.HDCODE)
    
    case {'FMRADG'}
        fseek(fid, 4, 'cof'); 
        NTL  = fread(fid,1,'int16'); % No of fractional series
        KSIG = fread(fid,1,'int16');
        KWIN = fread(fid,1,'int16');
        NX   = fread(fid,1,'int16');
        NY   = fread(fid,1,'int16');
        I2NN = fread(fid,1,'int16'); 
     fseek(fid, 4, 'cof'); 
        NC=NX*NY; 
     fseek(fid, 4, 'cof'); 
        LAT   = fread(fid,NC,'float32')*180/pi;
     fseek(fid, 4, 'cof');
     fseek(fid, 4, 'cof'); 
        LON   = fread(fid,NC,'float32')*180/pi;
     fseek(fid, 4, 'cof'); 
     fseek(fid, 4, 'cof'); 
        KUR   = fread(fid,NC,'int32')  ;   % Number of estimates in a grid
     fseek(fid, 4, 'cof'); 
     fseek(fid, 4, 'cof'); 
        URM   = fread(fid,NC,'float32');   % Sum i=1 to Kur of radial_vel(i) x SNR(i)
     fseek(fid, 4, 'cof'); 
     fseek(fid, 4, 'cof'); 
        URS   = fread(fid,NC,'float32');   % Sum i=1 to Kur of radial_vel(i)^2 x SNR(i)
     fseek(fid, 4, 'cof'); 
     fseek(fid, 4, 'cof'); 
        USM   = fread(fid,NC,'float32');   % Sum i=1 to Kur of SNR(i)
     fseek(fid, 4, 'cof'); 
     fseek(fid, 4, 'cof'); 
        PWR   = fread(fid,NC,'float32');   % Mean of Power
     fseek(fid, 4, 'cof'); 
     fclose(fid);
     
    case {'FUTURE'}
        % To add future cases/formats
        fclose(fid);
        
    otherwise
    disp('Header FMRADG was expected but not found');
    fclose(fid);
    return
    
end
%--------------------------------------------------------------------------
%                   Convert Data to Radial Velocities
%--------------------------------------------------------------------------
% ---------------- Set up the variables -----------------------------------
UR   = ones(size(LAT))*NaN;   
SIG  = UR;
SIGm = UR;
%--------------------------------------------------------------------------
Ik         = find(KUR>0);  % 
UR(Ik)     = URM(Ik)./USM(Ik);                 % Radial Velocity (m/s)
SIG(Ik)    = sqrt(URS(Ik)./USM(Ik)-UR(Ik).^2); % Std. Deviation (m/s)
In         = find(SIG(Ik)<0);
SIG(Ik(In))= SIG(Ik(In))*0;      
SIGm       = SIG./sqrt(KUR);                   % Accuracy (m/s)
%----- Use the WGS84 Ellipsoid to estimate ranges and angles --------------
[range,Ang] = WGS84v(LAT,LON,LATo,LONo);
Ang         = 90 - Ang;     % Convert to angle from East axis positive anticlockwise
%--------------------------------------------------------------------------
%---- Geographical to orthogonal UTM coordinate transformation ------------
if strcmp(GEO,'UTM')==1
    [X,Y] = geog2utm(LON,LAT,LONo,LATo);  % Geographic 2 Orthogonal
    X     = X*1000;                       % Convert from km to m
    Y     = Y*1000;                       % Convert from km to m
    x     = reshape(X,NX,NY);
    x     = x(SELx,SELy);
    y     = reshape(Y,NX,NY);
    y     = y(SELx,SELy);
else
    x     = 0;
    y     = 0;
end
%--------------------------------------------------------------------------
ig   = find(isnan(UR)==0);
dm   = median(UR(ig));
stdm = std(UR(ig));
%--------------------------------------------------------------------------
u    = reshape(UR,NX,NY);       u = u(SELx,SELy);
lat  = reshape(LAT,NX,NY);    lat = lat(SELx,SELy);
lon  = reshape(LON,NX,NY);    lon = lon(SELx,SELy);
pwr  = reshape(PWR,NX,NY);    pwr = pwr(SELx,SELy);
uvar = reshape(SIG,NX,NY);   uvar = uvar(SELx,SELy);
uacc = reshape(SIGm,NX,NY);  uacc = uacc(SELx,SELy);
ang  = reshape(Ang,NX,NY);   ang  = ang(SELx,SELy);      
Range= reshape(range,NX,NY);Range = Range(SELx,SELy);
%
end
