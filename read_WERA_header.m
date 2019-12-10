%% read_WERA_header.m
function [WERA]=read_WERA_header(head)
%% [WERA]=read_WERA_header(head)
%  
%  Function used to parse the header from any WERA file. It is called
%  from any function that tries to read any WERA file. All formats have
%  been accounted for but not all of them have been tested. The header is
%  created by the Fortran science package used in the WERA systems developed 
%  by Dr. Klaus-Werner Gurgel at the University of Hamburg and commercially 
%  available by Helzel-Messtechnic GmbH.
%
%% Input  
%  head  = the 512 byte header of the file as a string 
%      
%% Output  
%  WERA.xxx,  a structure array that xxx contains the variables:
%             TL, NFTL, NX, NY, DGT,RAN_OFF, RAN_OFF_m, N_COVERED, SSR, MT,
%             NTL, DOP_RX_OFF, RAN_RX_OFF, HDCODE, LAT, LON, HDG, NORD
%
%  If a file format (stored in WERA.HDCODE) does not contain all variables
%  then the structure value is empty at these locations.
%
%  NOTE: Currently the function can deal with the following code formats:
%
%  FMSOCO, FMSOCN, FMUNSO, FMUNSN, SORT, FMSPC2,  FMSPC3, FMSPEC, FMCALI, 
%  _MSPEC, FMCASW, FMWRAD and FMRADG
%  
%  but it has been tested only for
%     FMSOCO (SORT files) 
%     FMSOCN (RFI files)
%     FMSPC2 (SPEC files)
%     FMSPC3 (SPEC files)
%     FMWRAD (WRAD files)
%     FMRADG (CRAD files)
%
%% Uses
%  parseHeaderSORT.m (Internal function)
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
%% MOodifications
%
%      (1) March 13, 2011. reading the header of SPEC, SPEC2 and SPEC3, the
%      correct station coordinates are now assigned to LAT and LON
%      (2) September 10, 2011. the header for SORT and RFI files is
%      completed.
%      (3) September 26, 2018. Some scaling errors were found in checking
%       if VHF for reading SORT files and in the clock units.
%
% echo off
%------------------------- Parse the header string ------------------------
header = parseHeaderSORT(head);
% ------------------Start Internal Calculations required ------------------
% ---------------------- Parameter Initialization -------------------------
Mon=['JAN';'FEB';'MAR';'APR';'MAY';'JUN';'JUL';'AUG';'SEP';'OCT';'NOV';'DEC'];
c_light    = 299.792458;               % Speed of light
i_clock_90 = 8;
i_CHPdwn   = 128;
WERA.NTL        = [];     WERA.NFTL       = [];    WERA.RATE     = [];
WERA.NX         = [];     WERA.NY         = [];    WERA.NRRANGES = [];
WERA.DGT        = [];     WERA.RAN_OFF    = [];    WERA.PHASECOR = [];
WERA.RAN_OFF_m  = [];     WERA.N_COVERED  = [];
WERA.SSR        = [];     WERA.MT         = [];
WERA.NTL        = [];     WERA.DOP_RX_OFF = [];
WERA.RAN_RX_OFF = [];     WERA.HDCODE     = [];
WERA.LATo       = [];     WERA.LONo       = [];
WERA.LATg       = [];     WERA.LONg       = [];
WERA.HDG        = [];     WERA.NORD       = [];
%--------------------------------------------------------------------------    
% ------------- Information common for all file types  --------------------     
%--------------------------------------------------------------------------
HDCODE      = header.code;                       % Header type
WERA.HDCODE = HDCODE;
WERA.FREQ   = str2double(header.freq);           % Radar Frequency in MHz
qFREQ       = str2double(header.freq);
WERA.fbragg = sqrt(9.81*WERA.FREQ/(c_light*pi)); % Bragg Frequency
WERA.T_chirp= str2double(header.T_chirp);        % Chirp duration in sec
qT_chirp    = WERA.T_chirp;
fs          = 1./WERA.T_chirp;                   % Doppler Spec Sampling Frequency
Last_RC     = str2double(header.nran);           % No of Range Cells
WERA.NSER   = str2double(header.nser);           % Number of samples
WERA.RHF    = str2double(header.rhf);            % Range resolution in km 
WERA.NAnt_SORT = str2double(header.ant);         % No of Antennas in file
WERA.COMMENTS  = header.comment;
%
if ~isempty(header.month) % empty in RAW file
    ao=0; mo=0;
    while ao==0
        mo=mo+1;
        ao=strcmpi(header.month,Mon(mo,:));
    end
    WERA.DATE      = [str2double(header.year)+2000,mo,str2double(header.day)];
    WERA.TIME      = [str2double(header.hh), str2double(header.mm), 0];
    WERA.TIMEBASE  = header.time_base;
    WERA.SITE      = header.site_name;
end
%--------------------------------------------------------------------------
switch upper(HDCODE)

    case {'FMSPEC','FMSPC2','FMSPC3','FMRADG','FMWRAD'}   % All gridded data files
          disp(['Header is ',HDCODE])                     % Specs, CRads and WRAds
          WERA.HDCODE= HDCODE;
          WERA.NTL   = str2double(header.ntl);
          WERA.NFTL  = str2double(header.nftl);
          WERA.NX    = str2double(header.nx);
          WERA.NY    = str2double(header.ny);
          RAN_OFF    = str2double(header.ran_off);        % Range cell offset
          RAN_OFF_m  = RAN_OFF *WERA.RHF*1000.0;
	      I_100      = round(0.02*WERA.RAN_OFF_m)*50;
	      WERA.RAN_OFF    = I_100/(WERA.RHF*1000);
	      WERA.RAN_OFF_m  = RAN_OFF*WERA.RHF*1000.0d0;    % Range offset in meters
          WERA.N_COVERED  = str2double(header.n_covered);
          WERA.LATo  = str2double(header.lat);
          WERA.LONo  = str2double(header.lon);
          WERA.NORD  = str2double(header.nord);
          WERA.LATg  = str2double(header.lat_grid);
          WERA.LONg  = str2double(header.lon_grid);
          WERA.DGT   = str2double(header.dgt);

%     case {'FMRADG','FMWRAD'}   % Gridded Radial and Wave files
%           disp(['Header is ',HDCODE])
%           WERA.HDCODE= HDCODE;
%           WERA.NTL   = str2double(header.ntl);
%           WERA.NFTL  = str2double(header.nftl);
%           WERA.NX    = str2double(header.nx);
%           WERA.NY    = str2double(header.ny);
%           RAN_OFF    = str2double(header.ran_off);      % Range cell offset
%           RAN_OFF_m  = RAN_OFF *WERA.RHF*1000.0;
% 	        I_100      = round(0.02*RAN_OFF_m)*50;
% 	        WERA.RAN_OFF    = I_100/(WERA.RHF*1000);
% 	        WERA.RAN_OFF_m  = RAN_OFF*WERA.RHF*1000.0d0;   % Range offset in meters
%           WERA.N_COVERED  = str2double(header.n_covered);
%           WERA.LATo  = str2double(header.lat);
%           WERA.LONo  = str2double(header.lon);
%           WERA.NORD  = str2double(header.nord);
%           WERA.DGT   = str2double(header.dgt);

    case {'FMSOCO','FMSOCN'}           % SORT & RFI Files case
%           disp(['Header is ',HDCODE])
          WERA.SAMPZ  = str2double(header.nser);
          WERA.HDCODE = HDCODE;
          WERA.SSR    = str2double(header.ss);         % seconds
          WERA.MT     = str2double(header.mt);         % samples
          WERA.NTL    = str2double(header.ntl);
          RAN_OFF     = str2double(header.ran_off);    % Range cell offset
          RAN_OFF_m   = RAN_OFF *WERA.RHF*1000.0;
	      I_100       = round(0.02*RAN_OFF_m)*50;
	      WERA.RAN_OFF   = I_100/(WERA.RHF*1000);
	      WERA.RAN_OFF_m = RAN_OFF*WERA.RHF*1000.0d0;  % Range offset in meters
	      RXOFFSET    = str2double(header.rxoffset);
          WERA.PHASECOR = str2double(header.phase_cor);% Phase correction for each ant
          WERA.LAT    = str2double(header.lat);
          WERA.LON    = str2double(header.lon);
          WERA.NORD   = str2double(header.nord);
          WERA.PWR    = str2double(header.power);
          WERA.RATE   = str2double(header.rate);
          WERA.NRRANGES = str2double(header.nrranges);
          mode        = str2double(header.mode);
          binar       = '000000000000';
          binar1      = dec2bin(hex2dec(header.mode));
          nnn         = length(binar1);
          for i=12:-1:12-nnn+1
              binar(i)=binar1(12-i+1);
          end
       
          if binar(11)=='0'
              WERA.mes(2,:)='Line power for power amplifier is switched OFF';
          else
              WERA.mes(2,:)='Line power for power amplifier is switched ON ';
          end
          if binar(10)=='0'
              WERA.mes(2,:)='Internal calibration mode is switched OFF     ';
          else
              WERA.mes(2,:)='Internal calibration mode is switched ON      ';
          end
          
          if binar(9)=='0'
              WERA.mes(3,:)='IfM internal use (jitter)                     ';
          else
              WERA.mes(3,:)='IfM internal use (jitter)                     ';
          end
          
          if binar(8)=='0'
              WERA.mes(4,:)='90.742153846154 MHz system clock              ';
          else
              WERA.mes(4,:)='90.000000000000 Mhz system clock              ';
          end
          
          if binar(7)=='0'
              WERA.mes(5,:)='Leave RX DDS at DC                            ';
          else
              WERA.mes(5,:)='Set RX DDS to WERA working frequency          ';
          end
          
          if binar(6)=='0'
              WERA.mes(6,:)='IfM internal use (LAH)                        ';
          else
              WERA.mes(6,:)='IfM internal use (LAH)                        ';
          end
          
          if binar(5)=='0'
              WERA.mes(7,:)='Master (chirp up)                             ';
          else
              WERA.mes(7,:)='Slave (chirp down)                            ';
          end
          
          if binar(2)=='0'
              WERA.mes(10,:)='Disable FMICW                                 ';
          else
              WERA.mes(10,:)='Enable FMICW                                  ';
          end
          
          if binar(3)=='0'
              if binar(4)=='0'
                  WERA.mes(9,:)='FMICW sequence #0                             ';
              else
                  WERA.mes(9,:)='FMICW sequence #1                             ';
              end
          else
              if binar(4)=='1'
                  WERA.mes(9,:)='FMICW sequence #2                             ';
              else
                  WERA.mes(9,:)='FMICW sequence #3                             ';
              end
          end
          
          if binar(1)=='0'
              WERA.mes(11,:)='Chirps up or down                             ';
          else
              WERA.mes(11,:)='Triangular Chirps                             ';
          end
          
          %  Check system slave mode
          CL = bitand(uint8(mode),uint8(i_CHPdwn));
          if (CL==i_CHPdwn)
              RXOFFSET = -1 * RXOFFSET;   % System in Slave mode 
          end
          %   Check system clock
          CL = bitand(uint8(mode),uint8(i_clock_90));
          if (CL==i_clock_90)
              Fclk = 10^6*180.0;
	          Divc = 4.0 * 7618.0;	
          else
              Fclk = 10^6*181.48430768;
	          Divc = 4.0 * 7680.0;
          end
          WERA.RXOFFSET = RXOFFSET;
          ADCclk   =  Fclk/Divc;			         % calculate ADC clock
          Tc_exact = str2double(header.mt)/ADCclk;   % calculate exact chirp length
          %
          if (abs((Tc_exact-qT_chirp)) > 0.01*WERA.T_chirp)   % Check if VHF mode
              ADCclk       = 2* ADCclk;
              Tc_exact     = str2double(header.mt)/ADCclk;
              VHF_Sampling = true;
          end
          WERA.T_chirp  = Tc_exact;
          qT_chirp      = Tc_exact;
         % ------------    
              DOP_RX_OFF = RXOFFSET*Fclk/(2^48);
	          RAN_RX_OFF = 0.0;
              if (DOP_RX_OFF > 0.0)
                  while (DOP_RX_OFF > 0.5*Tc_exact)
                      RAN_RX_OFF = RAN_RX_OFF + 1;
                      DOP_RX_OFF = DOP_RX_OFF - 1/Tc_exact;
                  end
                  
              else
                  while (DOP_RX_OFF < -0.5*Tc_exact)
                      RAN_RX_OFF = RAN_RX_OFF - 1;
                      DOP_RX_OFF = DOP_RX_OFF + 1/Tc_exact;
                  end
              end
          WERA.DOP_RX_OFF=DOP_RX_OFF;
          WERA.RAN_RX_OFF=RAN_RX_OFF;
          
      case {'FMCASW','FMCALI'}           % RAW & CAL Files case
%           disp(['Header is ',HDCODE])
          WERA.SAMPZ  = str2double(header.nser);
          WERA.HDCODE = HDCODE;
          WERA.SSR    = str2double(header.ss);        % seconds
          WERA.MT     = str2double(header.mt);        % samples
          WERA.NTL    = str2double(header.ntl);
          RAN_OFF     = str2double(header.ran_off);   % Range cell offset
          RAN_OFF_m   = RAN_OFF *WERA.RHF*1000.0;
	      I_100       = round(0.02*RAN_OFF_m)*50;
	      WERA.RAN_OFF   = I_100/(WERA.RHF*1000);
	      WERA.RAN_OFF_m = RAN_OFF*WERA.RHF*1000.0d0; % Range offset in meters
	      RXOFFSET    = str2double(header.rxoffset);
          WERA.PHASECOR = str2num(header.phase_cor);  % Phase correction for each ant
          WERA.LAT  = str2double(header.lat);
          WERA.LON  = str2double(header.lon);
          WERA.NORD = str2double(header.nord);
          WERA.PWR  = str2double(header.power);
          WERA.RATE = str2double(header.rate);
          WERA.NRRANGES = str2double(header.nrranges);
          mode        = str2double(header.mode);
          binar       = '000000000000';
          binar1      = dec2bin(hex2dec(header.mode));
          nnn         = length(binar1);
          for i=12:-1:12-nnn+1
              binar(i)=binar1(12-i+1);
          end
       
          if binar(11)=='0'
              WERA.mes(2,:)='Line power for power amplifier is switched OFF';
          else
              WERA.mes(2,:)='Line power for power amplifier is switched ON ';
          end
          if binar(10)=='0'
              WERA.mes(2,:)='Internal calibration mode is switched OFF     ';
          else
              WERA.mes(2,:)='Internal calibration mode is switched ON      ';
          end
          
          if binar(9)=='0'
              WERA.mes(3,:)='IfM internal use (jitter)                     ';
          else
              WERA.mes(3,:)='IfM internal use (jitter)                     ';
          end
          
          if binar(8)=='0'
              WERA.mes(4,:)='90.742153846154 MHz system clock              ';
          else
              WERA.mes(4,:)='90.000000000000 Mhz system clock              ';
          end
          
          if binar(7)=='0'
              WERA.mes(5,:)='Leave RX DDS at DC                            ';
          else
              WERA.mes(5,:)='Set RX DDS to WERA working frequency          ';
          end
          
          if binar(6)=='0'
              WERA.mes(6,:)='IfM internal use (LAH)                        ';
          else
              WERA.mes(6,:)='IfM internal use (LAH)                        ';
          end
          
          if binar(5)=='0'
              WERA.mes(7,:)='Master (chirp up)                             ';
          else
              WERA.mes(7,:)='Slave (chirp down)                            ';
          end
          
          if binar(2)=='0'
              WERA.mes(10,:)='Disable FMICW                                 ';
          else
              WERA.mes(10,:)='Enable FMICW                                  ';
          end
          
          if binar(3)=='0'
              if binar(4)=='0'
                  WERA.mes(9,:)='FMICW sequence #0                             ';
              else
                  WERA.mes(9,:)='FMICW sequence #1                             ';
              end
          else
              if binar(4)=='1'
                  WERA.mes(9,:)='FMICW sequence #2                             ';
              else
                  WERA.mes(9,:)='FMICW sequence #3                             ';
              end
          end
          
          if binar(1)=='0'
              WERA.mes(11,:)='Chirps up or down                             ';
          else
              WERA.mes(11,:)='Triangular Chirps                             ';
          end
          
          %  Check system slave mode
          CL = bitand(uint8(mode),uint8(i_CHPdwn));
          if (CL==i_CHPdwn)
              RXOFFSET = -1 * RXOFFSET;   % System in Slave mode 
          end
          %   Check system clock
          CL = bitand(uint8(mode),uint8(i_clock_90));
          if (CL==i_clock_90)
              Fclk = 10^6*180.0;
	          Divc = 4.0 * 7618.0;	
          else
              Fclk = 10^6*181.48430768;
	          Divc = 4.0 * 7680.0;
          end
          WERA.RXOFFSET = RXOFFSET;
          ADCclk   =  Fclk/Divc;			         % calculate ADC clock
          Tc_exact = str2double(header.mt)/ADCclk;      % calculate exact chirp length
          %
          WERA.T_chirp_or  = WERA.T_chirp;
          if (abs((Tc_exact-qT_chirp)) > 0.01*WERA.T_chirp)   % Check if VHF mode
              ADCclk       = 2* ADCclk;
              Tc_exact     = str2double(header.mt)/ADCclk;
              VHF_Sampling = true;
          end
          WERA.T_chirp  = Tc_exact;      
          qT_chirp      = Tc_exact;
         % ------------    
              DOP_RX_OFF = RXOFFSET*Fclk/(2^48);
	          RAN_RX_OFF = 0.0;
              if (DOP_RX_OFF > 0.0)
                  while (DOP_RX_OFF > 0.5*Tc_exact)
                      RAN_RX_OFF = RAN_RX_OFF + 1;
                      DOP_RX_OFF = DOP_RX_OFF - 1/Tc_exact;
                  end 
              else
                  while (DOP_RX_OFF < -0.5*Tc_exact)
                      RAN_RX_OFF = RAN_RX_OFF - 1;
                      DOP_RX_OFF = DOP_RX_OFF + 1/Tc_exact;
                  end
              end
          WERA.DOP_RX_OFF=DOP_RX_OFF;
          WERA.RAN_RX_OFF=RAN_RX_OFF;
     otherwise
          disp(['Header is ',HDCODE,' which is not programmed yet'])
end  
   
end
%
%--------------------------------------------------------------------------
%% parseheader.m
function header=parseHeaderSORT(head0)
%
%  Function that reads the header of a WERA file and parses it into a
%  structure called header that contains on information as a string
%  character. The Format CODE of the file is stored in the header.code
%  variable. It covers the following code formats:
%      FMSOCO, FMSOCN, FMUNSO, FMUNSN,   SORT, FMSPC2
%      FMSPC3, FMSPEC, FMCALI, _MSPEC, FMCASW, FMWRAD and FMRADG 
%  but it has been tested only for
%   FMSOCO (SORT files) 
%   FMSOCN (RFI files)
%   FMSPC2 (SPEC files)
%   FMWRAD (WRAD files)
%   FMRADG (CRAD files)
%   FMCASW (RAW files) no tested  the SCAN files have info further below the
%  comments
%
CDE= ['FMSOCO'; 'FMSOCN'; 'FMUNSO'; 'FMUNSN'; '  SORT'; 'FMSPC2'; 'FMSPC3';
      'FMSPEC'; 'FMCALI'; '_MSPEC'; 'FMCASW'; 'FMRADG'; 'FMWRAD'];
%
N = size(head0,1);
if N~=1
    head0=head0';
end
%
head = char(head0);
ix   = strfind(head,'SAMPLES');
if isempty(ix)==0
    header.nser = head(1:ix-2);
else
    header.nser =-999;
end
%-----------------DECODE HEADER CODE -------------------------------------
header.code      = head(49:54);
if ( strcmp(header.code,CDE(1,:))==1 || strcmp(header.code,CDE(2,:))==1 || strcmp(header.code,CDE(3,:))==1 ...
  || strcmp(header.code,CDE(4,:))==1 || strcmp(header.code,CDE(9,:))==1 || strcmp(header.code,CDE(10,:))==1 ...
  || strcmp(header.code,CDE(11,:))==1)
      header.ant  = head(224:225);
elseif strcmp(header.code,CDE(5,:))==1
      header.ant = '4';      % direction finding sort data;
else
      header.ant = '999';
end
%------------------------------------------------------------------------
if strcmp(header.code, CDE(11,:))==1     % FMCASW
      header.year  = ''; 
      header.day   = '';
      header.month = '';
      header.hh    = '';
      header.mm    = '';
  header.time_base = '';
  header.site_name = '';
      header.freq  = '';
      header.rhf   = '';
      header.nord  = '';
    header.T_chirp = '';
else
      header.year  = head(23:24); 
      header.day   = head(16:17);
      header.month = head(19:21);
      header.hh    = head(26:27);
      header.mm    = head(29:31);
  header.time_base = head(32:34);
  header.site_name = head(37:48);
      header.freq  = head(65:71);
      if strcmp(head(90:94),'RANGE')==1
        header.rhf  = head(97:103);
      else
        header.rhf  = num2str(str2double(head(101:102))*0.15*8);
      end
    header.nord = head(118:120);
 header.T_chirp = head(132:139);
end
   
header.nran        = head(151:153);
header.lon_grad    = head(170:172);
header.lon_min     = head(174:175);
header.lon_sec     = head(177:178);
header.east_or_west= head(180);
header.lat_grad    = head(188:190);
header.lat_min     = head(192:193);
header.lat_sec     = head(195:196);
header.north_or_south = head(198);
lat = str2double(header.lat_grad)+str2double(header.lat_min)/60+str2double(header.lat_sec)/60/60;
if strcmp(header.north_or_south,'S')==1
    header.lat=num2str(-lat);
else
    header.lat=num2str(lat);
end 
lon = str2double(header.lon_grad)+str2double(header.lon_min)/60+str2double(header.lon_sec)/60/60;
if strcmp(header.east_or_west,'W')==1
    header.lon=num2str(-lon);
else
    header.lon=num2str(lon);
end 
%
if strcmp(head(200:202),'MT:')==1
	header.mt = head(203:208);
else
	header.mt = num2str(768 * 2);
end

header.power = head(215:218);

if strcmp(head(239),'T')==0
     header.ran_off = head(239:244);
elseif head(240)~=':'
     header.ran_off = head(240:244);
else
     header.ran_off = head(241:244);
end
 
header.rxoffset    = head(254:263);
 
if strcmp(head(264:265),'SS')==1 || strcmp(head(264:265),'ss')==1
     header.ss = head(267:274);
else
     header.ss = '0.0';
end
 
if (strcmp(header.code, CDE(1,:))==1 || strcmp(header.code, CDE(2,:))==1 || strcmp(header.code, CDE(3,:))==1 ...
  || strcmp(header.code, CDE(4,:))==1 || strcmp(header.code, CDE(9,:))==1 || strcmp(header.code, CDE(10,:))==1 ...
  || strcmp(header.code, CDE(11,:))==1  )
   if strcmp(head(227:229),'MD:')==1
         header.mode = head(230:232);
   else
         header.mode = '';
   end
    %
   if strcmp(head(245:253),'RXOFFSET:')==1
         header.rxoffset = head(254:263);
   else
         header.rxoffset = '0.0';
   end
end
 
if (strcmp(header.code, CDE(6,:))==1 || strcmp(header.code, CDE(7,:))==1  || strcmp(header.code, CDE(8,:))==1 ...
  || strcmp(header.code, CDE(1,:))==1 || strcmp(header.code, CDE(2,:))==1  || strcmp(header.code, CDE(12,:))==1 ...
  || strcmp(header.code, CDE(13,:))==1 )
 
   if strcmp(head(201:203),'NTL')==1
        header.ntl = head(205:208);
   else
        header.ntl = '13';
   end
    
   if strcmp(head(210:213),'NFTL')==1
         header.nftl = head(215:218);
   else
         header.nftl = '512';
   end
end
    
if (strcmp(header.code, CDE(6,:))==1  || strcmp(header.code, CDE(7,:))==1 || strcmp(header.code, CDE(8,:))==1 || ...
    strcmp(header.code, CDE(12,:))==1 || strcmp(header.code, CDE(13,:))==1 ) 
    header.nx   = head(223:225);
    header.ny   = head(230:232);
    header.lat_grid  = head(384+(86:95));
    header.lon_grid  = head(384+(101:111));
    header.dgt       = head(384+(117:127));
end

if (strcmp(header.code, CDE(1,:))==1  || strcmp(header.code, CDE(2,:))==1)
      phase_cor = head(384+(1:128));
      header.phase_cor = reshape(phase_cor,8,16)';
      header.rate      =head(132:139);
      header.nrranges  =head(151:153);
end
if ((strcmp(header.code, CDE(6,:))==1  || strcmp(header.code, CDE(7,:))==1 || strcmp(header.code, CDE(12,:))==1 || ...
     strcmp(header.code, CDE(13,:))==1) && strcmp(head(384+(65:69)),'NCOV ') == 1)
      header.n_covered = head(384+(70:79));
else
      header.n_covered='';    
end
 header.hd          = head(279:281);
 header.comment     = head(283:384);
end