%% Bragg.m
function [k_bragg,f_bragg]  = Bragg(f_radar,h)
%%  Calculate the Bragg Frequency and Wavelength of a radar system.
%
%  function [k_bragg,f_bragg]  = Bragg(f_radar, [h] )
%   
%  Function to calculate the Bragg Wavelength of a radar system.
%  If the water depth is given, the Frequency of the the Ocean Wave
%  that matches the Bragg Wavelength is calculated using the dispersion 
%  relationship. Otherwise deep water conditions are assumed.
%  George Voulgaris, Sept 21, 2012
%
%% Inputs
%          f_radar = Frequency of radar system (in MHz)
%          h       = [optional] water depth (in m)
%                     if no h is given a default value of 1000m is used
%                     that makes f to correspond to deep water solution
%% Outputs
%          k_bragg = Wavenumber of Bragg ocean wave (rads/m)
%          f_bragg = Frequency of Bragg ocean wave (in Hz)
%
%% Example:
% [k_bragg,f_bragg]  = Bragg(8,30)
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
freq    = f_radar(:);                  % radar frequency
c_light = 299.792458;                  % Speed of light (x10^6 m/s)
g       = 9.81;                        % Gravitational acceleration (m/s2)
%
L_bragg = c_light./(2*freq);           % Wavelength of Bragg ocean wave (=0.5*Lem)
k_bragg = 2*pi./L_bragg;               % Wavenumber of Bragg ocean wave
%
if nargin<2
    h=1000;                            % Deep waters
end
h = h(:);
%             
f_bragg = sqrt((g*freq/(c_light*pi)).*k_bragg.*tanh(k_bragg.*h));  % Frequency of Bragg ocean wave at depth h
%
end

