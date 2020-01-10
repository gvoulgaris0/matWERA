%% rads2uv.m 
function [U,VAR,n]=rads2uv(ur,sur,theta)
%% [U,VAR,n]=rads2uv(ur,sur,theta)
%
% Function to combine an arbritary number of radial velocities (n>=2)
% to estimate the 2-D current vector. This is the algorithm used in radar 
% applications when radials from 2 or more radars are available.
%
% The algorithm used is based on:
% 
% Gurgel, K.-W., 1994. Shipborne measurement of surface current fields
% by HF radar (extended version), L'Onde Electrique,74: 54-59.
%   
% It is described in detail in Appendix B of:
%
% Barth, A., A. Alvera-Azcarate, K-W. Gurgel, J. Staneva, A. Port, J-M Beckers
% and E.V. Stanev, 2010. Ensemble perturbation smoother for optimizing
% tidal boundary conditions by assimilation of High-Frequency radar surface
% currents' application to the German Bight. Ocean Sci., 6, 161-178.
%
%% Input 
%
%  ur(1:n)    = vector of n radial velocities (n>=2)
%  sur(1:n)   = corresponding standard deviation (sqrt(variance)) of each radial 
%               velocity (this is an output in WERA systems)
%  theta(1:n) = angle of radial from the horizontal (x-axis, east axis) 
%               measured counterclockwise (mathematical convention)
%
%% Output
%  U          = (u,v) two components of velocity along the x (East) and y (North) axes
%  VAR        = (su2,sv2) corresponding variances of the two velocity components
%  n          = number of radials used for the solution
%
%   
%  Schematic of geometric convention used (assuming 2 radials):
%
%        ^
%  ur2   |th2 /ur1
%    \  .|.  /
%     \. | ./\  
%      \ | /. \th1
% ------\-/------------> x(u)  
%        |
%        | 
%        y(v)
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
% Also note that a python version is available written by Douglas Cahl,
% University of South Carolina.
% 
ur    = ur(:);
sur   = sur(:);
if isempty(find(isnan(ur)))==0
 VAR=[nan;nan];
 U  =[nan;nan];
 n  = nan;
 return
end   
% -------- Check for zero variance ----------------------------------------
iz    = find(sur==0);
if isempty(iz)==0
    sur(iz)=iz*0+1e-3;
end
%--------------------------------------------------------------------------
theta = theta(:);
nu    = length(ur);
ns    = length(sur);
nt    = length(theta);
if (nu~=ns || nu~=nt || nt~=ns)
    disp('All variables must have the same length');
    U   = [nan; nan];
    VAR = U; 
    n   = 0;
    return
end
if (nu<2 || nt<2 || ns<2)
   % disp('You need a minimum of two values to get results');
    U   = [nan; nan];
    VAR = U;
    n   = 1;
    return
end
n   = nu;
A   = [cosd(theta)./sur, sind(theta)./sur];
b   = (ur./sur);
ATA = A'*A;
ATb = A'*b;
C   = ATA.^-1;
U   = ATA\ATb;
VAR = [C(1,1);C(2,2)];
end
