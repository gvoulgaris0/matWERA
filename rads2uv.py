import numpy as np
def rads2uv(ur,sur,theta):
# %
# %  Function to combine an arbritary number of radial velocities (n>1)
# %  to estimate 2D current vectors. This is the algorithm used in radar 
# %  applications and WERA in particular.
# %
# %  INPUTS 
# %
# %  ur(1:n)    = n radial velocities
# %  sur(1:n)   = corresponding st. dev (sqrt(variance)) of the radial 
# %               velocities
# %  theta(1:n) = angles of radial from the x-axis measured counterclockwise
# %
# %  OUTPUT
# %  U          = (u,v) two components of velocity along x and y axis
# %  VAR        = (su2,sv2) corresponding variances of velocity components
# %  n          = number of radials used for the solution
# %
# %   Based on:
# % 
# %   Gurgel, K.-W., 1994. Shipborne measurement of surface current fields
# %   by HF radar (extended version), L'Onde Electrique,74: 54-59.
# %   
# %   Described in Appendix B of:
# %
# %   A.Barth,A.Alvera-Azcarate,K-W.Gurgel,J.Staneva,A.Port,J-M Beckers
# %   and E.V. Stanev, 2010. Ensemble perturbation smoother for optimizing
# %   tidal boundary conditions by assimilation of High-Frequency radar surface
# %   currents' application to the German Bight. Ocean Sci., 6, 161-178.
# %   
# %
# %    y(v)^
# %  ur2   |th2 /ur1
# %    \  .|.  /
# %     \. | ./\  
# %      \ | /. \th1
# % ------\-/------------> x(u)  Schematic showing the orientation of variables  
# %        |
# %        | 
# %==========================================================================
# %
# %% Copyright 2019, George Voulgaris, University of South Carolina
# %
# % This file is part of matWERA.
# %
# % matWERA is free software: you can redistribute it and/or modify
# % it under the terms of the GNU General Public License as published by
# % the Free Software Foundation, either version 3 of the License, or
# % (at your option) any later version.
# % 
# % This program is distributed in the hope that it will be useful,
# % but WITHOUT ANY WARRANTY; without even the implied warranty of
# % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# % GNU General Public License for more details.
# % 
# % You should have received a copy of the GNU General Public License
# % along with this program.  If not, see <https://www.gnu.org/licenses/>.
# %
# % If you find an error please notify G. Voulgaris at gvoulgaris@geol.sc.edu
# % 
# % 
# %% Ported to python 1/9/2020 by Douglas Cahl dcahl@geol.sc.edu
	#% -------- Check for zero variance ----------------------------------------
	iz = [i for i, x in enumerate(sur) if x == 0]
	for i in iz:
		sur[i] = 0.001 # minimum variance
	nu    = len(ur)
	ns    = len(sur)
	nt    = len(theta)
	if (nu!=ns or nu!=nt or nt!=ns):
		print('All variables must have the same length')
		U   = [-999, -999]
		VAR = U
		n   = 0
		return U,VAR,n;
	
	if (nu<2 or nt<2 or ns<2):
		print('You need a minimum of two values to get results')
		U   = [-999, -999]
		VAR = U
		n   = 1
		return U,VAR,n;
	
	n   = nu
	ur = np.array(ur)
	sur = np.array(sur)
	theta = np.array(theta)
	theta = theta*np.pi/180
	A   = np.matrix([np.cos(theta)/sur,np.sin(theta)/sur])
	b   = (ur/sur)
	ATA = A*A.T
	ATb = A.dot(b)
	C 	= np.power(ATA,-1)
	U 	= np.linalg.solve(ATA,ATb.T)
	U 	= np.squeeze(np.array(U.T))
	VAR = np.array([C[0,0],C[1,1]])
	return U,VAR,n;

	
# Example, 4 radials
U,VAR,n = rads2uv([1,2,2,1],[1,1.1,1,1],[0,45,90,1]) # Assign returned tuple 
print(U) 
print(VAR)
print(n) 