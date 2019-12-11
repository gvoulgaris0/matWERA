%% crad
clear
fn = 'C:\mat_wera_v3\4PUBv2\examples\20100401851_hat.crad';
[Time,lat,lon,x,y,u,uvar,uacc,pwr,ang,Range]=read_WERA_crad(fn);
figure, quiver(lon,lat,u,u)
pause
%% spec with SORT file for proper scaling
clear
fn = 'C:\mat_wera_v3\4PUBv2\examples\20100410851_hat.spec';
[Time,LON,LAT,X,Y,freq,fbragg,PXY]=read_WERA_spec(fn);
figure
plot(PXY{30,30})
pause
%% VHF spec without SORT
clear
fn = 'C:\mat_wera_v3\4PUBv2\examples\20172810041_ocr.spec'; 
[Time,LON,LAT,X,Y,freq,fbragg,PXY]=read_WERA_spec(fn);
figure
plot(PXY{100,100})
pause
%% spec without SORT
clear
fn = 'C:\mat_wera_v3\4PUBv2\examples\20161800023_gtn.spec';
[Time,LON,LAT,X,Y,freq,fbragg,PXY]=read_WERA_spec(fn);
pause
figure
plot(PXY{30,30})
%% SORT
clear
fn = 'C:\mat_wera_v3\4PUBv2\examples\20100410851_hat.SORT';
[WERA,t,R,I,Q]=read_WERA_sort(fn);
figure
plot(squeeze(I(4,20,:)))
hold on
plot(squeeze(Q(4,20,:)))
pause
%% RFI
clear
fn = 'C:\mat_wera_v3\4PUBv2\examples\20100410851_hat.RFI';
[WERA,t,R,I,Q]=read_WERA_sort(fn);
figure
plot(squeeze(I(4,20,:)))
hold on
plot(squeeze(Q(4,20,:)))
pause
%% RAW
clear
fn = 'C:\mat_wera_v3\4PUBv2\examples\20193221553_gtn.RAW'
[WERA,t,tc,I,Q,Fs]=read_WERA_raw(fn);
figure
plot(squeeze(I(4,20,:)))
hold on
plot(squeeze(Q(4,20,:)))
pause