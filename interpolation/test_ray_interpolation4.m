%%
addpath('D:\\projects\\WavesCNR\\matlab\\dynamic\\dual_quaternions');
close all;
clear;

%%
Rp = [1    1  0;
      0   0  0];
  
Rv = [0   0  1;
      0.4472   0  0.8944];
  
nrays = size( Rp, 1);
for ii=1:nrays
    Rv(ii,:)=Rv(ii,:)/norm(Rv(ii,:));
end 

[dq, midpoint ] = rays_min_transform(Rp,Rv);
dq
midpoint