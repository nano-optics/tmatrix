
path(pathdef); % clear previous path changes
addpath(genpath('~/Documents/nano-optics/smarties/'));
addpath(genpath('~/Documents/nano-optics/easyh5/'));
clearvars;


%% actual data

% prolate Au spheroid in water (n=1.33), vacuum wavelength 633nm, 
% epsAu = -11.4+1.181i, semi-axes a=b=20nm, c=40nm, Lmax=3
wavelength = 633;
epsilon=-11.4+1.181i;%epsAu(wavelength);
medium=1.33;
stParams.a=20; stParams.c=40;
stParams.N=3; stParams.nNbTheta=50;

stOptions.bGetR = false;
stOptions.Delta = 0;
stOptions.NB = 0;
stOptions.bGetSymmetricT = false;

stParams.k1=medium*2*pi/wavelength; 
stParams.s=sqrt(epsilon) / medium;
[~, stT] = slvForT(stParams,stOptions);
[T, q, qp] = exportTmatrix( stT, true, [], [] );
qmax = 2*(stParams.N*(stParams.N + 1) + stParams.N );
tmatrix = zeros(qmax, qmax);
ind = sub2ind(size(tmatrix),q,qp);
tmatrix(ind) = T(:,7) + 1i*T(:,8);
l = T(:,3);
m = T(:,5);
s = T(:,1);
% tmatrix2 = sparse(q,qp,T(:,7) + 1i*T(:,8));
% sum(sum(tmatrix==tmatrix2))

%% fake data to export
% wavelength = 633.0;
% tmatrix = [1+1i, 2i; 1, 2];
% modes = struct('l', 1:30,'m', 1:30, 'polarization', repmat(["electric","magnetic"],1,15));
polars = ["electric","magnetic"];
modes = struct('l', l,'m', m, 'polarization', polars(s));

epsilon = struct('embedding', 1.33^2, 'Au', -11.4+1.181i);

geometry = struct('description',  'prolate spheroid', ...
    'shape', 'spheroid','radiusxy', 20.0, 'radiusz', 40.0);

computation = struct('method','EBCM',...
    'software','SMARTIES',...
    'version','1.1',...
    'unit','nm', ...
    'Ntheta', 40, ...
    'accuracy','1e-10');

comments = struct('name', 'Au prolate spheroid in water',...
    'description', 'Computation using SMARTIES, a numerically robust EBCM implementation for spheroids',...
    'keywords', 'gold, spheroid, ebcm', ...
    'script', 'test_dummy.m');


[f, uuid] = tmatrix_hdf5('exporting.h5', tmatrix, modes, wavelength, epsilon, geometry, computation, comments)