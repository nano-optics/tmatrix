
path(pathdef); % clear previous path changes
addpath(genpath('~/Documents/nano-optics/smarties/'));
addpath(genpath('~/Documents/nano-optics/easyh5/'));
clearvars;


%% actual data

% prolate Au spheroid in water (n=1.33), vacuum wavelength 633nm, 
% epsAu = -11.4+1.181i, semi-axes a=b=20nm, c=40nm, Lmax=3
wavelength = 600:50:800;
wavelength = wavelength(:);
epsilon=epsAu(wavelength);
medium=1.33;
stParams.a=20; stParams.c=40;
stParams.N=3; stParams.nNbTheta=50;

stOptions.bGetR = false;
stOptions.Delta = 0;
stOptions.NB = 0;
stOptions.bGetSymmetricT = false;
Nl = length(wavelength);
qmax = 2*(stParams.N*(stParams.N + 1) + stParams.N );
tmatrix = zeros(qmax, qmax, Nl);
for i=1:Nl
stParams.k1=medium*2*pi/wavelength(i); 
stParams.s=sqrt(epsilon(i)) / medium;
[~, stT] = slvForT(stParams,stOptions);
[T, q, qp] = exportTmatrix( stT, true, [], [] );
ind = sub2ind(size(tmatrix),q,qp,q*0+i);
tmatrix(ind) =  T(:,7) + 1i*T(:,8);
end

vecq = 1:qmax;
vecs = [repmat(1,1, qmax/2), repmat(2,1, qmax/2)];
vecp = vecq - (vecs - 1) * qmax/2;
vecl = floor(sqrt(vecp));
vecm = vecp - vecl.*(vecl + 1);

%% data to export

polars = ["electric","magnetic"];
modes = struct('l', int16(vecl),'m', int16(vecm), 'polarization', polars(vecs));

epsilon = struct('embedding', medium^2, 'Au', epsilon);

geometry = struct('description',  'prolate spheroid', ...
    'shape', 'spheroid','radiusxy', stParams.a, 'radiusz', stParams.c);

computation = struct('method','EBCM',...
    'software','SMARTIES',...
    'version','1.1',...
    'unit','nm', ...
    'Lmax', stParams.N, ...
    'Ntheta', stParams.nNbTheta, ...
    'accuracy','1e-10');

comments = struct('name', 'Au prolate spheroid in water',...
    'description', 'Computation using SMARTIES, a numerically robust EBCM implementation for spheroids',...
    'keywords', 'gold, spheroid, ebcm', ...
    'script', 'test_export.m');


[f, uuid] = tmatrix_hdf5('exporting_multiple.h5', tmatrix, modes, wavelength, epsilon, geometry, computation, comments)