path(pathdef); % clear previous path changes
addpath(genpath('~/Documents/nano-optics/smarties/'));
addpath(genpath('~/Documents/nano-optics/easyh5/'));
clearvars;

%% example
% prolate Au spheroid in water
% semi-axes a=b=20nm, c=40nm
wavelength = 600:50:800; wavelength = wavelength(:); 
Nl = length(wavelength);
epsilon=epsAu(wavelength);
medium=1.33;

% constant simulation parameters
stParams.a=20;
stParams.c=40;
stParams.N=3; 
stParams.nNbTheta=90;

% internal options
stOptions.bGetR = false;
stOptions.Delta = 0;
stOptions.NB = 0; % NB will be estimated automatically
stOptions.bGetSymmetricT = false;
stOptions.bOutput = false; % verbosity

%% calculation for all wavelengths

% allocate 3D array for all results
qmax = 2*(stParams.N*(stParams.N + 1) + stParams.N); % size of full T-matrix
tmatrix = zeros(qmax, qmax, Nl);

% indices in TERMS order (2x2 blocks)
vecq = 1:qmax;
vecs = [repmat(1,1, qmax/2), repmat(2,1, qmax/2)];
vecp = vecq - (vecs - 1) * qmax/2;
vecl = floor(sqrt(vecp));
vecm = vecp - vecl.*(vecl + 1);

% loop over wavelengths
for i=1:Nl
    stParams.k1=medium*2*pi/wavelength(i);
    stParams.s=sqrt(epsilon(i)) / medium;

    [stCoa, stT] = slvForT(stParams,stOptions);

    [T, q, qp] = exportTmatrix( stT, true, [], [] );
    % q, qp are the row,col indices in TERMS convention
    % converted to u, up for treams conventions
    [u] = treams_indexing(q, qmax);
    [up] = treams_indexing(qp, qmax);
    ind = sub2ind(size(tmatrix),u,up,u*0+i); 
    % which linear index for elements (u,u',lambdai)
    tmatrix(ind) =  T(:,7) + 1i*T(:,8);

end

% analytical zeros are those entries SMARTIES did not bother to compute
nonzeros = sub2ind([qmax,qmax],q,qp);
zeros = setdiff((1:qmax^2)', nonzeros);
[zerosq, zerosqp] = ind2sub([qmax,qmax], zeros);
zeros  = struct('q', int64(zerosq), 'qp', int64(zerosqp));


%% data to export

polars = ["electric","magnetic"];
modes = struct('l', int64(vecl),'m', int64(vecm), 'polarization', polars(vecs));

epsilon = struct('embedding', medium^2, 'spheroid', epsilon);

geometry = struct('radiusxy', stParams.a, 'radiusz', stParams.c);

method_parameters = struct('Lmax', int64(stParams.N), ...
                           'Ntheta', int64(stParams.nNbTheta));
computation = struct('method_parameters', method_parameters, ...
                     'analytical_zeros', zeros);

comments = struct('name', 'Au prolate spheroid in water',...
    'description', 'Computation using SMARTIES, a numerically robust EBCM implementation for spheroids',...
    'material_reference', 'Au from Raschke et al 10.1103/PhysRevB.86.235147', ...
    'material_spheroid', 'gold', ...,
    'material_embedding', 'water', ...,
    'keywords', 'gold, spheroid, ebcm', ...
    'script', [mfilename '.m']);


[f, uuid] = tmatrix_hdf5('smarties_spectrum.tmat.h5', tmatrix, modes, wavelength, epsilon, geometry, computation, comments)
