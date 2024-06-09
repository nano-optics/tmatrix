path(pathdef); % clear previous path changes
addpath(genpath('~/Documents/nano-optics/smarties/'));
addpath(genpath('~/Documents/nano-optics/easyh5/'));
clearvars;

%% example

% requested precision (OA Cext)
accuracy = 1e-10;

% prolate Au spheroid in water
% semi-axes a=b=20nm, c=50nm
wavelength = 350:5:850; wavelength = wavelength(:); 
Nl = length(wavelength);
epsilon=epsilon_Au(wavelength);
medium=1.33;

% constant simulation parameters
stParams.a=20;
stParams.c=40;

% internal options
stOptions.bGetR = false;
stOptions.Delta = 0;
stOptions.NB = 0; % NB will be estimated automatically
stOptions.bGetSymmetricT = false;
stOptions.bOutput = false; % verbosity

%% first, figure out the maximum N's needed

globalN = 1;
globalnNbTheta = 1;

for i=1:Nl
    stParams.k1=medium*2*pi/wavelength(i);
    stParams.s=sqrt(epsilon(i)) / medium;
    % Estimated convergence params
    [N, nNbTheta] = sphEstimateNandNT(stParams, stOptions, accuracy);
    stParams.N=N; stParams.nNbTheta=nNbTheta;
    % Increase params to test accuracy
    stParams2=stParams;
    stParams2.N=stParams2.N+5;
    stParams2.nNbTheta=stParams2.nNbTheta+5;

    [stCoa, stT] = slvForT(stParams,stOptions);
    [stCoa2, stT2] = slvForT(stParams2,stOptions);

    if(stOptions.bOutput)
        fprintf('Convergence testing... lambda = %.5g\n', wavelength(i));
        fprintf('<Cext> = %.10g,   relative error: %.2g\n', stCoa.Cext, abs(stCoa.Cext./stCoa2.Cext-1));
        fprintf('<Csca> = %.10g,   relative error: %.2g\n', stCoa.Csca, abs(stCoa.Csca./stCoa2.Csca-1));
        fprintf('<Cabs> = %.10g,   relative error: %.2g\n', stCoa.Cabs, abs(stCoa.Cabs./stCoa2.Cabs-1));
    end

    if(abs(stCoa.Cext./stCoa2.Cext-1) > 5*accuracy) % some leaway, just order of magnitude
        warning('requested precision was not achieved')
    end

    globalN = max(globalN, stParams.N);
    globalnNbTheta = max(globalnNbTheta, stParams.nNbTheta);
end

%% now redo calculations for all wavelengths with these fixed params

% allocate 3D array for all results
qmax = 2*(globalN*(globalN + 1) + globalN );
tmatrix = zeros(qmax, qmax, Nl);

for i=1:Nl
    stParams.k1=medium*2*pi/wavelength(i);
    stParams.s=sqrt(epsilon(i)) / medium;

    stParams.N=globalN; stParams.nNbTheta=globalnNbTheta;

    [stCoa, stT] = slvForT(stParams,stOptions);

    [T, q, qp] = exportTmatrix( stT, true, [], [] );
    ind = sub2ind(size(tmatrix),q,qp,q*0+i);
    tmatrix(ind) =  T(:,7) + 1i*T(:,8);

end

%% data to export
vecq = 1:qmax;
vecs = [repmat(1,1, qmax/2), repmat(2,1, qmax/2)];
vecp = vecq - (vecs - 1) * qmax/2;
vecl = floor(sqrt(vecp));
vecm = vecp - vecl.*(vecl + 1);

% analytical zeros are those entries SMARTIES did not bother to compute
nonzeros = sub2ind([qmax,qmax],q,qp);
zeros = setdiff((1:qmax^2)', nonzeros);
[zerosq, zerosqp] = ind2sub([qmax,qmax], zeros);
zeros  = struct('q', int16(zerosq), 'qp', int16(zerosqp));

polars = ["electric","magnetic"];
modes = struct('l', int16(vecl),'m', int16(vecm), 'polarization', polars(vecs));

epsilon = struct('embedding', medium^2, 'Au', epsilon);

geometry = struct('description',  'prolate spheroid', ...
    'shape', 'spheroid','radiusxy', stParams.a, 'radiusz', stParams.c);

computation = struct('method','EBCM',...
    'software','SMARTIES',...
    'version','1.1',...
    'unit','nm', ...
    'Lmax', int16(globalN), ...
    'Ntheta', int16(globalnNbTheta), ...
    'analytical_zeros', zeros, ...
    'accuracy', accuracy);

comments = struct('name', 'Au prolate spheroid in water',...
    'description', 'Computation using SMARTIES, a numerically robust EBCM implementation for spheroids',...
    'keywords', 'Au, prolate spheroid, ebcm', ...
    'sources', 'Au from Raschke et al 10.1103/PhysRevB.86.235147',...
    'script', [mfilename '.m']);


[f, uuid] = tmatrix_hdf5('../export/smarties_Au_water_20_50.tmat.h5', tmatrix, modes, wavelength, epsilon, geometry, computation, comments)