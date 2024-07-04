path(pathdef); % clear previous path changes
addpath(genpath('~/Documents/nano-optics/smarties/'));
addpath(genpath('~/Documents/nano-optics/easyh5/'));
clearvars;

%% example

% requested precision (OA Cext)
accuracy = 1e-10;

% prolate Au spheroid in water
% semi-axes a=b=20nm, c=40nm
% wavelength = 600:50:800; wavelength = wavelength(:);
wavelength = 633:50:800; wavelength = wavelength(:);
Nl = length(wavelength);
epsilon=epsAu(wavelength);
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

    if(abs(stCoa.Cext./stCoa2.Cext-1) > 1.1*accuracy)
        warning('requested precision was not achieved')
    end

    globalN = max(globalN, stParams.N);
    globalnNbTheta = max(globalnNbTheta, stParams.nNbTheta);
end

%% now redo calculations for all wavelengths with these fixed params
globalN = 3;
% allocate 3D array for all results
qmax = 2*(globalN*(globalN + 1) + globalN );
tmatrix = zeros(Nl, qmax, qmax);

for i=1:Nl
    stParams.k1=medium*2*pi/wavelength(i);
    stParams.s=sqrt(epsilon(i)) / medium;

    stParams.N=globalN; stParams.nNbTheta=globalnNbTheta;

    [stCoa, stT] = slvForT(stParams,stOptions);

    [T, q, qp] = exportTmatrix( stT, true, [], [] );

    % convert these indices to (u,up) for treams convention
    [u] = treams_indexing(q, qmax);
    [up] = treams_indexing(qp, qmax);
    % and find linear indices in 3D array to insert these elements
    ind = sub2ind(size(tmatrix),q*0+i, u, up);
    tmatrix(ind) =  T(:,7) + 1i*T(:,8);

end

%% data to export
% vecq = 1:qmax;
% vecs = [repmat(1,1, qmax/2), repmat(2,1, qmax/2)];
% vecp = vecq - (vecs - 1) * qmax/2;
% vecl = floor(sqrt(vecp));
% vecm = vecp - vecl.*(vecl + 1);
% 
% 
% polars = ["electric","magnetic"];
% modes = struct('l', int64(vecl),'m', int64(vecm), 'polarization', polars(vecs));

epsilon = struct(...
    'embedding', medium^2,...
    'particle', epsilon, ...
    'embedding_name', 'H2O, Water', ...
    'embedding_keywords', 'non-dispersive',...
    'embedding_reference', 'constant', ...
    'material_name', 'Au, Gold', ...
    'material_reference', 'Au from Raschke et al 10.1103/PhysRevB.86.235147',...
    'material_keywords', 'dispersive, plasmonic');

geometry = struct('description',  'prolate spheroid', ...
    'shape', 'spheroid','radiusxy', stParams.a, 'radiusz', stParams.c);


computation = struct('description', 'Computation using SMARTIES, a numerically robust EBCM implementation for spheroids', ...
    'accuracy', accuracy, ...
    'Lmax', globalN, ...
    'Ntheta', globalnNbTheta);

comments = struct('name', 'Au prolate spheroid in water',...
    'description', 'Computation using SMARTIES, a numerically robust EBCM implementation for spheroids',...
    'keywords', 'gold, spheroid, ebcm', ...
    'script', [mfilename '.m']);


[f] = tmatrix_hdf5('smarties_spectrum.tmat.h5', tmatrix, wavelength, epsilon, geometry, computation, comments)
