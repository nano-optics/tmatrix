% possibly multiple wavelengths
wavelength  = (400:50:800)';
Nl = length(wavelength);

Lmax = 3;
qmax = 2*(Lmax*(Lmax+1)+Lmax); % T-matrix size

% dummy 30x30 matrix values for each wavelength
% note the transpose due to HDF5 expecting
% row-major ordering vs matlab's default column-major
tdata = transpose(reshape((1:qmax^2) + 1i*(1:qmax^2), [qmax,qmax]));

tmatrix = zeros(Nl,qmax,qmax);
for i=1:Nl
    tmatrix(i,:,:) = tdata;
end

squeeze(tmatrix(1,1:3,1:3))

% modes, but note that polarization is turned into strings separately
i=1;
for l=1:3
    for m=-l:l
        for s=1:2
            modes.l(i) = int64(l);
            modes.m(i) = int64(m);
            modes.s(i) = int64(s);
            i=i+1;
        end
    end
end
polars = ["electric","magnetic"];
polarization = polars(modes.s);
modes = rmfield(modes,'s');

% dummy 'analytical zeros' for e.g. EBCM methods
% [zerosq, zerosqp] = ndgrid(1:2:30, 1:2:30);
% zeros  = struct('q', zerosq, 'qp', zerosqp);

% materials
embedding = struct('relative_permeability', 1.0, ...
                   'relative_permittivity', 1.33^2);
particle = struct('relative_permeability', 1.0, ...
                  'relative_permittivity', repmat(-11.4+1.181i, [Nl,1]));

% geometry
geometry = struct('radiusxy', 20.0, 'radiusz', 40.0);

scatterer = struct('material', particle, ...
                  'geometry', geometry);

% details about computation

method_parameters = struct('Lmax', int64(3), ...
                           'Ntheta', int64(100));

computation = struct('method_parameters', method_parameters);
% 'analytical_zeros', zeros can be added here

script = convertCharsToStrings(fileread('export_matlab.m'));

% combined (almost all) information into one struct
s = struct('tmatrix', tmatrix, ...
    'vacuum_wavelength', wavelength, ...
    'embedding', embedding,...
    'scatterer', scatterer, ...
    'modes', modes, ...
    'computation', computation);

addpath(genpath('../easyh5/'));

f = 'am.tmat.h5';
saveh5(s, f, 'ComplexFormat', {'r','i'}, 'rootname', '', 'Compression', 'deflate'); 

% deal with string objects manually
h5create(f,'/computation/files/script', size(script), 'Datatype', 'string')
h5write(f,'/computation/files/script', script)

h5create(f,'/modes/polarization', size(polarization), 'Datatype', 'string')
h5write(f,'/modes/polarization', polarization)

% root attributes
h5writeatt(f, '/', 'name', 'Au prolate spheroid in water');
h5writeatt(f, '/', 'storage_format_version', 'v0.01'); 
h5writeatt(f, '/','description', ...
    'Computation using SMARTIES, a numerically robust EBCM implementation for spheroids');
h5writeatt(f, '/','keywords', 'gold, spheroid, ebcm, passive, reciprocal, czinfinity, mirrorxyz');

% object and group attributes
h5writeatt(f, '/vacuum_wavelength', 'unit', 'nm');

h5writeatt(f, '/embedding', 'keywords', 'non-dispersive');
h5writeatt(f, '/embedding', 'name', 'H2O, Water');

h5writeatt(f, '/scatterer/material', 'name', 'Au, Gold');
h5writeatt(f, '/scatterer/material', 'reference', 'Au from Raschke et al 10.1103/PhysRevB.86.235147');
h5writeatt(f, '/scatterer/material', 'keywords', 'dispersive, plasmonic');

h5writeatt(f, '/scatterer/geometry', 'name', 'homogeneous spheroid with symmetry axis z');
h5writeatt(f, '/scatterer/geometry', 'unit', 'nm');
h5writeatt(f, '/scatterer/geometry', 'shape', 'spheroid')

h5writeatt(f, '/computation', 'description', 'Computation using SMARTIES, a numerically robust EBCM implementation for spheroids');
h5writeatt(f, '/computation', 'name', 'SMARTIES');
h5writeatt(f, '/computation', 'method', 'EBCM, Extended Boundary Condition Method');
h5writeatt(f, '/computation', 'software', 'SMARTIES');
h5writeatt(f, '/computation', 'version', '1.1');
