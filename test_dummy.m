addpath(genpath('../easyh5/'));
clearvars;


% possibly multiple wavelengths
wavelength  = (400:50:800)';
Nl = length(wavelength);

% dummy 30x30 matrix values for each wavelength
tdata = reshape((1:900) + 1i*(1:900), [30,30]);
tmatrix = repmat(tdata, [1,1,Nl]);

% modes, but note that polarization is handled separately
modes = struct('l', int16(1:30), 'm', int16(1:30)); 
polarization = repmat(["electric","magnetic"], 1, 15);

% dummy 'analytical zeros' for e.g. EBCM methods
[zerosq, zerosqp] = ndgrid(1:2:30, 1:2:30);
zeros  = struct('q', zerosq, 'qp', zerosqp);

% materials
embedding = struct('relative_permeability', 1.0, ...
                   'relative_permittivity', 1.33^2);
particle = struct('relative_permeability', 1.0, ...
                  'relative_permittivity', repmat(-11.4+1.181i, [1,Nl]));
materialname = 'Au'; 
materials = struct('embedding', embedding, materialname, particle);

% geometry
geometry = struct('shape', 'spheroid','radiusxy', 20.0, 'radiusz', 40.0);

% details about computation, including full script
computation = struct('method','EBCM',...
    'software','SMARTIES',...
    'version','1.1',...
    'Ntheta', 40, ...
    'accuracy','1e-10', ...
    'analytical_zeros', zeros, ...
    'script', fileread('test_dummy.m'));

% combined (almost all) information into one struct
a = struct('tmatrix', tmatrix, ...
    'vacuum_wavelength', wavelength, ...
    'embedding', embedding,...
    'materials', materials, ...
    'geometry', geometry, ...
    'modes', modes, ...
    'computation', computation, ...
    'uuid', char(matlab.lang.internal.uuid()));    


%% save to file
f = 'am.tmat.h5';
saveh5(a,f, 'ComplexFormat', {'r','i'}, 'rootname', '', 'Compression', 'deflate'); 

h5create(f,'/modes/polarization', size(polarization), 'Datatype', 'string')
h5write(f,'/modes/polarization', polarization)
h5writeatt(f, '/', 'name', 'Au prolate spheroid in water');
h5writeatt(f, '/', 'created_with', 'Matlab easyh5');
[major,minor,rel] = H5.get_libversion();
h5writeatt(f, '/', 'storage_format_version', sprintf('%d.%d.%d',major,minor,rel));

h5writeatt(f, '/','description', ...
    'Computation using SMARTIES, a numerically robust EBCM implementation for spheroids');
h5writeatt(f, '/','keywords', 'gold, spheroid, ebcm');
h5writeatt(f, '/vacuum_wavelength', 'unit', 'nm');
h5writeatt(f, '/uuid', 'version', '4');
h5writeatt(f, '/geometry', 'name', 'prolate spheroid');

