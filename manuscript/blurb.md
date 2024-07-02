
== SMARTIES

The original T-matrix method, devised by Waterman [ref], introduced alongside a specific calculation scheme -- the Extended Boundary Condition Method. This technique has strong analytical roots, requiring no meshing of the particle's volume or surface and, instead, computes T-matrix elements via analytical formulas which reduce to Mie theory for spherical particles [ref]. For other shapes, the computation requires integration over the particle surface, by numerical quadrature. For axi-symmetric particles, the method is remarkably efficient as the matrix elements are obtained via simple one-dimensional integrals. The EBCM method is particularly popular for simple geometrical shapes, where it typically provides the fastest and most accurate way to calculate a T-matrix [ref]. 

SMARTIES is a Matlab implementation of the EBCM to simulate the optical properties of oblate and prolate spheroidal particles, with comparable speed, convenience and accuracy as Mie theory for spheres. SMARTIES is only applicable to spheroidal particles, for which it uses an improved algorithm that overcomes some of the numerical difficulties related to loss of precision faced by EBCM in the case of large and elongated particles [ref]. The code may be useful to researchers seeking a fast, accurate and reliable tool to simulate the near-field and far-field optical properties of elongated particles, but can also appeal to other developers of light-scattering software seeking a reliable benchmark for non-spherical particles with a challenging aspect ratio and/or refractive index contrast.

We provide below an example script to output the T-matrix of a gold spheroid, used in the calculation of Figure XX.

```
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

% analytical zeros are those entries SMARTIES did not need to compute
nonzeros = sub2ind([qmax,qmax],q,qp);
zeros = setdiff((1:qmax^2)', nonzeros);
[zerosq, zerosqp] = ind2sub([qmax,qmax], zeros);
zeros  = struct('q', int64(zerosq), 'qp', int64(zerosqp));

%% data to export
polars = ["electric","magnetic"];
modes = struct('l', int64(vecl),'m', int64(vecm), 'polarization', polars(vecs));

epsilon = struct('embedding', medium^2, 'spheroid', epsilon);
geometry = struct('spheroid', struct('radiusxy', stParams.a, 'radiusz', stParams.c));
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
```


== TERMS

TERMS is a Fortran program based on the superposition T-matrix method, designed to simulate the near-field and far-field optical properties of collections of particles [ref]. It was developed primarily to model relatively compact clusters of resonant scatterers, such as plasmonic particles often requiring large multipolar orders [ref]. TERMS implements several independent algorithms, with complementary strengths and weaknesses, to describe the self-consistent electromagnetic interaction between multiple scatterers, and from there compute far-field optical properties such as absorption, scattering, extinction, circular dichroism, as well as near-field intensities and the local degree of optical chirality. By describing the incident and scattered fields in a basis of spherical waves the T-matrix framework lends itself to analytical formulas for orientation-averaged quantities such as far-field cross-sections and near-field quantities, greatly reducing the computational time needed to simulate particles and systems of particles in random orientation [ref].

Each scatterer is described by a T-matrix, which is computed internally for spherical particles (including layered spheres), or using external files computed with any other method.
 
TERMS computations are divided into three main modes:

- Far-field quantities (absorption, scattering, extinction, circular dichroism) for multiple wavelengths and angles of incidence, as well as orientation-averages
- Near-field calculations for multiple wavelengths and incident angles, also computing the local degree of chirality, as well as orientation-averages
- Stokes parameters and differential scattering cross-sections for multiple incidence or scattering angles


The program's documentation and website offer many examples of use [ref]; for the purpose of this work we only illustrate the import of an external T-matrix in the `tmat.h5` format. The input file for the simulation reproduced below considers two gold spheroids in water, separated by 100 nm and rotated by 45 degrees to form a chiral structure. 

```
ModeAndScheme 2 3
MultipoleCutoff 5
Wavelength 400 800 200
Medium 1.7689 # water

# dimer of Au spheroids
Scatterers 2
TF1 0 -50 0.0 50 0.0 0.0 0.0  2.5
TF1 0  50 0.0 50 0.0 0.7853982 0  2.5
```

The simulation is run with the command

````
terms input > log
````

and outputs cross-sections in the file `results.h5`, displayed in Figure XX. For comparison, the same simulation was run with a T-matrix produced by SMARTIES (Script XX) for the same geometry. 


= Tools for conversion

A number of open-source programs are available to compute T-matrices, but many of them do not (yet) implement the output format presented herein. While we encourage the community to add this functionality in order to fully benefit from interoperability between programs, it can also be useful, as a short-term or one-time workaround, to _convert_ T-matrix data stored in a different form. One example is the "long format" used to store T-matrix entries in earlier versions of SMARTIES [ref], or Scuff-EM [ref], or PXTAL [ref]. We include example scripts at [url] to reshape such data and produce a standard h5 format. 

```
d = read.table('data/tmat_Au20x40_Nmax3.tmat')
names(d) = c('s','sp','l','lp','m','mp','Tr','Ti')
head(d)
  s sp l lp  m mp            Tr            Ti
1 1  1 1  1 -1 -1 -6.049214e-05 -4.266526e-04
2 1  1 1  1  0  0 -3.331557e-05 -3.932179e-04
3 1  1 1  1  1  1 -6.049214e-05 -4.266526e-04
4 1  1 1  3 -1 -1 -2.374705e-07 -1.995117e-06
5 1  1 1  3  0  0 -1.110299e-07 -1.278537e-06
6 1  1 1  3  1  1 -2.374705e-07 -1.995117e-06
```

Conversion to wide format, and export as `.tmat.h5`, can then be done by adding the required geometry and material information to make a complete entry. Basic export scripts are available in 4 different languages (R, Julia, Matlab, Python) to serve as examples for similar conversion tasks.

