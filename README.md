A Benchmark for Surface Reconstruction
--------------------------------------
This repository contains a working version of the code for the paper [A Benchmark for Surface Reconstruction](http://vgc.poly.edu/files/berger/recon_bench/paper/tog.pdf) by Berger et. al.
This repository modifies the original source code of the paper to compile on Linux and MacOS without changing the functionality. The repository includes the source code for surface modeling, sampling, reconstruction, evaluation, and plotting results. Here we provide details regarding the various executables necessary for each of these tasks. Throughout the description, we have provided example data and instructions on how to process the data to produce error plots.

## Compilation
This version of the reconstruction benchmark uses CMake to generate build files. It's beeen tested on OSX and Ubuntu 18.04.
The benchmark requires that the following dependencies be installed on the system for compilation.
- libpng12
- lapack/blas
- ffmpeg
- gnuplot
- epstopdf
- cmake
- OpenEXR
- libtiff

### Installing Dependencies on Ubuntu
To install the dependencies on Ubuntu, type
```
sudo apt-get install libpng12-dev liblapack-dev libblas-dev ffmpeg gnuplot texlive-font-utils libtiff-dev openexr cmake
```

### Building the project on Mac and Linux
To build the project on Mac and Linux you can type the following from the root directory of the repository:
```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```


## Surface Modeling
We allow for the creation of polygonal MPU implicit surfaces, obtained through approximating
triangulated surface meshes:
```
./build/mesh_to_implicit surface_mesh implicit_surface min_samples fit_epsilon covering
```
where:
- `surface_mesh` can be an off, obj, or ply mesh.
- `implicit_surface` is the name of the output MPU surface -> needs to have .mpu as an extension, for later use.
- `min_samples` specifies the minimum number of triangles necessary to fit a shape function. It also has an effect on CSG operations. A value of 6 tends to be a good trade-off.
- `fit_epsilon` specifies the quality of a fit, as a percentage of the bounding box of the input mesh -> if satisfied, subdivision is terminated in the octree. This is largely dependent on the complexity of the shape, and the tesselation. Typical values range from 0.005 - 0.01. Larger values result in details being smoothed out.
- `covering` specifies the radius of the sphere which occupies each octree cell, specified as a fraction of half the bounding box diagonal. Typical values range from 1.0 - 1.25.

We note that creating MPU polygonal surfaces from triangle meshes can be a trial-and-error process, as
certain surface meshes may be difficult to fit shape functions to. Hence to facilitate this, we have
provided a marching cubes implementation so that one may quickly observe the resulting zero-set of
the implicit function:
```
./build/isosurface implicit_surface resolution output_surface
```
where:
- `implicit_surface` is the aforementioned MPU surface.
- `resolution` specifies the marching cubes grid resolution. This should be rather high (i.e. 256), in order to stay within the bounds of the MPU surface definition.
- `output_surface` is the resulting isosurface, either off or obj format.

Included in `modeling/models` is an example implicit surface, bumps, used in the "simple shapes" part
of the benchmark.


## Synthetic Scanning
From the MPU surfaces we next allow for synthetically scanning the surfaces, simulating the process of
an optical triangulation-based scanner. We break up point cloud generation into generation of configuration
files, followed by executing the coniguration files to obtain the point cloud. To generate configuration files:

```
./build/pc_generator implicit_surface config_base ([param value])* ([param range min_value max_value number])*
```
where:
- `implicit_surface` is an MPU surface file - may be specified as absolute path, or if in reconbench directory,
can specify as relative.
- `config_base` represents the base file name at which the configuration, and subsequently point cloud, files
will be stored. Depending on the number of parameters set, files will be labeled config_base_0.cnf, config_base_1.cnf,
etc.. May be specified as absolute path, or if in reconbench directory, can specify as relative
- `param value`: assign a scanner parameter to value. Please see sampler/pc_generator.cpp, as well as the original paper,
for further description on the parameters.
- `param range min_value max_value number`: for a given parameter, assign a range of values, starting from min_value
to max_value, in uniform increments, with number being the amount generated.

It is necessary to at least supply the image resolution and number of scans. All other parameters are optional, with
defaults already at hand - see sampler/UniformSampler.cpp for default parameters.

To generate the point cloud, from within the reconbench directory:

```
python scripts/RunSampler.py config_file
```
where `config_file` is the aforementioned configuration file generated through pc_generator. The result is a `.npts` file, as well as a `.mov` file, which is a movie of all scans and laser stripes taken through the scanning simulation.

We have provided some example config files for the bumps shape, varying in increasing resolution, found in `data/pcs/bumps`. Give the above a try to produce the point clouds for these config files.


## Reconstruction
We have provided a script to more easily run reconstruction algorithms on the generated point clouds. As
every reconstruction algorithm has its own set of parameters, we require the user to provide a script to run
their algorithm, and to modify `scripts/scripts_recon.py` to support their algorithm. As an example, we have included
Poisson Surface Reconstruction and its associated script, in the recon directory. Paths can be either absolute,
or relative to the reconbench directory. See `data/meshes/bumps/recon_config.cnf` to see how to set parameters
to your algorithm. Once all set, algorithms may be run in batch by:

```
python scripts/scripts_recon.py config_file
```

We suggest compiling Poisson Surface Reconstruction, and running the above command on the "bumps" point clouds produced above to get a feel for the reconstruction script and configuration.


## Evaluation
Evaluation requires: the MPU implicit surface, a dense uniform sampling of the surface and the output
reconstructed mesh. We have provided dense uniform samplings used in our benchmark in the modeling/models
directory. However, if you have generated your own implicits, you must generate these samplings yourself.
We have provided an executable to do so:

```
./build/implicit_uniform implicit_surface num_samples
```
where:
- `implicit_surface` is the MPU surface file.
- `num_samples` is the number of samples to distribute on the surface. This should be a sufficiently large number to guarantee that all surface features are covered by the sampling. Of course this depends on the shape, so please use best judgement in determining a sufficient number of samples.

Evaluation may then be performed as follows:
```
./build/run_evaluation reconstructed_mesh implicit_surface dense_sampling output_base write_correspondences
```
where:
- `reconstructed_mesh` is the mesh output from the reconstruction algorithm.
- `implicit_surface` is the MPU surface file.
- `dense_sampling` is the dense uniformly sampled point cloud.
- `output_base` is the base file from which the reconstruction results will be written to. For instance, if 'results' is specified, then results.dist, results.recon, and optionally results.i2m and results.m2i will be output.
- `write_correspondences` is a flag indicating whether or not (1 or 0) the implicit to mesh and mesh to implicit point correspondences are to be written out (the .i2m and .m2i files).

The .dist file contains the individual distributions of the positional and normal error metrics: min, lower quartile,
median, upper quartile, max, and mean. The .recon file contains topological information about the mesh, see
`evaluator/GlobalStats.cpp` for more information. The m2i and i2m files may be read in via
`evaluator/ShortestDistanceMap.cpp`.

We suggest running evaluation on the surfaces produced via Poisson Surface Reconstruction. Please see
`modeling/models/bumps` for the reference point cloud produced via the particle system.


## Plotting Results
From the .dist file(s) generated through evaluation, we allow for two different options in plotting the results.
To generate a distribution over a single point cloud:
```
./build/single_distribution dist_file output_base
```
where:
- `dist_file` is the .dist file generated through `build/run_evaluation`.
- `output_base` is the base file name from which two plots, in pdf, will be generated. One is a box plot of the positional error distribution, and the other is a box plot of the normal error distribution. This is with respect to a single point cloud.

To generate distribution plots over a collection of point clouds:
```
./build/aggregate_distribution dist_base num_pcs output_base
```
where:
- `dist_base` is the base file name which the .dist files over all reconstruction evaluations reside. They must be numbered as dist_base_0.dist, dist_base_1.dist, ... dist_base_(num_pcs-1).dist.
- `num_pcs` is the number of point clouds over the distribution.
- `output_base` is the base file name from which four plots, in pdf, will be generated. These are mean distance distribution, Hausdorff distance distribution, mean normal deviation distribution, and max normal deviation distribution.

We suggest running the plotting executables on the running example, both individual and aggregate distributions,
to get a feel for the plotting.
