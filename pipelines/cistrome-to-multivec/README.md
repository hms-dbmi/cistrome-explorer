
# pipelines/cistrome-to-multivec

Generate multivec files from CistromeDB bigWig files.

```sh
conda activate cistrome-to-multivec-pipeline
# to generate multivec outputs
snakemake --cores 2 --config filetype=mv5
# or, to generate zarr outputs
snakemake --cores 2 --config filetype=zarr
# or, if on O2 (replace with your username)
./submit.sh mv5 my_username
```

# Setup

## Conda environment

```sh
conda env create -f environment.yml
conda activate cistrome-to-multivec-pipeline
```

## Snakemake cluster config

```sh
mkdir -p ~/.config/snakemake/cistrome-explorer
cp ./cluster-profile.yml ~/.config/snakemake/cistrome-explorer/config.yaml
```


## Using parallel hdf5 via h5py and mpi4py

*The following info is outdated, since h5py does not yet work with the parallel version of hdf5 installed on the o2 cluster. In the meantime we can do parallelization by submitting many simultaneous snakemake jobs for each output bigwig file.*

https://docs.h5py.org/en/latest/build.html#building-against-parallel-hdf5

### On O2

```sh
module load gcc/6.2.0
module load openmpi/3.1.0
module load hdf5/1.12.0

which mpicc
which h5cc # doesn't work for some reason - but the hdf5 dir is /n/app/hdf5/1.12.0.parallel

# CC="mpicc" HDF5_MPI="ON" HDF5_DIR=/n/app/hdf5/1.12.0.parallel pip install --no-binary=h5py h5py # doesn't work since pip h5py not compatible with 1.12.0

cd path/to/h5py-parent
# Clone h5py so that the latest code with support for hdf5 v1.12.0 (since not yet on pip).
git clone git@github.com:h5py/h5py.git

cd h5py
python setup.py configure --hdf5=/n/app/hdf5/1.12.0.parallel
python setup.py configure --mpi
python setup.py install
```

### On macOS

Download hdf5 1.10.6 source code from https://www.hdfgroup.org/downloads/hdf5/source-code/ and un-tar-gz

```sh
brew install openmpi

# Make a directory in which hdf5 can be installed.
mkdir -p ~/software/hdf5

brew info openmpi # Use this to find the CC value for the next line.

# In the downloaded hdf5-1.10.6 source directory:
CC=/usr/local/Cellar/open-mpi/4.0.3/bin/mpicc ./configure --enable-parallel --enable-shared --prefix=$HOME/software/hdf5
make
export NPROCS=3 # https://github.com/open-mpi/ompi/issues/6497
make check
make install

cd path/to/h5py-parent
# Clone h5py so that the latest code with support for hdf5 v1.12.0 (since not yet on pip).
git clone git@github.com:h5py/h5py.git

cd h5py
export CC=/usr/local/Cellar/open-mpi/4.0.3/bin/mpicc
python setup.py configure --hdf5=$HOME/software/hdf5
python setup.py configure --mpi
python setup.py install
cd ..
python
>>> import h5py
>>> h5py.get_config().mpi # Should return True if MPI has been enabled
```
