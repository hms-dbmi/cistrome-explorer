
```
snakemake --cores 2 --config filetype=mv5
# or
snakemake --cores 2 --config filetype=zarr
```

## Using pybbi with `summary="sum"`

Until version `0.2.3` is pushed to PyPI, install from source to get the changes added in pull request https://github.com/nvictus/pybbi/pull/12

```sh
git clone git@github.com:nvictus/pybbi.git
cd pybbi
pip install -e .
```

## Using parallel hdf5 via h5py and mpi4py

https://docs.h5py.org/en/latest/build.html#building-against-parallel-hdf5

### On O2

```sh
module load gcc/6.2.0
module load openmpi/3.1.0
module load hdf5/1.12.0

which mpicc
which h5cc # doesn't work for some reason - but the hdf5 dir is /n/app/hdf5/1.12.0.parallel

# CC="mpicc" HDF5_MPI="ON" HDF5_DIR=/n/app/hdf5/1.12.0.parallel pip install --no-binary=h5py h5py # doesn't work since pip h5py not compatible with 1.12.0

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
>>> h5py.get_config().mpi # Check whether MPI has been enabled, should return True
```