# Frontend and Groth16 system

Run the algorithms in Groth16 ZK backend.

## Usage

Make sure you have cloned and update the submodule already.

First install essential packages:

``` bash
sudo apt-get install build-essential cmake git libgmp3-dev libprocps4-dev python-markdown libboost-all-dev libssl-dev
```

Then install `flatbuffer` through these cmds:

``` bash
git clone https://github.com/google/flatbuffers.git
# under `flatbuffers` folder
cmake -G "Unix Makefiles"
make
make install
```

Now build the repo using:

``` bash
# under `groth16` folder
$ mkdir build && cd build
# under `build` folder
$ cmake ..
$ make
```

Run the algotithms using given datasets, e.g.:

``` bash
# still under `build` folder
$ ./src/zkTI ../../datasets/d_duck_identification/answer.csv ../../datasets/d_duck_identification/truth.csv 30 10
```

The program will run Groth16 for algorithms and export the constraint system into .zkif files. Then use the files for Spartan system.