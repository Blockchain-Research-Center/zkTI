# Frontend and Groth16 system

Run the algorithms in Groth16 ZK backend.

## Usage

Make sure you have cloned and update the submodule already.

Build the repo using:

``` bash
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