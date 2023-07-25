# ZKTI

Zero-knowledge proofs (ZK) for truth inference (TI) algorithms.

## Clone

Clone the repository along with its submodules.

``` bash 
$ git clone --recursive https://github.com/Blockchain-Research-Center/zkTI.git
```

## Usage

We write the circuit of our scheme with libSNARK which supports the Groth16 system. The circuit can be imported to a Spartan system with our exporter.

To run this program, first to build the C++ code inside the $groth16$ folder, then use the Spartan system under $spartan$ folder.

Run the programs according to the README.md instructions inside each folder.