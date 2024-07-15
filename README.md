# Encoding Genome Reads as MPS

This is a simple script for encoding genome reads as matrix
product states (MPS) using the `TensorNetwork` and `numpy`
python packages.

## Install

To install the dependencies, run:

```bash
pip3 install -r requirements.txt
```

## Usage

We plan to add more functionality soon, this is just a proof-of-
concept. Modify the variables in the `main` method in
`genome2mps.py` and run with:

```bash
python3 genome2mps.py
```

This will print the tensors and the edges between them in the MPS

## TODO
[] Implement contractions
[] Port to a quantum simulator
[] Add more ?
