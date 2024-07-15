"""
genome2mps : decomposing nucleotide reads into matrix product states
"""
import numpy as np
import tensornetwork as tn
import networkx as nx
import matplotlib.pyplot as plt


def nuc_tensor(base, bond_dim, left=False, right=False, name="node"):
    """
    Return a tensor for each nucleotide base

    Parameters:
    -----------
     base : string
         Nucleotide base, one of A,T,C,G
     bond_dim : int
         Bond dimension to link MPS nodes
     left : bool, optional, default False
         True if the current base is the first in the sequence
     right : bool, optional, default False
         True if the current base is the last in the sequence
     name : string, optional, default "node"
         Name to assign the node

    Returns:
    --------
     tn.Node(tensor, name=name) : TN Node object
         The nucleotide base encoded as an MPS node
    """

    # Initialise a tensor of correct shape for left, right, and centre
    #
    # Left node:   |     Centre node:   |      Right node:    |
    #              o-                  -o-                   -o
    #
    # Edge dimesnions:
    # - has dimension `bond_dim`
    # | has dimension 2, the number of qubits needed for the alphabet {A,T,C,G}
    tensor = np.zeros((2, bond_dim)) if left else (np.zeros((bond_dim, 2)) if right else np.zeros((bond_dim, 2, bond_dim)))

    # Map nucleotides to their positions for a unique map
    base_map = {'A': [0,0], 'C': [0,1], 'G': [1,0], 'T': [1,1]}
    pos = base_map[base]
    
    # Set the appropriate position depending on node location
    # Logic: Generalise a unique encoding for each dimension of bond_dim
    #        such that the 2x2 structure is propagated
    for idx in range(bond_dim):
        if left:
            tensor[:, idx] = pos
        elif right:
            tensor[idx, :] = pos
        else:
            tensor[idx, :, idx] = pos
    
    return tn.Node(tensor, name = name)

def mps_from_seq(seq, bond_dim):
    """
    Return a MPS for a sequence of nucleotide bases

    Parameters:
    -----------
     seq : string
         Nucleotide sequence
     bond_dim : int
         Bond dimension to link MPS nodes

    Returns:
    --------
     mps : List
         List of TN Node objects in the MPS
     connected_edges : List
         List of TN Edge objects displaying connections in the MPS
    """
    mps = (
        [nuc_tensor(seq[0], bond_dim, left=True, name="node 0")] + 
        [
            nuc_tensor(nuc, bond_dim, name=f"node {idx + 1}") 
            for idx, nuc in enumerate(seq[1:-1])
        ] + 
        [nuc_tensor(seq[-1], bond_dim, right=True, name="end node")]
    )

    connected_edges = []
    conn = mps[0][1] ^ mps[1][0]
    connected_edges.append(conn)
    for idx in range(1,len(seq)-1):
        conn = mps[idx][2] ^ mps[idx+1][0]
        connected_edges.append(conn)

    return mps, connected_edges 

def main():
    """
    Main method to demonstrate functionality.
    Alter variables to show uniqueness and generality of encoding.
    """
    # Create tensors for the nucleotide sequence
    sequence = 'ACTAGTCGGT'
    bond_dimension = 2
    mps, edges = mps_from_seq(sequence, bond_dimension)

    for idx, node in enumerate(mps):
        print(node.tensor)
        if idx < len(edges):
            print(edges[idx])

###################################################################
if __name__ == """__main__""":
   main()
