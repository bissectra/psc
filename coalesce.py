from numpy.random import default_rng

class Node:
    """
    A node in a gene genealogy. Has the following fields:
    left: pointer to left child, or None if there are no children
    right: same for right child
    parent: same for parent
    age: age of current node in generations
    """

    def __init__(self, left=None, right=None, age=0.0):
        self.parent = None
        self.left = left
        self.right = right
        self.age = age
        if self.left != None:
            self.left.parent = self
        if self.right != None:
            self.right.parent = self

    """
    Print self and all children
    """
    def print(self, indent=0):
        if self == None:
            return
        print(indent * "-", self.age, sep="")
        if self.left != None:
            self.left.print(indent+1)
            self.right.print(indent+1)

    """
    Length of gene genealogy: total length of all branches.
    """
    def length(self):
        if self.parent != None:
            brlen = self.parent.age - self.age
        else:
            brlen = 0.0
            
        if self.left != None:
            brlen += self.left.length()

        if self.right != None:
            brlen += self.right.length()

        return brlen
    """
    Return the number of descendants of a Node
    """
    def ndescendants(self):
        if self.left == None:
            # Tip node has 1 descendant: itself
            return 1
        else:
            return self.left.ndescendants() + self.right.ndescendants()

"""
One iteration of a coalescent simulation. This function allocates
an array of k gene copies, and then runs a coalescent simulation, with
diploid population size n, to link them into a gene genealogy. It
returns a pointer to the root Node. 
"""        
def coalesce(rng, n = 10000, k=5):
    age = 0.0
    fourN = 4*n

    # Array of nodes
    nodes = [Node() for i in range(k)]
    
    # Each pass through loop deals with one coalescent interval.
    while k > 1:
        # The mean time until the next coalescent event is the
        # inverse of the hazard.
        tmean = fourN/(k*(k-1))
        age += rng.exponential(tmean)  # time until next coalescent event

        # Choose 2 random nodes, ensuring that i < j.
        i = rng.integers(k)
        j = rng.integers(k-1)
        if j == i:
            j = k-1
        if j < i:
            i, j = j, i

        # Join the two nodes, putting their parent into position i.
        nodes[i] = Node(nodes[i], nodes[j], age)

        # Shorten array; forget about the original nodes[j]
        if j < k-1:
            nodes[j] = nodes[k-1]
        
        k -= 1

    # Return root node
    return nodes[0]

"""
Expected site frequency spectrum under selective neutrality and constant
population size. k is the haploid sample size, and theta = 4*N*mu,
where N is diploid population size and mu is the mutation rate per
generation. 
"""
def espectrum(theta, k=5):
    return [theta/i for i in range(1, k)]

"""
Mean squared difference between two lists. (This would be faster in
Numpy.)
"""
def msqdiff(a, b):
    n = len(a)
    assert n == len(b)
    return sum([(x-y)**2 for x, y in zip(a,b)])/n

def main(n = 10000, k=5):
    # Random number generator.
    # rng.poisson(3) # x is a sample from Poisson dist w/ mean 3
    # rng.integers(3) # x is a sample from {0, 1, 2}
    # rng.exponential(3) # x is a sample from an exponential w/ mean 3
    rng = default_rng()

    theta = 1         # theta = 4*n*mu. Keep this constant
    mu = theta/(4*n)  # mutation rate

    # expected spectrum
    espec = espectrum(theta, k)

    nreps = 10     # make this big after everything works
    mean_msq = 0.0 # mean squared deviation btw observed and expected
              
    for rep in range(nreps):
        # There are lots of print statements inside the loop. Get rid
        # of them once you get your code working.
        root = coalesce(rng, n, k)

        # For after writing your spectrum function:
        # spec = [0 for i in range(k-1)]
        # spectrum(root, rng, mu, spec)
        # msq = msqdiff(spec, espec)
        # mean_msq += msq

        print(f"Repetition {rep}")
        print(f"  tree depth         : {root.age}")
        print(f"  total branth length: {root.length()}")
        print(f"  number of tips     : {root.ndescendants()}")

        # For comparing observed and expected spectra
        # print("spectrum:", end="")
        # for i in range(k-1):
        #    print(f" {spec[i]}", end="")
        # print()
        # print("espectrum:", end="")
        # for i in range(k-1):
        #     print(f" {espec[i]}", end="")
        # print()
        
        # print(f"  msq deviation      : {msq}")

    # For after writing your spectrum function.
    # mean_msq /= nreps
    # print("mean msq:", mean_msq)
