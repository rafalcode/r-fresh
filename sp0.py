#!/usr/bin/env python3
from simuPOP import *

pop = Population(size=1000, loci=[2])
pop.evolve(
    initOps = [
        InitSex(),
        InitGenotype(genotype=[1, 2, 2, 1])],  
    matingScheme=RandomMating(ops=Recombinator(rates=0.01)),
    postOps = [
        Stat(LD=[0, 1]),
        PyEval(r"'%.2f\n' % LD[0][1]", step=10),
    ],
    gen=100
)
