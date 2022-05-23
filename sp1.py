#!/usr/bin/env python3
import simuPOP as sim

pop = sim.Population(size=16, loci=5)
# above: ploidy=2 .. the default is diploid, so not necessary.
sim.initGenotype(pop, freq=[0.3, 0.7])
sim.dump(pop)
