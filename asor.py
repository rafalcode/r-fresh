def simuAssortativeMating(w, size, gen, vsp=[0, 4]):
    '''
        w       proportion of general random mating.
        size    population size
        gen     how many generation to run
        vsp     virtual subpopulations for assortative mating.
    '''
    pop = Population(size, loci=[1])
    # define four virtual subpopulations. individuals in the first three virtual
    # subpopulation.have genotype (0, 0), (0, 1) or (1, 0), and (1, 1) respectively,
    # and have at leat one mutant (allele 1) in the last virtual subpopulation.
    pop.setVirtualSplitter(GenotypeSplitter(loci=0,
        alleles=[[0, 0], [0, 1], [1, 1], [0, 0, 0, 1], [0, 1, 1, 1]]))

    pop.evolve(
        initOps = [
            InitSex(),
            InitGenotype(freq=[0.5, 0.5]),
            PyExec('AaNum=[]'),  # initialize a list in population's local dictionary
            ],
            # calculate virtual population sizes, and allele frequency at locus 0.
        preOps = Stat(popSize=True, alleleFreq=[0], subPops=[(0,0), (0,1), (0,2)]),
        # Negative weight means fixed size (weight * current subpopulation size).
        # In the case of no positive weight, zero weights means proportional to
        # parental (virtual) subpopulation size.
        matingScheme = HeteroMating([RandomMating(weight = -1*w),
            RandomMating(subPops=[(0, x) for x in vsp], weight = 0)]),
        postOps = [
            # print size of virtual populations and allele frequency
            PyEval(r"'#inds with genotype AA %4d, Aa %4d, aa %4d, freq of A: %.1f\n' % "
                "(subPopSize[0], subPopSize[1], subPopSize[2], alleleFreq[0][0]*100)"),
            # append number of individuals with genotype Aa to list AaNum
            PyExec(r"AaNum.append(subPopSize[1])")
        ],
        gen = gen
    )
    return pop.dvars().AaNum
