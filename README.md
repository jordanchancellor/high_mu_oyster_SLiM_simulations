# SLiM-Oysters-Neutral

//neutral test script for oyster

//burn-in
initialize()
{

//non WF model type:
initializeSLiMModelType("nonWF");

//define life table constant for high mortality in junvenile stage
//set up life table assuming oysters won’t live past six years in an aquaculture-type environment 
//constant L and probability of mortality at each age
defineConstant("L", c(0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0));

//define carrying capacity
defineConstant("K", 1000);

//define output path 
defineConstant("output", "output/");

//overall mutation rate per base position per generation across the whole chromosome from Nate
initializeMutationRate(1e-5);

//mutation type “m1”: neutral mutation (0.5), fixed distribution of fitness effects given by 0.0 for “f”
initializeMutationType("m1", 0.5, "f", 0.0);

//sets mutations to convert to substitutions automatically becuase all mutations are completely neutral
m1.convertToSubstitution=T;

//genomic element types: specifies type of chromosomal region
initializeGenomicElementType("g1", m1, 1.0);

//genomic element type from above, that stretches from base position 0 to 99999
//chromosome with length 100,000 base pairs 
initializeGenomicElement(g1, 0, 99999);

//recombination rate
initializeRecombinationRate(1e-8);

//simulate separate sexes with a sex ratio default of 0.5
initializeSex("A");

}

//reproduction callbacks are used for each individual in the model for the generation of offspring in nonWF models
//called once per female, at the beginning of each generation
reproduction(p1, "F")
{

//subpopulation methods requesting new offspring with bi-parental mating 
//will draw one random individual as a mate for individuals older than age 1
//if (individual.age > 1) {
//mate = subpop.sampleIndividuals(1, replace = T, sex = "M", minAge=1); 
//generate litters drawing from a Poisson distribution with a mean of one million
litterSize = rpois(1, 10e6);
//add new offspring to the population

if (individual.age > 1) {
subpop.addCrossed(individual, subpop.sampleIndividuals(1, replace = T, sex = "M", minAge=1));
//disable callback for this generation
self.active =0;
	}
}

//execution of early Eidos events: viability/survival selection
1 early() 
{

//set initial population size, after generation1 will vary over time by births/deaths
//with random age from discrete uniform distribution from ages 0 to 6
sim.addSubpop("p1", 1000);
p1.individuals.age = rdunif(1000, min=0, max=6);

}

early()
{
//life table based age-related mortality
inds = p1.individuals;
ages = inds.age;
mortality = L[ages];
survival = 1 - mortality;
//survival-based fitness effect for each individual based on age
inds.fitnessScaling = survival;

//density-dependent mortality factoring in survival rate due to age-related mortality
p1.fitnessScaling = K / (p1.individualCount * mean(survival));

}

//execution of late Eidos events: output events reflecting the final state at the end of a generation
1 late()

{
sim.outputFull();
}


50 late()

{
sim.outputFull();
}

100 late()

{
sim.outputFull();
}

150 late()

{
sim.outputFull();
}

200 late()

{
sim.outputFull();
}
