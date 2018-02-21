# Simulation for Gibson, Stoy and Lively 2018 "Bloody-minded parasites and sex: the effects of fluctuating virulence"

# based upon 2010 Epidemiological Model of Virulence, CM Lively
# Addition of variation in parasite virulence
# via variable birth rate of infected individuals

# This is the NBAR simulation, to be used in cases of negative autocorrelation (rho<0)
# PBAR from McKenzie 1985 Management Science 31

binomnbar<- function(bu,alpha,beta,rho){ 
  
  #  THE FUNCTION
  # these are the parameters you'll need to specify when you call this function:
  # bu  = number of offspring made by uninfected females
  # alpha = shape parameter for beta distribution to determine modification of bu with infection
  # beta = shape parameter for beta distribution to determine modification of bu with infection
  # rho = temporal autocorrelation coefficient, rho<0, rho>=-1
  
  # NUMBER OF INFECTED OFFSPRING
  # The number of offspring made by infected females (bi) is some fraction that made by uninfected females
  # The ratio bi/bu is temporally stochastic and is drawn from a beta distribution
  # The beta distribution is defined by alpha, beta, and rho
  
  # McKenzie's PBAR function
  # We establish two beta distribution, Q and W, with 2100 sampled values
  # See McKenzie 1985 for more information
  p= (rho*(alpha+beta))/(rho-1) # p is just a parameter derived from alpha, beta and rho; makes things cleaner
  Q=rbeta(2100, shape1=alpha, shape2=beta-p, ncp = 0)
  W=rbeta(2100, shape1=p, shape2=alpha-p, ncp = 0) 
  
  # Q and W are combined, resulting in a negatively autocorrelated beta distribution bounded from 0 to 1
  # The joint distribution will be stored as x
  
  x<- rep(0,2100)
  x[1] = 1
  # We then use this distribution to define the per-trial probability
  # of a binomial distribution, resulting in autocorrelated, beta-binomially
  # distributed numbers of offspring for infected moms
  
  bi<-rep(0,2100)
  bi[1] = rbinom(n=1,size=10,prob=x[1]) 
  # this is a way of showing what we'll do in the simulation 
  # we're sampling 1 value (n=1) from a binomial distribution
  # This value is the number of successes that we get in 10 trials (size=10), 
  # with a per-trial probability of success equal to x[1] (prob=x[1]). In this case, x[1]=1, so
  # we end up with infected females making 10 babies in the first generation, because the
  # probability of successfully having a baby in 10 tries is 100% per try 
  
  # OTHER PARAMETERS
  # Reproduction is density dependent: au and ai tell us to the degree to which offspring number is sensitive to population density
  # you'll see later that the ultimate number of offspring (Wu and Wi) is the result of bu,bi and au,ai and N
  au<-0.0001
  ai<-0.0001
  # here au=ai
  # if au and ai were different, infected individuals would differ in their sensitive to density
  
  
  # transmission parameter: number of propagules produced by each infection
  # I'm calling this z rather than Beta (term used in Lively 2010) because beta is being used in the beta distribution
  z<-10
  
  # cost of males 
  # this is equivalent to the sex ratio - the proportion of a female's offspring that are male. 
  # Here, we have an equal sex ratio - a mom makes 50% sons
  # the whole fitness function of sexuals gets multiplied by 1-s to give us the number of sexual females in the next generation
  s<-0.500
  
  # recombination - can range from 0 to 0.5
  r<-0.2
  
  # probability of migration of one infected individual
  # 2% chance of an infected migrant
  mi <- 0.02
  # probability of migration of one uninfected individual
  # 10% chance of an uninfected migrant
  mu <- 0.1
  
  # INITIATING POPULATION VECTORS
  # Here, we establish a bunch of empty vectors using rep(0,2100) with space to fill 2100 generations' worth of data
  # Data entered in these vectors will be stored. 
  
  # The population follows this sequence: allele frequencies to genotype frequencies; numbers of individuals per genotype and reproductive mode;
  # selection event: determination of infection rates, then realized birth rates, then reproduction of infected and uninfected individuals
  # then there's a recombination event in the sexual subpopulation.
  # We start all this outside the loop by establishing the population in the first generation
  
  # Population size - establish an empty vector that you'll fill in the loop below
  N<-rep(0,2100)
  
  # Starting population size
  N[1] = 8000
  
  # 1: allele frequencies
  # three alleles at two loci each (haploid)
  # locus 1
  A<-rep(0,2100)
  B<-rep(0,2100)
  C<-rep(0,2100)
  # locus 2
  X<-rep(0,2100)
  Y<-rep(0,2100)
  Z<-rep(0,2100)
  
  # Allele frequencies 
  # for locus 1 ( random deviates from 1/3 frequency.)
  A[1]<-1/3-(1/2-runif(1))/6
  B[1]<-1/3-(1/2-runif(1))/6
  C[1]<-1-A[1]-B[1]
  # locus 2
  X[1]<-1/3-(1/2-runif(1))/6
  Y[1]<-1/3-(1/2-runif(1))/6
  Z[1]<-1-X[1]-Y[1]
  
  # 2: genotype frequencies
  # sexual genotype frequencies 
  AX<-rep(0,2100)
  BX<-rep(0,2100)
  CX<-rep(0,2100)
  
  AY<-rep(0,2100)
  BY<-rep(0,2100)
  CY<-rep(0,2100)
  
  AZ<-rep(0,2100)
  BZ<-rep(0,2100)
  CZ<-rep(0,2100)
  
  # 3: numbers of individuals per genotypes
  # Total number of individuals in each genotype 
  # using the allele frequencies for locus 1 - ABC - and locus 2 - XYZ - above
  AX[1]<-A[1]*X[1]*N[1]
  BX[1]<-B[1]*X[1]*N[1]
  CX[1]<-C[1]*X[1]*N[1]
  
  AY[1]<-A[1]*Y[1]*N[1]
  BY[1]<-B[1]*Y[1]*N[1]
  CY[1]<-C[1]*Y[1]*N[1]
  
  AZ[1]<-A[1]*Z[1]*N[1]
  BZ[1]<-B[1]*Z[1]*N[1]
  CZ[1]<-C[1]*Z[1]*N[1]
  
  # and for each reproductive modes 
  # SEX: the total number of sexual individuals is the sum of numbers of each sexual genotypes
  sex<-rep(0,2100)
  sex[1]<-AX[1]+AY[1]+AZ[1]+BX[1]+BY[1]+BZ[1]+CX[1]+CY[1]+CZ[1]
  
  # ASEX: the number of asexuals is determined within the simulation 
  # the clone will be introduced at generation 1000
  asex<-rep(0,2100)
  
  # 4: selection event - Infection 
  # establish vectors for numbers of infected and uninfected individuals
  AXi<-rep(0,2100)
  BXi<-rep(0,2100)
  CXi<-rep(0,2100)
  
  AYi<-rep(0,2100)
  BYi<-rep(0,2100)
  CYi<-rep(0,2100)
  
  AZi<-rep(0,2100)
  BZi<-rep(0,2100)
  CZi<-rep(0,2100)
  
  # total sexual infected
  sexi<-rep(0,2100)
  
  # total asexual infected 
  asexi<-rep(0,2100)
  
  # uninfected individuals
  AXu<-rep(0,2100)
  BXu<-rep(0,2100)
  CXu<-rep(0,2100)
  
  AYu<-rep(0,2100)
  BYu<-rep(0,2100)
  CYu<-rep(0,2100)
  
  AZu<-rep(0,2100)
  BZu<-rep(0,2100)
  CZu<-rep(0,2100)
  
  # total asexual uninfected 
  asexu<-rep(0,2100)
  
  #Starting number of infected individuals 
  AXi[1]<-1
  BXi[1]<-1
  CXi[1]<-1
  
  AYi[1]<-1
  BYi[1]<-1
  CYi[1]<-1
  
  AZi[1]<-1
  BZi[1]<-1
  CZi[1]<-1
  
  # Uninfected individual numbers
  # This is the total number of individuals in a genotype minus the number of infected
  AXu[1]<-AX[1]-AXi[1]
  BXu[1]<-BX[1]-BXi[1]
  CXu[1]<-CX[1]-CXi[1]
  
  AYu[1]<-AY[1]-AYi[1]
  BYu[1]<-BY[1]-BYi[1]
  CYu[1]<-CY[1]-CYi[1]
  
  AZu[1]<-AZ[1]-AZi[1]
  BZu[1]<-BZ[1]-BZi[1]
  CZu[1]<-CZ[1]-CZi[1]
  
  # 5: selection event: determine realized birth rates, or fitness, based on density and infection rates
  # Realized birth rate calculation: use intrinsic birth rates corrected by population size and density sensitivity
  Wu<-rep(0,2100)
  Wi<-rep(0,2100)
  
  # birth rates at gen 1
  Wu[1]<-bu/(1+(au*N[1]))
  Wi[1]<-bi[1]/(1+(ai*N[1])) 
  
  # these values define virulence 
  # calculated as the proportional decrease in birth rate due to infection
  V<-rep(0,2100)
  V[1]<-(Wu[1]-Wi[1])/Wu[1]
  
  # 6: selection event - reproduction 
  # These are the product of the selection event (post-selection) - number of individuals born to each genotype
  AXps<-rep(0,2100)
  AYps<-rep(0,2100)
  AZps<-rep(0,2100)
  
  BXps<-rep(0,2100)
  BYps<-rep(0,2100)
  BZps<-rep(0,2100)
  
  CXps<-rep(0,2100)
  CYps<-rep(0,2100)
  CZps<-rep(0,2100)
  
  # starting number offspring  - selection
  # sexuals
  # cost of males multiplied by mean fitness of sexuals
  # mean fitness is number of parents * fitness (offspring number), with differentation according to infected or not
  AXps[1]<-(1-s)*((AXu[1]*Wu[1])+(AXi[1]*Wi[1]))
  AYps[1]<-(1-s)*((AYu[1]*Wu[1])+(AYi[1]*Wi[1]))
  AZps[1]<-(1-s)*((AZu[1]*Wu[1])+(AZi[1]*Wi[1]))
  
  BXps[1]<-(1-s)*((BXu[1]*Wu[1])+(BXi[1]*Wi[1]))
  BYps[1]<-(1-s)*((BYu[1]*Wu[1])+(BYi[1]*Wi[1]))
  BZps[1]<-(1-s)*((BZu[1]*Wu[1])+(BZi[1]*Wi[1]))
  
  CXps[1]<-(1-s)*((CXu[1]*Wu[1])+(CXi[1]*Wi[1]))
  CYps[1]<-(1-s)*((CYu[1]*Wu[1])+(CYi[1]*Wi[1]))
  CZps[1]<-(1-s)*((CZu[1]*Wu[1])+(CZi[1]*Wi[1]))
  
  # total number of sexual offspring
  sexps<-rep(0,2100)
  sexps[1]<-AXps[1]+AYps[1]+AZps[1]+BXps[1]+BYps[1]+BZps[1]+CXps[1]+CYps[1]+CZps[1]
  
  # asexual vector
  # asexuals enter the simulation 1000 gens in
  # 0 until gen 999
  asexps<-rep(0,2100)
  asexps[999]<-1
  
  # 7: recombination
  # frequencies of genotypes post-selection
  fAXps<-rep(0,2100)
  fBXps<-rep(0,2100)
  fCXps<-rep(0,2100)
  
  fAYps<-rep(0,2100)
  fBYps<-rep(0,2100)
  fCYps<-rep(0,2100)
  
  fAZps<-rep(0,2100)
  fBZps<-rep(0,2100)
  fCZps<-rep(0,2100)
  
  # genotype frequency within sexual population
  # necessary for recombination calculations
  fAXps[1]<-AXps[1]/sexps[1]
  fBXps[1]<-BXps[1]/sexps[1]
  fCXps[1]<-CXps[1]/sexps[1]
  
  fAYps[1]<-AYps[1]/sexps[1]
  fBYps[1]<-BYps[1]/sexps[1]
  fCYps[1]<-CYps[1]/sexps[1]
  
  fAZps[1]<-AZps[1]/sexps[1]
  fBZps[1]<-BZps[1]/sexps[1]
  fCZps[1]<-CZps[1]/sexps[1]
  
  
  ######################################################        
  # SIMULATION
  gen<-rep(0,2100)
  gen[1]<-1
  
  # the loop
  for(i in 2:2100) {
    gen[i]<-i
    
    # 1: allele frequencies - we're connecting step 6 (selection event) to step 1 by using the genotype frequencies post-selection to calculate allele frequencies
    # based upon results of first generation
    # initiates generation 2
    
    # allele frequencies
    A[i]<-(AXps[i-1]+AYps[i-1]+AZps[i-1])/sexps[i-1]
    B[i]<-(BXps[i-1]+BYps[i-1]+BZps[i-1])/sexps[i-1]
    C[i]<-(CXps[i-1]+CYps[i-1]+CZps[i-1])/sexps[i-1]
    
    X[i]<-(AXps[i-1]+BXps[i-1]+CXps[i-1])/sexps[i-1]
    Y[i]<-(AYps[i-1]+BYps[i-1]+CYps[i-1])/sexps[i-1]
    Z[i]<-(AZps[i-1]+BZps[i-1]+CZps[i-1])/sexps[i-1]
    
    # 2: genotype frequencies - again, connecting step 7 (recombination) to step 2 using the recombination equation to get our new genotype frequencies
    
    # proportion (1-r) stay unrecombined, while proportion r recombine
    AXr<-(1-r)*fAXps[i-1]+(r*A[i]*X[i])
    BXr<-(1-r)*fBXps[i-1]+(r*B[i]*X[i])
    CXr<-(1-r)*fCXps[i-1]+(r*C[i]*X[i])
    
    AYr<-(1-r)*fAYps[i-1]+(r*A[i]*Y[i])
    BYr<-(1-r)*fBYps[i-1]+(r*B[i]*Y[i])
    CYr<-(1-r)*fCYps[i-1]+(r*C[i]*Y[i])
    
    AZr<-(1-r)*fAZps[i-1]+(r*A[i]*Z[i])
    BZr<-(1-r)*fBZps[i-1]+(r*B[i]*Z[i])
    CZr<-(1-r)*fCZps[i-1]+(r*C[i]*Z[i])
    
    # 3: numbers of individuals per genotype in next generation, following recombinations
    # genotype numbers post-recombination (frequency * # sexuals)
    # Calculate numbers by using genotype frequencies multiplied by known number of sexual offspring post-selection
    AX[i]<-AXr*sexps[i-1]
    BX[i]<-BXr*sexps[i-1]
    CX[i]<-CXr*sexps[i-1]
    
    AY[i]<-AYr*sexps[i-1]
    BY[i]<-BYr*sexps[i-1]
    CY[i]<-CYr*sexps[i-1]
    
    AZ[i]<-AZr*sexps[i-1]
    BZ[i]<-BZr*sexps[i-1]
    CZ[i]<-CZr*sexps[i-1]
    
    # Clone numbers
    # clone doesn't enter until generation 1000
    # at generation 999, 1 clone of genotype AX has appeared post-selection, giving us 1 clone in asex[] vector
    if(i<1000){
    }else {
      asex[i]<-asexps[i-1]
    }
    
    # population size
    # sexual population size, post-recombination
    sex[i]<-AX[i]+AY[i]+AZ[i]+BX[i]+BY[i]+BZ[i]+CX[i]+CY[i]+CZ[i]
    # total population size, sexual plus clone
    N[i]<-AX[i]+AY[i]+AZ[i]+BX[i]+BY[i]+BZ[i]+CX[i]+CY[i]+CZ[i]+asex[i]
    
    # 4 - selection event. Infection
    # number of infected (i) and uninfected individuals (u)
    # the term 
    #number per genotype/total population size = poisson mean number of exposures of hosts and parasites with matching genotypes
    # exp(): zero class, i.e. probability of not encountering a parasite; 1-exp(), probability of encountering 1 or more parasites
    # stochastic migration term
    
    AXi[i]<-AX[i]*(1-exp(-z*(AXi[i-1]+asexi[i-1])/N[i]))+ if (runif(1)<mi) # probability of encountering a parasites
    {1} else {0} # ifelse term gives us stochastic migration term
    AXu[i]<-AX[i]*(exp(-z*(AXi[i-1]+asexi[i-1])/N[i]))+ if (runif(1)<mu) # probability of not encountering a parasite
    {1} else {0}
    
    AYi[i]<-AY[i]*(1-exp(-z*AYi[i-1]/N[i]))+ if (runif(1)<mi)
    {1} else {0}
    AYu[i]<-AY[i]*(exp(-z*AYi[i-1]/N[i]))+ if (runif(1)<mu)
    {1} else {0}
    
    AZi[i]<-AZ[i]*(1-exp(-z*AZi[i-1]/N[i]))+ if (runif(1)<mi)
    {1} else {0}
    AZu[i]<-AZ[i]*(exp(-z*AZi[i-1]/N[i]))+ if (runif(1)<mu)
    {1} else {0}
    
    BXi[i]<-BX[i]*(1-exp(-z*BXi[i-1]/N[i]))+ if (runif(1)<mi)
    {1} else {0}
    BXu[i]<-BX[i]*(exp(-z*BXi[i-1]/N[i]))+ if (runif(1)<mu)
    {1} else {0}
    
    BYi[i]<-BY[i]*(1-exp(-z*BYi[i-1]/N[i]))+ if (runif(1)<mi)
    {1} else {0}
    BYu[i]<-BY[i]*(exp(-z*BYi[i-1]/N[i]))+ if (runif(1)<mu)
    {1} else {0}
    
    BZi[i]<-BZ[i]*(1-exp(-z*BZi[i-1]/N[i]))+ if (runif(1)<mi)
    {1} else {0}
    BZu[i]<-BZ[i]*(exp(-z*BZi[i-1]/N[i]))+ if (runif(1)<mu)
    {1} else {0}
    
    CXi[i]<-CX[i]*(1-exp(-z*CXi[i-1]/N[i]))+ if (runif(1)<mi)
    {1} else {0}
    CXu[i]<-CX[i]*(exp(-z*CXi[i-1]/N[i]))+ if (runif(1)<mu)
    {1} else {0}
    
    CYi[i]<-CY[i]*(1-exp(-z*CYi[i-1]/N[i]))+ if (runif(1)<mi)
    {1} else {0}
    CYu[i]<-CY[i]*(exp(-z*CYi[i-1]/N[i]))+ if (runif(1)<mu)
    {1} else {0}
    
    CZi[i]<-CZ[i]*(1-exp(-z*CZi[i-1]/N[i]))+ if (runif(1)<mi)
    {1} else {0}
    CZu[i]<-CZ[i]*(exp(-z*CZi[i-1]/N[i]))+ if (runif(1)<mu)
    {1} else {0}
    
    # asex infection and non-infection
    # operates only in generation 1000 and on
    if(i<1000){
    }else {
      asexi[i]<-asex[i]*(1-exp(-z*(AXi[i-1]+asexi[i-1])/N[i]))
      asexu[i]<-asex[i]*(exp(-z*(AXi[i-1]+asexi[i-1])/N[i]))
    }
    
    # total number of infected sexuals
    sexi[i]<-AXi[i]+AYi[i]+AZi[i]+BXi[i]+BYi[i]+BZi[i]+CXi[i]+CYi[i]+CZi[i]
    
    # 5 - selection event: realized birth rate
    # determine intrinsic offspring number for infected individuals in this generation
    # then we have to calculate realized number offspring (fitness), based upon population size at this time step.
    
    # draw your per-trial probability - bi/bu = x - for the 2100 generations of the simulation
    x[i]=1-(U[i]*(1-W[i]*x[i-1])) # this merging of U and W creates a beta distribution with shape parameters alpha and beta
    # with autocorrelation defined by rho
    bi[i]=rbinom(1,size=10,prob=x[i]) # then draw your infected birth number from a binomial distribution
    # 10 trials, per-trial probability equal to x
    
    # we include these ifelse statements as additional way to ensure that birth rate stays above zero
    if ((bu/(1+(au*N[i])))<0){
      Wu[i]<-0
    } else {
      Wu[i]<-(bu/(1+(au*N[i])))
    }
    
    if (bi[i]/(1+(ai*N[i]))<0){
      Wi[i]<-0
    } else {
      Wi[i]<-(bi[i]/(1+(ai*N[i])))
    } # note key difference from deterministic simuation - bi now changes with generation. bi is not a fixed parameter
    
    # these equations determine virulence
    V[i]<-(Wu[i]-Wi[i])/Wu[i]
    
    # 6 selection event: reproduction 
    AXps[i]<-(1-s)*((AXu[i]*Wu[i])+(AXi[i]*Wi[i])) # number of offspring produced by each genotype, based upon infection rate and fitness
    AYps[i]<-(1-s)*((AYu[i]*Wu[i])+(AYi[i]*Wi[i]))
    AZps[i]<-(1-s)*((AZu[i]*Wu[i])+(AZi[i]*Wi[i]))
    
    BXps[i]<-(1-s)*((BXu[i]*Wu[i])+(BXi[i]*Wi[i]))
    BYps[i]<-(1-s)*((BYu[i]*Wu[i])+(BYi[i]*Wi[i]))
    BZps[i]<-(1-s)*((BZu[i]*Wu[i])+(BZi[i]*Wi[i]))
    
    CXps[i]<-(1-s)*((CXu[i]*Wu[i])+(CXi[i]*Wi[i]))
    CYps[i]<-(1-s)*((CYu[i]*Wu[i])+(CYi[i]*Wi[i]))
    CZps[i]<-(1-s)*((CZu[i]*Wu[i])+(CZi[i]*Wi[i]))
    
    
    sexps[i]<-AXps[i]+AYps[i]+AZps[i]+BXps[i]+BYps[i]+BZps[i]+CXps[i]+CYps[i]+CZps[i]
    
    #asexuals
    if(i<1000){
    }else {
      asexps[i]<-((asexu[i]*Wu[i])+(asexi[i]*Wi[i]))
    }
    
    # 7 - recombination
    # make genotype frequencies pre-recombination
    fAXps[i]<-AXps[i]/sexps[i]
    fBXps[i]<-BXps[i]/sexps[i]
    fCXps[i]<-CXps[i]/sexps[i]
    
    fAYps[i]<-AYps[i]/sexps[i]
    fBYps[i]<-BYps[i]/sexps[i]
    fCYps[i]<-CYps[i]/sexps[i]
    
    fAZps[i]<-AZps[i]/sexps[i]
    fBZps[i]<-BZps[i]/sexps[i]
    fCZps[i]<-CZps[i]/sexps[i]
  }
  
  # some summary results
  par(mfrow=c(3,1))
  plot(gen[950:1200],sex[950:1200],type="l",col="red",xlab="Generation",ylab="Number",ylim=c
       (0,50000),xlim=c(950,1200))
  points(gen[950:1200],asex[950:1200],type="l",col="blue")
  plot(gen[950:1200],sexi[950:1200]/sex[950:1200],type="l",col="red",xlab="Generation",ylab="Proportion 
       Infected",ylim=c(0,1),xlim=c(950,1200))
  points(gen[950:1200],asexi[950:1200]/asex[950:1200],type="l",col="blue")
  plot(gen[950:1200],V[950:1200],type="l",xlab="Generation",ylab="Virulence",col="green",ylim=c(0,1), xlim=c(950,1200))
  
  Num_asex<-mean(asex[2000:2100])
  Num_sex<-mean(sex[2000:2100])
  Percent_sex<-mean(sex[2000:2100]/N[2000:2100])
  Virulence<-mean(V[2000:2100])
  vector<-c(Num_asex,Num_sex,Percent_sex,Virulence) # returns a vector with mean number of asexuals, sexuals, mean percent sex, and mean virulence
  # from generation 2000-2100
  return(vector)
  
}