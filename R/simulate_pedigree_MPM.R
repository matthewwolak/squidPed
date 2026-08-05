

#' @aliases simulate_pedigree
#' @rdname simulate_pedigree
#' @param A,F,U Matrices representing the Population, Fertility, and Survival
#'   components of the matrix projection model
#' @param mat_sex character indicating the population projection matrices
#'   either characterize females specifically, or there are sex-specific matrices
#' @param psr_p_female numerical proportion of female at conception (primary sex-
#'   ratio)

#' @export


## what do we output - how do we incorporate alive individuals that don't breed or multiple measurements
## output pedigree and data structure?

	# years = 5
	# n_females = 10

	# p_breed = 1
	# juv_surv = 0.25
	# adult_surv = 0.5
	# immigration = 0
	# afr=1

	# fecundity = 4
	# p_sire = 1
	# p_retain = 0.8
	# polgyny_rate = 0

	# constant_pop = TRUE 
	# known_age_structure = FALSE

	# 		set.seed(85)
	# years = peds_param_reduced[k,"generations"]
	# n_females = peds_param_reduced[k,"n_females"]
	# fecundity = peds_param_reduced[k,"fecundity"]
	# fixed_fecundity = TRUE
	# p_sire = peds_param_reduced[k,"p_sire"]
	# p_polyandry=1
	# p_breed=1
	# juv_surv = c(peds_param_reduced[k,"juv_surv_f"],peds_param_reduced[k,"juv_surv_m"])
	# adult_surv = 0
	# immigration = c(peds_param_reduced[k,"immigration_f"],peds_param_reduced[k,"immigration_m"])
	# constant_pop = TRUE
	# 	p_retain = 0
	# afr=1
	# verbose=TRUE


##  todo
# sex specific rates
# juvenile survival

simulate_pedigree_MPM <- function(
        # A, F, U: NO defaults, therefore must be supplied to run (A = F + U)
	A,  #<-- matrix population model (a square projection matrix) 
	F,  #<-- square projection matrix reflecting transitions due to reproduction
	U,  #<-- square projection matrix reflecting survival-related transitions 
	mat_sex = c("Female", "sex-specific"),
	psr_p_female = 0.5, #<-- proportion female at conception (primary sex ratio)
	years = 5,
	n_females = 50, #<-- total number of sexually mature females to begin
#	afr=1,
#	p_breed = 1,
#	fecundity = 4,
	fixed_fecundity = TRUE,
	p_polyandry = 0,
	p_sire = 1,
	p_retain = 0,
#	juv_surv = 0.25,
#	adult_surv = 0.5,
#	immigration = 0,
	# polgyny_rate = 0,
#	constant_pop = TRUE ,
#	known_age_structure = FALSE,
	verbose = FALSE){

#XXX CHANGES denoted (KEY to markings) XXX
## [+++ +++] addition of entire chunks/lines of code between these lines
## [--- ---] deletion of entire chunks/lines of between these lines
## [||| |||] partial addition/deletion to existing code (a Change)


#[+++ generate generic variable used a lot
  ## Total number of stages
  stgn <- nrow(A)
#+++] 

#[+++ whether matrices female or sex-specific (allow partial matching)
##TODO See "StudiedSex" in com(p)adre database entry (e.g., `RIdbt$StudiedSex`)
## XXX Then determine how expand or not matrices and stage counts if sex-specific
### matrices (i.e., both matrices specified). Otherwise FORCE a default=50:50
#### (because if no information given about males then they cannot
#### have different U matrix and only way to have sex ratio != 50:50 requires
#### different U
#### also assume matrix=female specific=male-specific matrix except for psr_p_female
  mat_sex <- match.arg(mat_sex)
  if(mat_sex == "Female"){
    if(psr_p_female != 0.5){
      psr_p_female <- 0.5
      warning("Resetting primary sex ratio (psr_p_female) to 0.5. If matrices are specific to females then no information about males requires sex ratio to be 50:50\n")
    }
  }
#+++]

#[+++ new checks
  # Make sure matrices are OK
  ## re-assign row/column names of matrices to generic numbers
  ###XXX TODO THOUGH watch out for sex-specific names and how to handle sex-specific
  #### matrices
  dimnames(A) <- list(as.character(seq.int(stgn)), as.character(seq.int(stgn)))
  dimnames(F) <- list(as.character(seq.int(stgn)), as.character(seq.int(stgn)))
  dimnames(U) <- list(as.character(seq.int(stgn)), as.character(seq.int(stgn)))
  #TODO FIXME XXX what else?
#+++]  


  options(stringsAsFactors=FALSE) # as long as later version of R - dont need

	## get environment, for transform det function
	Renv <- environment()


	# check and generate all sex specific variables

#[--- No longer needed since matrix model
#	general_check("afr", env=Renv, rate=FALSE, sex_specific=TRUE)
#	lapply(c("p_breed", "adult_surv","juv_surv","immigration"),general_check, env=Renv, rate=TRUE, sex_specific=TRUE)
#---]


#[||| removed n_females since part of A matrix
## no longer need to lapply() since left with only 1 function eval
	general_check("years", env=Renv, rate=FALSE, sex_specific=FALSE)
#|||]

#	if(fixed_fecundity){
#		lapply("fecundity", general_check, env=Renv, rate=FALSE, sex_specific=FALSE)
#	}else{
#		if(x<0) stop("fecundity must be >=0")
#	}

	lapply(c("p_retain", "p_sire","p_polyandry"), general_check, env=Renv, rate=TRUE, sex_specific=FALSE)



#[---
#	det_growth_rate_f <- (p_breed_f * juv_surv_f * fecundity)/2 + adult_surv_f + immigration_f
	# male growth rate probably doesn't matter as long as all eggs are fertilised? although at some point would run out of males
	# det_growth_rate_m <- (p_breed_m *juv_surv_m * fecundity)/2 + adult_surv_m + immigration_m
#---]
#[+++ determine the growth rate, particularly to make sure number of years of data
## requested is likely possible:
#TODO need argument above about sex-specific matrix or not
## XXX currently ASSUMING female matrix as A
	det_growth_rate_f <- Re(eigen(A, symmetric = FALSE,
	                                      only.values = TRUE)$values[1]) 
#+++]
#[||| when growth rate >1 changed message from warning to just a printed message
	if(!isTRUE(all.equal(det_growth_rate_f,1))){
		if(det_growth_rate_f>1) cat("growth rate is more than 1\n") 
	  if(det_growth_rate_f<1) warning("growth rate is less than 1") 	
	}
#|||]	
	# v_as <- adult_surv * (1-adult_surv)
	# v_imm <- immigration * (1-immigration)
	# v_js <- juv_surv * (1-juv_surv)
	# v_fec <- fecundity
	
	# v_rec <- fecundity^2*v_js + juv_surv^2*v_fec + v_js*v_fec

	# stoch_growth_rate <- det_growth_rate - (v_as + v_rec + v_imm)/(2*n_females)

	# immigration <- 1- stoch_growth_rate

#[--- Can delete all below and not make starting population but ENTIRE pedigree
###
# MAKE STARTING POPULATION
###

	# make pedigree for base population
	# pedigree <- data.frame(
	# 	animal = paste0("0_",1:(n_females*2)),
	# 	dam = NA,
	# 	sire = NA,
	# 	sex = rep(c("F","M"),each=n_females),
	# 	## starting age structure for female and males, based on constant survival rate
	# 	## we could delete age/cohort of these individuals in the output as that would be realistic to a real pedigree?
	# 	cohort= -1 * c(
	# 		if(adult_surv_f==0){
	# 			rep(0,n_females)
	# 	  }else{
	# 		  rgeom(n_females,adult_surv_f)
	# 	  },
	# 	  if(adult_surv_m==0){
	# 	  	rep(0,n_females)
	# 	  }else{
	# 	  	rgeom(n_females,adult_surv_m)
	# 	  })
	# )

## this needs to be sex specific with JS
# 	yearly_recruits <- n_females*p_breed*fecundity*juv_surv/2
#	starting_n <- c(n_females,rep(yearly_recruits,afr-1))

#	pedigree <- data.frame(
#		animal = paste0("0_",1:sum(starting_n*2)),# or could code them with their cohort. using 0 means all founder are the same
#		dam = NA,
#		sire = NA,
		# equal sex ratio to start
#		sex = rep(c("F","M"),sum(starting_n)),
		## starting age structure for female and males, based on constant survival rate
		## we could delete age/cohort of these individuals in the output as that would be realistic to a real pedigree?
#		cohort= rep(-1*(afr-1):0,c(n_females,rep(yearly_recruits,afr-1))*2)
		## need to start sims with some age structure, otherwise with afr>1 the population size will drop after first year, as new offspring wont have recruited yet. So this way has founder that are recruited for the first afr-1 years. the new number of recruits per year is n_females*p_breed*fecundity*juv_surv for each sex
#	)
## could make some stochasticity in the starting number if constant_pop=FALSE

#---]





#[+++  determine number of individuals each age/stage for entire simulation
         # create generic function to do this
         
#########################################################
# Define FUNCTION to determine number of individuals in each stage each timestep
## do generically so can use same code/function for females and males

#TODO/FIXME: if mat_sex != "Female" then alter below

stg_sim_n <- function(A, F, n, ...){
  #TODO: should population be modeled at stable stage/age distribution?
  # Seed population count vector  #<--TODO better/more official name for this
  ## Determine first stage/age/index at reproduction
  stgfr <- min(which(F[1, ] > 0))
  # Create matrix of stage abundances for each timestep
  nmat <- matrix(NA, nrow = stgn, ncol = years)
    dimnames(nmat) <- list(dimnames(A)[[1L]], paste0("n", seq.int(years) - 1))

  # Initialize population in timestep 0
  ## decide if all adults in earliest reproductive stage or else how distributed?
  ###FIXME: for now (and the RIdbt) stick them all in last/only reproductive stage
  #FIXME not doing  nmat[, "n0"] <- matrix(c(rep(0, stgfr - 1), n), ncol = 1)


  ###FIXME: for now determine stable stage distribution 
  ###XXX ASSUME population is at Stable Stage Distribution XXX ###
  ssdist <- popdemo::eigs(A, what = "ss")
  ##### number individuals in nmat[, n0] <-- rounding=slightly different than
  ###### n*ssdist
  #### Then use n divided among reproductive stages (according to ssdist)
  prop_dist_repro <- sum(ssdist[stgfr:stgn])
  prop_dist_nonrepro <- 1 - prop_dist_repro
  n_nonrepro <- (n * prop_dist_nonrepro) / prop_dist_repro
  nmat[stgfr:stgn, "n0"] <- round(n * (ssdist[stgfr:stgn] / prop_dist_repro), 0)
  #### Also use n to determine number of non-reproductive individuals in each stage
  nmat[which(seq.int(stgn) < stgfr), "n0"] <- round(n_nonrepro *
                                                ssdist[seq.int(stgn) < stgfr], 0)
  # Iterate through years to project individuals in each stage across years
  invisible(sapply(seq(from = 0, to = years - 2, by = 1),
    FUN = function(t){
      nmat[, paste0("n", t + 1)] <<- round(A %*% nmat[, paste0("n", t)], 0)
    }))  #<-- end sapply


  ## See how many Stage 1 individuals produced from each stage
  # Create matrix of stage contributions to stage 1 (whole individuals) 
  ## ROW is the stage contributing stage 1 individuals in a given column
  ## COLUMN is the year
  ## ENTRIES are the  number of stage 1 individuals produced by each
  ### previous year's stage (ROW) in a given year (COLUMN)
  ###TODO check above assumption that F matrix only has entries in first row
  #### (i.e., fertility only contributes to stage 1)
  ## break down matrix multiplication (inner product):
  ### matrix multiplication is trace of outer product
  #### (simplified in code since vectors always dimensions: 1 x n %o% n x 1)
  ### if do not do the summing step of the trace, then get each stage's contribution
  rpro_stg2stg1 <- cbind(matrix(NA, nrow = stgn, ncol = 1),
      sapply(seq(from = 0, to = years - 2, by = 1),
        FUN = function(t){
          round(diag(F[1, ] %o% nmat[, paste0("n", t)]), 0)
      }))  #<-- end function definition, sapply, and cbind
    dimnames(rpro_stg2stg1) <- list(dimnames(A)[[1L]],
                                    paste0("n", seq.int(years) - 1))


  # Total number of individuals in simulation
  ## All individuals in initial population (n0) plus offspring produced (top row)
  #FIXME if immigrants allowed, then will need to factor these in as well
  if(any(diag(F)[-1] != 0.0)){  #<-- check for immigrants through F
    stop("F matrix contains non-zeroes on diagonals")
  }
  N <- sum(nmat[, 1]) + sum(nmat[1, -1])

 return(list(nmat = nmat, N = N, stgfr = stgfr, rpro_stg2stg1 = rpro_stg2stg1))
}
# END FUNCTION DEFINITION    ############################ 
#########################################################

  nmat_N_female <- stg_sim_n(A, F, n_females)
  # TODO/FIXME: assume males same as females EXCEPT:
  ## male fecundity  = fecundity (in F) * ((1 - psr_p_female) / psr_p_female)
  ###                = fecundity (in F) scaled by sex ratio
  F_male <- F * ((1 - psr_p_female) / psr_p_female)
  U_male <- U  #<-- FIXME: now assumes same survival; TODO sex-specific survival
  nmat_N_male <- stg_sim_n(U_male + F_male, F_male, n_females) #<--FIXME assume n_female=n_male

  N <- nmat_N_female$N + nmat_N_male$N
  nmat_female <- nmat_N_female$nmat
  nmat_male <- nmat_N_male$nmat
  
  stgfr_f <- nmat_N_female$stgfr
  stgfr_m <- nmat_N_male$stgfr

  rpro_stg2stg1_f <- nmat_N_female$rpro_stg2stg1
  rpro_stg2stg1_m <- nmat_N_male$rpro_stg2stg1
#+++]	





#[+++
  pedigree <- data.frame(animal = seq.int(N),
    dam = NA,
    sire = NA,
    sex = NA,
    cohort = NA)    
  # initialize pedigree with founders
  ## (note founders can be from different cohorts)    
  ## reverse order: oldest show up at top and young-of-the year are last in line
  n0_alt_FM_by_stg <- c(t(cbind(rev(nmat_female[, 1]), rev(nmat_male[, 1]))))
  nchrt_ped <- sum(n0_alt_FM_by_stg)  #<-- number of individuals in this block of pedigree
  endchrt_ped <- nchrt_ped  #<-- initialize this before for loop
  pedigree$sex[1:nchrt_ped] <- unlist(mapply(rep, c("F", "M"), n0_alt_FM_by_stg,
                                                 USE.NAMES = FALSE))
  pedigree$cohort[1:nchrt_ped] <- unlist(mapply(rep,
  		rep(seq.int(stgn), each = 2) - stgn, n0_alt_FM_by_stg,
  		USE.NAMES = FALSE))
                                                  
#+++]	


#[|||  NOTE: comment below changed to: list stores who is ALIVE each year
## NOTE: XXX changing object name from `dat` to `census`
	# make list that stores who is alive in each year
## also change since we know how many years to produce	
#	dat <- list()
	census <- vector("list", length = years)
#|||]	
	# dat[[1]] <-  data.frame(animal = pedigree$animal, sex = pedigree$sex, age=NA, year=1)
	census[[1]] <-  data.frame(
		animal = pedigree$animal[1:sum(n0_alt_FM_by_stg)], 
		sex = pedigree$sex[1:sum(n0_alt_FM_by_stg)], 
#[|||
## NOTE: here stage indicates STAGE of all alive in population at this time 
#|||]
		stage = 1-pedigree$cohort[1:sum(n0_alt_FM_by_stg)], 
		year = 1)

	if(verbose) cat("starting pop created \n")


	## stores male-female pairings across years
#[|||	
##  since we know how many years to produce	
#	pairs <- list()
	pairs <- vector("list", length = years)
#|||]
	



      for(year in 1:(years - 1)){

	if(verbose) cat("year ",year,": ")

	####
	# BREEDING FEMALE AND MALE
	####

		# get vectors of females and males available to breed[---, accounting for the probability of females breeding---]
#		females <- subset(dat[[year]],sex=="F")$animal
#[|||   replaced age at first reproduction with stage at first reproduction
# also no p_breed (since part of matrix model) so females=breeding_females
## Can also get rid of "females" and "males" objects; not used for adult survival
		breeding_females <- subset(census[[year]], sex=="F" &
		                                  stage >= stgfr_f)$animal
#|||]		
#[---		## maybe just dont include these individuals in dat?
#---]
#[---
#		breeding_females <- females[as.logical(stats::rbinom(length(females),1,p_breed_f))]
#---]
		# dat[[year]][(match(females,dat[[year]]$animal)),"age"]
		
#[|||   replaced age at first reproduction with stage at first reproduction
# also no p_breed (since part of matrix model) so males=breeding_males
		breeding_males <- subset(census[[year]], sex=="M" &
		                                stage >= stgfr_m)$animal
#|||]		
#[---
#		breeding_males <- males[as.logical(stats::rbinom(length(males),1,p_breed_m))]
#---]
		if(verbose) cat("Breeding Individuals, ")
### maybe need something that says that if no breeding female skip offspring creation?

	####
	# PAIRING
	####

		#work out number of pairs that can be formed 
		n_pair <- length(breeding_females)#min(length(breeding_females),length(males))
		# print(length(males))
		# print(length(breeding_females))


		### assign 'social' male
		
#[||| removed condition if adult_surv==0 since in A, but hope works without
		if(p_retain==0 || year==1){ 
#|||]		
			## when there is no pairing to be done, just take random males, :
			social_male <- sample_males(breeding_females,breeding_males)
		}else{
			## MAKE VECTOR OF THOSE PAIRING UP AGAIN
			social_male<-sapply(breeding_females,function(bf){
				if(bf %in% pairs[[year-1]]$female
				# female bred in the last year 
				& pairs[[year-1]]$male[match(bf,pairs[[year-1]]$female)] %in% breeding_males
				# and paired male is breeding this year
				& stats::rbinom(1,1,p_retain)==1
				# and they retain each other
				 ){
				  pairs[[year-1]]$male[match(bf,pairs[[year-1]]$female)]
				}else{
					NA
				}
			})
			## Assign remaining males 
			social_male[is.na(social_male)] <- sample_males(
				females = breeding_females[is.na(social_male)], 
				males = breeding_males, 
				unpaired_males = breeding_males[!breeding_males%in%social_male]
				)
			# for the ones that havent been assigned, sample the males. 
			#there is a potential issue here if the number of females is larger than the males, then the males being chosen to sire the unpaired females will get lots of pairings, but the ones already paired wont get any extra.
		}
		## what do do with probability of breeding and mate retention?! If mate retention is 1, and one member of the pair doesnt breed once, then they will end up swapping, so there will be some level of divorce
		
		
# if female existed last year, who was male. else sample new male
		# if male alive, then rbinom(1,1,p_retain), 
		  # if retain then male, 
		  # else sample new male
		pairs[[year]] <- data.frame(
			female = breeding_females,
			male = social_male
			)		

		if(verbose) cat("Breeding pairs created, ")

	####
	# FEMALE FECUNDITY
	####
#[---    Need NEW approach to simulate fecundity	
#https://stat.ethz.ch/pipermail/r-help/2005-May/070680.html
# rztpois <- function(N, lambda) stats::qpois(stats::runif(N, stats::dpois(0, lambda), 1), lambda)

		# number of offspring per female
#		n_juv <- if(fixed_fecundity) {
#			rep(fecundity,n_pair)
#		}else{
#			stats::rpois(n_pair,fecundity)
#		}

#---]


#[+++   NEW APPROACH:
##	With a projection model, already have the number of offspring produced
##	Have the number of pairs and need to assign offspring numbers to pairs
##	Need to assign Poisson distributed number of offspring so each pair
###	  given a Poisson-distributed value, but such that the sum of these draws
###	  over all pairs equals the total number of offspring produced in the model
##	When there are several independent Poisson random variables, their JOINT
###	  behavior conditioned on their sum = MULTINOMIAL DISTRIBUTION
##	`rmultinom()` generates random counts from multinomial distribution, that
###	  when using fixed total and equal probabilities (or proportional to
###	  desired Poisson scale parameter) then function will draw counts
###	  conditioned on the total sum
##      The extra layer/complication is that each stage of reproductive females
###       are allowed to contribute different numbers of offspring on average.

  # number of offspring per female
  ##TODO this won't work/line-up if sex-specific stages of first reproduction
  total_juv_by_stg <- rpro_stg2stg1_f[stgfr_f:stgn, year + 1] +
                        rpro_stg2stg1_m[stgfr_m:stgn, year + 1]
  # number of offspring
  ## (`stg_sim_n()` uses F to create rpro_stg2stg1 so above should
  ### at least be a correct/pure *fertility* count)
  total_juv <- sum(total_juv_by_stg)

  # make general index for row of each breeding_female in census[[year]]
  brdng_f_ind <- match(breeding_females, census[[year]][, "animal"])
  # index for breeding females: which stage after stgfr does each belong to
  brdng_f_stgInd <- match(census[[year]][brdng_f_ind, "stage"], stgfr_f:stgn)
  # determine each breeding females stage-specific total juvenile production
  if(fixed_fecundity) { 
        n_juv_stg <- floor(total_juv_by_stg / table(brdng_f_stgInd))
        n_juv <- n_juv_stg[brdng_f_stgInd]
        # Check/adjust when number of per pair offspring was not an integer
        if((sum(n_juv) - total_juv) != 0){
          for(s in 1:length(total_juv_by_stg)){
            f_stg_s <- which(brdng_f_stgInd == s)
            diff_juv_s <- sum(n_juv[f_stg_s]) - total_juv_by_stg[s]
            if(diff_juv_s > 0) admin1 <- -1 else admin1 <- 1
            n_juv[f_stg_s[seq.int(abs(diff_juv_s))]] <-
                               n_juv[f_stg_s[seq.int(abs(diff_juv_s))]] + admin1
          }  #<-- end for s
        }  #<-- end if
        
    } else{    
        # create vector of probabilities for females
        ## different probability/assigned number offspring based on female's stage
        ### (but don't include stage not contributing any offspring)
        n_juv <- rep(0, length(breeding_females))  #<-- pre-allocate vector
        for(s in 1:length(total_juv_by_stg)){        
          if(total_juv_by_stg[s] == 0) next else{
            n_juv[which(brdng_f_stgInd == s)] <- as.vector(stats::rmultinom(n = 1,
                size = total_juv_by_stg[s],
                prob = rep(F[1, (stgfr_f + s - 1)], sum(brdng_f_stgInd == s))))
            }  #<-- end if/else 0 juv for stg s
        }  #<-- end for s  		
      }  #<-- end if/else fixed_fecundity
             
      
#+++]
#[---
#		total_juv <- sum(n_juv)
#---]
		if(verbose) cat("Eggs created, ")


	####
	# PATERNITY ASSIGNMENT
	####

### for each female
### start with 'social' male 
### how many eggs does he sire of total, with probability p_sire
### randomly choose a different male, with probability p_sire
### how many of the remaining eggs does he sire
### continue until all eggs are sired		

				# EPP
				### can make this more efficient - if p_polyandry=0 then rep(social_male,n_juv)
				#if p_polyandry==1 & p_sire==0 then sample(males,n_juv*length(breeding_females))

		sires <- if(p_polyandry==0) {
				rep(social_male,n_juv)
			}else if(p_polyandry==1 & p_sire==0){
				## should it use sample_male here, so that all males are sampled once, and then its random? because otherwise the p_breed doesn't work?
				sample(breeding_males,n_juv*length(breeding_females), replace=TRUE)
			}else{
				c(
					lapply(1:n_pair,function(i){
					## probability of any EPP
					polyandry <- stats::rbinom(1,1,p_polyandry)
					if(polyandry){
						## if there is EPP, how much
						## this is calculated by sampling how many of the offspring the paired male sired, and then giving the same probability to subsequent males. This means that extra pair males will be few, and have several offspring if p_sire if high - think this is more realistic
						if(p_sire==0) {
							## should each breeding male sire at least one? so use smaple_male?
							sample(breeding_males,n_juv[i])
						}else{
							within_sires <- fill_sires(n_juv[i],p_sire)
							if(length(within_sires)>0){
								c(social_male[i], sample(breeding_males,max(within_sires)-1,replace=TRUE))[within_sires]
								## have put replace=TRUE as if there are few males this might not work - some males might get chosen twice, but this will likely only happen when N is low, and so isn't unrealistic anyway
							}else{ NULL }
						}
						# n_sired <- rbinom(1,n_juv[i],p_sire)
						# c(rep(social_male[i], n_sired), sample(males,n_juv[i]-n_sired,replace=TRUE))
					}else{
						rep(social_male[i],n_juv[i])
					}
				})
					, recursive=TRUE)

			}

	if(verbose) cat("Paternity Assigned, ")


	####
	# PEDIGREE CREATION
	####
	# Set up starting and ending row indices for this cohort in the pedigree
        stchrt_ped <- endchrt_ped + 1 
        endchrt_ped <- stchrt_ped + total_juv - 1

		## make ped incorporating EPP and fecundity
#[|||		ped <- if(total_juv>0){  
#|||]
                pedigree[stchrt_ped:endchrt_ped, -1] <- if(total_juv > 0){

			data.frame(
#[---				#new animal ids 
#				animal = paste0(year,"_",1:total_juv),
#---]				
				# fecundity
				dam = rep(breeding_females,n_juv),
				sire = sires,
				#equal sex ratio
#[---				
#				sex = if(constant_pop){
				# in fixed pop, need to create exact amounts, otherwise it might not work!
					# bit convuluted, but makes sure there is the right number if there is an odd number of offspring
#					sample(rep(c("M","F"),ceiling(total_juv/2)),total_juv, replace=FALSE)
#				}else{
#					sample(c("M","F"),total_juv,replace=TRUE)
#				},
#---]
#[+++
				sex = sample(rep(c("M","F"), ceiling(total_juv/2)),
					total_juv, replace = FALSE),
#+++]
				cohort = year
			)
		}else{
			NULL
		}

		## create individuals present in the next year
#[--- matrix projection gives numbers of each stage class present in next year
## Use that to randomly sample IDs accordingly

	##----------
	## Juvenile Survival
	##----------
	
#		if(is.null(ped)){
#			next_year_juvF <- NULL
#			next_year_juvM <- NULL
#		}else{
#			pedM <- ped[ped[,"sex"]=="M",]
#			pedF <- ped[ped[,"sex"]=="F",]

#			next_year_juvF <- if(nrow(pedF)==0){
#				NULL	
#			}else if(constant_pop){
#				sample(pedF[,"animal"], round(juv_surv_f*fecundity*n_females/2), replace=FALSE)
#			}else{
#				pedF[as.logical(stats::rbinom(nrow(pedF),1,juv_surv_f)),"animal"]
#			}

			
#			next_year_juvM <- if(nrow(pedM)==0){
#				NULL
#			}else if(constant_pop){
				### need to ensure equal sex ratio of recruits, otherwise population size fluctuations
#				sample(pedM[,"animal"], round(juv_surv_m*fecundity*n_females/2), replace=FALSE)

#			}else{
#				pedM[as.logical(stats::rbinom(nrow(pedM),1,juv_surv_m)),"animal"]

#			}
#		}
		
#		if(verbose) cat("Juvenile survival, ")


		## or could take the new recruits from the whole pedigree, subset by cohort = year-afr - maybe this is better, because otherwise census always has a load of pre afr individuals in it

	##----------
	## Adult Survival
	##----------

#		next_year_AF <- if(adult_surv_f==0){
#			NULL
#		}else if(constant_pop){
#			sample(females, round(adult_surv_f*n_females), replace=FALSE)
#		}else{
#			females[as.logical(stats::rbinom(length(females),1,adult_surv_f))]
#		}
#		
#		next_year_AM <- if(adult_surv_m==0){
#	    NULL
#		}else if(constant_pop){
#			sample(males, round(adult_surv_m*n_females), replace=FALSE)
#		}else{
#			males[as.logical(stats::rbinom(length(males),1,adult_surv_m))]
#		}
#		
#---]


	if(verbose) cat("Adult survival, ")
#[+++ Determine transition of individual IDs from one stage into next according
## to the projection matrices for survival
### Already created number of newly created individuals for next year

 #FIXME: Handle this better in future
 if(year == 1){
   cat("\nsurvival assumes PRE-breeding census\n",
     "(excludes first year survival: this is incorporated into fertility)\n")
  }  
  # create storage matrix for where next year individuals come from
  ## ROW = the stage receiving n individuals
  ## COLUMN = the stage from which the n individuals currently reside
  surv_stg2stgFun <- function(s, U, nmat){
    round(diag(U[s, ] %o% nmat[, paste0("n", year - 1)]), 0)
  }
  surv_stg2stg_yr_female <- t(sapply(seq.int(stgn),
    FUN = surv_stg2stgFun,
    U = U, nmat = nmat_female))
    dimnames(surv_stg2stg_yr_female)[[1L]] <- dimnames(surv_stg2stg_yr_female)[[2L]]
  surv_stg2stg_yr_male <- t(sapply(seq.int(stgn),
    FUN = surv_stg2stgFun,
    U = U_male, nmat = nmat_male))
    dimnames(surv_stg2stg_yr_male)[[1L]] <- dimnames(surv_stg2stg_yr_male)[[2L]]

  # below gives index in census[[year]] for individuals "selected" to survive from
  ## each stage as well as their stage in the next year
  next_year_AF <- do.call(rbind, sapply(seq.int(stgn),
    FUN = function(s){
      ssub <- surv_stg2stg_yr_female[, s]
      if(all(ssub == 0)) return(NULL) else{
        yrIndx <- sample(with(census[[year]], which(sex == "F" & stage == s)),
                         size = sum(ssub),
                         replace = FALSE)
        nxtYrStg <- unlist(mapply(rep, which(ssub != 0), ssub[which(ssub != 0)]))
        return(structure(cbind(yrIndx, nxtYrStg),
                         dimnames = list(NULL, c("yrIndx", "nxtYrStg"))))
        }  #<-- end else
    })  #<-- end sapply
  )  #<-- end do.call  
  next_year_AF[] <- next_year_AF[order(next_year_AF[, "yrIndx"]), ]
  
  next_year_AM <- do.call(rbind, sapply(seq.int(stgn),
    FUN = function(s){
      ssub <- surv_stg2stg_yr_male[, s]
      if(all(ssub == 0)) return(NULL) else{
        yrIndx <- sample(with(census[[year]], which(sex == "M" & stage == s)),
                         size = sum(ssub),
                         replace = FALSE)
        nxtYrStg <- unlist(mapply(rep, which(ssub != 0), ssub[which(ssub != 0)]))
        return(structure(cbind(yrIndx, nxtYrStg),
                         dimnames = list(NULL, c("yrIndx", "nxtYrStg"))))
        }  #<-- end else
    })  #<-- end sapply
  )  #<-- end do.call    
  next_year_AM[] <- next_year_AM[order(next_year_AM[, "yrIndx"]), ]


#+++]


#[--- turning OFF immigration for now: In future, it will probably be represented in a transition matrix anyway (e.g., diagonal element of F)

	##----------
	## IMMIGRATION
	##----------

		## need to sort out with constant pop
		## might be worth giving the immigrants age=afr, so cohort=year-afr

#		n_imm <- if(constant_pop){
#			c(immigration_f,immigration_m)*n_females
#		}else{
#			stats::rbinom(2,n_females,c(immigration_f,immigration_m))
#		}
		
#		imm_females <- if(n_imm[1]>0){
#			data.frame(
#				animal=paste(year+1,"IF",seq_len(n_imm[1]),sep="_"),
#				dam=NA,
#				sire=NA,
#				sex="F",#rep(c("F","M"),c(n_imm)),
#				cohort=year-afr +1
#				)
#			# paste(year+1,"IF",seq_len(n_imm[1]),sep="_")	
#		}else{NULL}
		
#		imm_males <- if(n_imm[2]>0){
#			data.frame(
#				animal=paste(year+1,"IM",seq_len(n_imm[2]),sep="_"),
#				dam=NA,
#				sire=NA,
#				sex="M",#rep(c("F","M"),c(n_imm)),
#				cohort=year-afr +1
#				)
#			# paste(year+1,"IM",seq_len(n_imm[2]),sep="_")
#		}else{NULL}


#		immigrants <- if(year!=years){
#			rbind(imm_females,imm_males)
#			# data.frame(
#			# 	animal=c(imm_females,imm_males),
#			# 	dam=NA,
#			# 	sire=NA,
#			# 	sex=rep(c("F","M"),c(n_imm)),
#			# 	cohort=year-afr +1
#			# 	)
#		}else{
#			NULL
#		}
#
#	if(verbose) cat("immigration, ")
#---]


	##----------
	## OUTPUTS
	##----------
		
#[--- also removing immigrants
		#update pedigree
#		pedigree <- rbind(pedigree,ped)#,immigrants)
#---]
#[--- record this in dat instead
		# put together all individuals surviving to the next time step
#		next_year_ind <- data.frame(
#			animal = c(next_year_juvF,
#			 	next_year_juvM,
#			 	next_year_AF,
#			 	next_year_AM,
#			 	immigrants[,"animal"]),
#			sex= c(rep(c("F","M","F","M"),c(length(next_year_juvF),length(next_year_juvM),length(next_year_AF),length(next_year_AM))),immigrants[,"sex"])
#		)

#	  if(nrow(next_year_ind)==0){
#	  	message("Population went extinct in year ", year+1)
#	  	break
#	  }

		## maybe we don't need this - redundant with cohort?
#	  next_year_ind$age <- (year+1) - pedigree[match(next_year_ind$animal,pedigree$animal),"cohort"]
#	  next_year_ind$year <- year+1
#---]
#[|||	  dat[[year+1]] <- next_year_ind
          census[[year + 1]] <- rbind(data.frame(
                       animal = census[[year]][c(next_year_AF[, "yrIndx"],
                                           next_year_AM[, "yrIndx"]), "animal"],
                       sex = c(rep("F", nrow(next_year_AF)),
                               rep("M", nrow(next_year_AM))),
                       stage = c(next_year_AF[, "nxtYrStg"],
                                 next_year_AM[, "nxtYrStg"]),
                       year = rep(I(year + 1),
                                    I(nrow(next_year_AF) + nrow(next_year_AM)))),
                  data.frame(animal = pedigree$animal[stchrt_ped:endchrt_ped],
                       sex = pedigree$sex[stchrt_ped:endchrt_ped],
                       stage = 1,
                       year = year + 1))
	  
		
#	  if(length(unique(next_year_ind$sex))==1){
          if(length(unique(census[[year + 1]][, "sex"])) == 1){
#	  	if(unique(next_year_ind$sex)=="F"){
                if(unique(census[[year + 1]][, "sex"]) == "F"){
	  		message("Males went extinct in year ", year+1)	
	  	}
#	  	if(unique(next_year_ind$sex)=="M"){
                if(unique(census[[year + 1]][, "sex"]) == "M"){
	  		message("Females went extinct in year ", year+1)	
	  	}
	  	break
	  } 
#	  if(length(unique(next_year_ind$sex))==0){
          if(length(unique(census[[year + 1]][, "sex"])) == 0){
	  	message("Population went extinct in year ", year+1)
	  	break
	  }
#|||]

	if(verbose) cat("\n")

      }  #<-- end for year

#[+++ # intermediate step to assign ages in pedigree
    data_str <- do.call(rbind, census)
    ## find individuals in data_str with match
    ### take this year and subtract from it the stage of first appearance to get age
    #### XXX Assumes founding cohort stages all == 1 year in length
    ### match returns first match, so go bottom-to-top of data_str so match
    #### returns index of oldest occurrence of individual in data_str
    pedigree$minYr <- data_str$year[match(pedigree$animal, data_str$animal)] 
    pedigree$mxYr <- rev(data_str$year)[match(pedigree$animal,
                                               rev(data_str$animal))] 
    pedigree$lstObsStg <- ifelse(pedigree$mxYr == 10,
                                NA, 
                                rev(data_str$stage)[match(pedigree$animal,
                                                         rev(data_str$animal))])
    pedigree$longevity <- with(pedigree, ifelse(mxYr == 10, NA, mxYr - cohort))
    
       
#+++]    
#[|||
 return(list(pedigree=pedigree, data_str = data_str))
#|||]
}

