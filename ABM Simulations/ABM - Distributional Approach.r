###Libraries
library(plyr)
library(data.table)
library(ggplot2)
library(mixtools)
library(gridExtra)

###Shorthand
PARAMS = 1 ##PARAMETER LIST
TOKENS = 2 ##TOKEN SET

##PARAM ROW [[category, mean_values, cov, count]]
CAT = 1 ##Category Number
MEAN = 2 ##Mean Vector
SIG = 3 ##Covariance Matrix
COUNT = 4 ##Prior Probability

###Helper Functions
##Returns a list consisting of the nth element of each list in x
listslice <- function(x, n){
	return(lapply(x, '[[', n))
}

##ParamEM - Finds parameters for tokens
ParamEM <- function(df, tol=1e-10,cats=0){
    if (cats > 0){ ##Force function to report a specific number of categories 
        params <- llply(1:cats,function(x){
            ss <- df[df[[1]]==x,]
            if (nrow(ss)<1){
                means = rep(0,ncol(df)-1)
                d = length(means)
                names(means) <- paste0("value",1:d)
                return(list(1, means, diag(d),0))
            } else {
                return(list(1, unlist(colMeans(ss[,-c(1),drop=F])), ##Mean
                cov(ss[,-c(1),drop=F]), #Sigma
                nrow(ss)/nrow(df))) #Lambda
            }
        })
    } else { ##Otherwise estimate as many categories as appear in the data
        ##Initialize Values with Category Labels
        params <- dlply(df, c(colnames(df)[1]), function(x){
            return(list(1, unlist(colMeans(x[,-c(1),drop=F])), ##Mean
                cov(x[,-c(1),drop=F]), #Sigma
                nrow(x)/nrow(df))) #Lambda
        })    
    }
    
    ##Fix Degenerate Categories
    for(a in 1:length(params)){
        if (all(is.na(params[[a]][[3]])) || any(eigen(params[[a]][[3]])$values < tol)){ ##If NAs or non-positive definite matrix
            params[[a]][[3]] = diag(dim(params[[a]][[3]])[1]) #Make into Identity Matrix
            params[[a]][[2]] = rep(0,length(params[[a]][[2]])) #Set mean to 0 vector
            names(params[[a]][[2]]) <- paste0("value",1:length(params[[a]][[2]]))
            temp_lam = params[[a]][[4]]
            params[[a]][[4]] = 0 ##Adjust param to 0
            for(b in 1:length(params)){
                if(params[[b]][[4]]){
                    params[[b]][[4]] = params[[b]][[4]]/(1-temp_lam) ##Readjust prior probabilities accordingly
                }
            }
            
        }
    }   
    return(llply(1:length(params), function(x){
		return(list(x, params[[x]][[MEAN]], params[[x]][[SIG]], params[[x]][[COUNT]]))
	}))
}

##produceToken - generate n tokens from a gaussian mixture
produceToken <- function(params, n=1){
	tokens <- ldply(1:max(1,floor(n)), function(x){
		category <- sample(1:length(params), 1, prob=unlist(listslice(params, COUNT)))
		catinfo <- params[[category]]
		return(c(category, rmvnorm(1, catinfo[[MEAN]], catinfo[[SIG]])))
	})
	colnames(tokens) <- c("cat", paste0("value",1:length(listslice(params, MEAN)[[1]])))
    if(n < 1){
        return(tokens[F,])
    }
	return(tokens)
}

##perceiveToken - return conditional probability of each category given the token
perceiveToken <- function(token, params){
	priorc <- unlist(listslice(params,COUNT))
	xgivenc <- unlist(ldply(params, function(x) dmvnorm(unlist(token[-c(CAT),drop=F]), x[[MEAN]], x[[SIG]])))
	priorx <- sum(xgivenc*priorc)
	return(priorc*xgivenc/priorx)
}

##createAgents - initialize nagents agents with params parameters
createAgents <- function(nagents, params){
	results <- llply(1:nagents, function(i){
        agentTokens <- produceToken(params, 0)
        return(list(params, agentTokens))
	})
	return(list(lapply(results, '[[', PARAMS), lapply(results, '[[', TOKENS)))
}

###runSimulation - Runs Simulation
##steps: Number of iterations the simulation will run for
##globalparams: The initial distribution for the phonetic categories
##nagents: The number of learner agents per generation
##ntokens: The number of tokens a learner will categorize per generation
##misperc: In the event that a token is miscategorized, the likelihood this token will be stored
##speaker_cat: In the event that a token is miscategorized, the likeilhood it is stored as the category intended by the speaker
##cat_map: A list defining a mapping between categories (e.g., from vowel dependent to vowel independent)
runSimulation <- function (steps, globalparams, nagents, ntokens, misperc, speaker_cat, cat_map=NULL, genparams=NULL){
    cats=length(globalparams)
    feats = length(globalparams[[1]][[2]])
    
    ##Initalize Parameters
    temp_params = globalparams
    if(!is.null(cat_map) && !is.null(genparams)){
        temp_gen = genparams
    } else {
        temp_gen = NULL
    }
    
    speakers <- createAgents(nagents, temp_params)
    
    ##Save snapshot of current generation
    piece <- ldply(speakers[[1]], function(x){
        test <- cbind(1:length(x),ldply(lapply(listslice(x,3),det)),ldply(listslice(x,4)),ldply(listslice(x,2)))
        colnames(test)[1:3] <- c("cat","dets","prior")
        return(test)
    })
    piece$step <- 0
    snapshot <- piece
    snap_end = 3+feats
    colnames(snapshot)[1:snap_end] <- c("cat","dets","prior",paste0('value',1:feats))
    
    for(y in 1:steps){
        temp_tokens = produceToken(temp_params, 0)
        ##Initiate new learners
        listeners = createAgents(nagents, temp_params)
        ##Pair teachers with learners        
        ##Teacher token production

        for(x in 1:nagents){
            ##Teacher token production
            token <- produceToken(speakers[[PARAMS]][[x]],ntokens)
            
            ##Learner categorization
            if(!is.null(temp_gen)){
                outcome_gen <- as.data.frame(t(apply(token, 1, function(z) perceiveToken(z, temp_gen))))
                listener_categories_gen <- apply(outcome_gen,1,which.max)
            }
            outcome <- as.data.frame(t(apply(token, 1, function(z) perceiveToken(z, listeners[[PARAMS]][[x]]))))
            listener_categories <- apply(outcome,1,which.max)
                
            if(!is.null(temp_gen)){
                correct <- mapvalues(token[[1]], cat_map[[1]], cat_map[[2]], warn_missing = F) == listener_categories_gen ##Correctly Identified
            } else {
                correct <- token[[1]]==listener_categories ##Correctly Identified
            }
                
            listener_tokens <- llply(1:length(correct), function(y){
                if(correct[y]){
                    return(c(token[[1]][y], unlist(token[y,-c(CAT)])))
                } else if (runif(1) <= misperc){ ##If token still classified
                    if(runif(1) <= speaker_cat){ ##Classify it is the category intended by the speaker
                        return(c(token[[1]][y], unlist(token[y,-c(CAT)])))
                    } else {
                        return(c(listener_categories[y], unlist(token[y,-c(CAT)])))
                    }
                }
                return()
            })
            listener_tokens <- ldply(listener_tokens[!is.null(listener_tokens)])
            colnames(listener_tokens)[c(1,2)] <- c("cat",'value1')

            ##Re-estimation of listener parameters
            listeners[[PARAMS]][[x]] <- ParamEM(listener_tokens,cats=cats)
            temp_tokens <- rbind(temp_tokens, listener_tokens)
        }
        ##Update global distribution
        temp_params <- ParamEM(temp_tokens,cats=cats)
        if(!is.null(temp_gen)){
            temp_gen <- ParamEM(cbind(cat=mapvalues(temp_tokens$cat,cat_map[[1]], cat_map[[2]]),temp_tokens[,-c(1)],warn_missing=F),cats=cats)
        }
        
        ##Learners become teachers
        speakers <- listeners
        
        ##Save snapshot of current generation
        piece <- ldply(listeners[[1]], function(x){
            test <- cbind(1:length(x),ldply(lapply(listslice(x,3),det)),ldply(listslice(x,4)),ldply(listslice(x,2)))
            colnames(test)[1:3] <- c("cat","dets","prior")
            return(test)
        })
        piece$step <- y
        snapshot<-rbind(snapshot,piece)
        print(paste0("Finished step ",y,"."))
    }
    return(snapshot)
}
