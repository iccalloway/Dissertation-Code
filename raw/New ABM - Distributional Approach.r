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

##ParamEM - Finds parameters for Gaussian Mixture
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
probs <- adply(kt_f[,c(2,7,8)],1, function(x) perceiveToken(x,globalparams))


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






##SIMULATIONS FOR A SINGLE CONTEXT##

##Subset by vowel context
globaltokens <- kp_f[height=="low" & backness=="back",]
#globaltokens <- fth_f
features <- c('V12', 'V21')
features <- colnames(globaltokens)[(ncol(globaltokens)-1):ncol(globaltokens)]
#globalparams <- ParamEM(globaltokens[,c('cat',features)])


##Category related info
conditions <- globaltokens$segment
conditions_un <- unique(conditions)
globaltokens$cat <- mapvalues(conditions,conditions_un, 1:length(conditions_un))    
category_key <- unique(globaltokens[,c('segment','height','backness','cat')])

##Parameter Estimation
globalparams <- ParamEM(globaltokens[,c('segment',features),with=F])


#sub <- results[results$step %in% c(0,10) & results$cat==1,]
    


theme_set(theme_grey(base_size=21))
results<-runSimulation(40, globalparams, 20, 100,1,0)


colnames(results)[4:5] <-c("F12", "F21")


##Plotting
results$catname<-mapvalues(results$cat, 1:length(conditions_un),rev(conditions_un))
#colnames(results)[4:(4+length(features)-1)] <- features
results_m <- melt(results, measure.vars=c("F12", "F21"))

##Prior Evolution
ggplot(results_m,aes(x=step, y=prior,colour=as.factor(catname))) +
geom_smooth() +
ylim(c(0,1)) +
labs(x="Generation", y="Prior Probability", colour="Consonant") +
scale_color_manual(values=c("k"="#D55E00", "t"="#0072B2", "p" = "#009E73", "th" = "#E69F00", "f"="#56B4E9"))


t.test(prior ~ step, data=diffs[segment=="k" & step %in% c(0,40),])



##Featural Evolution
ggplot(results_m[results_m$prior > 0,],aes(x=step, y=value,colour=as.factor(catname))) +
geom_smooth() +
facet_grid(~variable) +
labs( x="Generation", y="Feature Value", colour="Consonant")+
scale_color_manual(values=c("k"="#D55E00", "t"="#0072B2", "p" = "#009E73", "th" = "#E69F00", "f"="#56B4E9"))


##Divergence
step=results[(1:(nrow(results)/2))*2,c('step')]
diffs<-cbind(results[(1:(nrow(results)/2))*2,c('F12','F21')]-results[(1:(nrow(results)/2))*2-1,c('F12','F21')],step)
ggplot(melt(diffs, measure.vars=c('F12','F21')), aes(x=step, y=value, colour=variable)) + geom_smooth()

t.test(F12 ~ step, data=diffs[step %in% c(0,40),])
t.test(F21 ~ step, data=diffs[step %in% c(0,40),])

##Determinant Evolution
ggplot(results_m[results_m$prior > 0,],aes(x=step, y=log10(dets),colour=as.factor(catname))) +
geom_smooth() +
labs(x="Generation", y="Log of Generalized Variance", colour="Consonant") +
scale_color_manual(values=c("k"="#D55E00", "t"="#0072B2", "p" = "#009E73", "th" = "#E69F00", "f"="#56B4E9"))

grid.arrange(p1,p2,p3,nrow=2,layout_matrix=rbind(c(1,1),c(2,3))) 

##SIMULATIONS WITH CONSONANT GENERAL PERCEPTION##

##Subset to include only point vowel contexts
contexts = data.table(height = c('high','high','low'), backness=c('front','back','back'))
globaltokens <- fth_f[contexts,on=c(height='height',backness='backness')]

##Vowel Context Specific Category Info
conditions <- interaction(globaltokens$segment,globaltokens$height,globaltokens$backness)
conditions_un <- unique(conditions)
globaltokens$cat <- as.numeric(droplevels(mapvalues(conditions,conditions_un, 1:length(conditions_un))))

##Consonant General Category Info
gen_conditions_un <- unique(globaltokens$segment)
globaltokens$gencat <- mapvalues(globaltokens$segment,gen_conditions_un, 1:length(gen_conditions_un))

##Mapping between vowel-specific and vowel-general categories
cat_un <- unique(globaltokens[,c('cat','gencat')])
cat_map <- list(as.numeric(as.character(cat_un$cat)), as.numeric(cat_un$gencat))
cat_info <- unique(globaltokens[,c('segment','height','backness','cat')])

##Initial Parameter Estimation
features <- c('PC1', 'PC2')

globalparams <- ParamEM(globaltokens[,c("cat","PC1", "PC2")])
results<-runSimulation(20, globalparams, 30, 400,1,1,cat_map)
colsend <- 3+length(features)
colnames(results)[4:colsend] <- features
results_merge <-merge(results,cat_info, by.x='cat',by.y='cat')
results_m <- melt(results_merge, measure.vars=features)

subs <- results_m[results_m$height=="high" & results_m$backness=="back",]

subs <-results_m
##Plotting
ggplot(subs[subs$prior > 0,],aes(x=step, y=value,colour=as.factor(interaction(segment,height,backness)))) +
geom_smooth() +
facet_grid(~variable) +
labs(title="Evolution of feature values over time", x="Simulation Iteration", y="Feature Value", colour="Category")
scale_color_manual(values=c("k"="#D55E00", "t"="#0072B2", "p" = "#009E73", "th" = "#E69F00", "f"="#56B4E9"))

ggplot(subs,aes(x=step, y=prior,colour=interaction(segment,height,backness))) +
geom_smooth() +
ylim(c(0,1)) +
labs(title="Evolution of prior probabilties over time", x="Simulation Iteration", y="Prior Probability", colour="Category")
scale_color_manual(values=c("k"="#D55E00", "t"="#0072B2", "p" = "#009E73", "th" = "#E69F00", "f"="#56B4E9"))

ggplot(subs[subs$prior > 0,],aes(x=step, y=dets,colour=as.factor(segment))) +
geom_smooth() +
labs(title="Evolution of generalized variance over time", x="Simulation Iteration", y="Generalized Variance", colour="Category")
grid.arrange(p1,p2,p3,nrow=2,layout_matrix=rbind(c(1,1),c(2,3))) 


