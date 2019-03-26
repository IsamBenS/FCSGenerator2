generate.pattern.codes <- function(nmb.pattern.codes, nmb.params, nmb.seuils)
{
    pattern.codes <- list()
    for (i in 1:nmb.pattern.codes) 
    {
        tmp.code <- c()
        while(length(tmp.code)==0 || list(tmp.code) %in% pattern.codes)
        {
            tmp.code <- sapply(1:nmb.params, function(i)
            {
                return(sample(1:nmb.seuils, 1))
            })
        }
        pattern.codes[[i]] <- tmp.code
    }
    
    return(pattern.codes)
}


generate.position.from.pattern <- function(pattern.code, lower.value, upper.value, nmb.pos=2)
{
    l <- length(pattern.code)
    range.unit <- (upper.value - lower.value)/nmb.pos
    m.list <- sapply(1:l, function(i)
    {
        pop.mean <- pattern.code[i]
        pop.mean <- runif(1, (pop.mean-2/3)*range.unit, (pop.mean-1/2)*range.unit)
        return(pop.mean)
    })
    sd.list <- sapply(1:l, function(i)
    {
        return(range.unit/l/3)
    })
    
    return(list(m.list,sd.list))
}


generate.pop.from.position <- function(pop.m.list, pop.sd.list, nmb.events,
                                      lower.value, upper.value, limited=T)
{
    new.pop <- matrix(nrow = nmb.events, ncol = length(pop.m.list))
    l <- length(pop.m.list)
    range.unit <- (upper.value - lower.value)/l
    
    for(i in 1:l)
    {
        new.pop[,i] <- rtruncnorm(nmb.events, lower.value, upper.value, pop.m.list[[i]], pop.sd.list[[i]])
    }
    return(new.pop)
}


extract.position.from.pop <- function(pop, unused.columns = NULL) #POP : matrice, row = events, col = params
{
    
    col.list <- 1:ncol(pop)
    if(!is.null(unused.columns))
    {
        col.list <- col.list[-unused.columns]
    }
    m.list <- sapply(1:length(col.list), function(i)
    {
        return(trunc(mean(pop[,col.list[[i]]])*1000)/1000)
    })
    sd.list <- sapply(1:length(col.list), function(i)
    {
        val <- trunc(sd(pop[,col.list[[i]]])*1000)/1000
        if(is.na(val))
        {
            val <- 0
        }
        return(val)
    })
    
    return(list(m.list, sd.list))
}


move.population <- function(pop, new.m.list, new.sd.list, m.ranges = list(0.1), sd.ranges = list(0.1), limited = T, 
                            max.val = 4.5, min.val = -0.5, unused.columns = NULL)
{
    col.list <- 1:ncol(pop)
    if(!is.null(unused.columns))
    {
        col.list <- col.list[-unused.columns]
    }
    for(i in 1:length(col.list))
    {
        me <- 0
        sd <- 0
        id <- col.list[[i]]
        ##NEW M--------------------------------------------------------------------------------
        if(length(m.ranges) == length(new.m.list))
        {
            me <- new.m.list[[i]] + runif(1, -m.ranges[[i]], m.ranges[[i]])
        }
        else
        {
            me <- new.m.list[[i]] + runif(1, -m.ranges[[1]], m.ranges[[1]])
        }
        
        ##NEW SD-------------------------------------------------------------------------------
        if(length(sd.ranges) == length(new.sd.list))
        {
            sd <- abs(new.sd.list[[i]] + runif(1, -sd.ranges[[i]], sd.ranges[[i]]))
        }
        else
        {
            sd <- abs(new.sd.list[[i]] + runif(1, -sd.ranges[[1]], sd.ranges[[1]]))
        }
        if(nrow(pop)==1)
        {
            sd <- 0
        }
        
        ##RANDOMIZE POP------------------------------------------------------------------------
        pop[,id] <- rtruncnorm(nrow(pop), min.val, max.val, me, sd)
    }
    
    return(pop)
}


reduce.population <- function(events, red.size) #red.size = % to remove
{
    events.to.remove <- sample(events, as.integer(red.size/100*length(events)))
    
    return(events.to.remove)
}


create.pop.from.pop <- function(pop, inc.size, unused.columns = NULL, limited = T, max.val = 4.5, min.val = -0.5) #inc.size = % to add
{
    col.list <- 1:ncol(pop)
    if(!is.null(unused.columns))
    {
        col.list <- col.list[-unused.columns]
    }
    
    new.pop <- matrix(nrow = as.integer(inc.size/100*nrow(pop)), ncol = ncol(pop))
    pop.pos <- extract.position.from.pop(pop, unused.columns = unused.columns)
    for(i in 1:length(col.list))
    {
        id <- col.list[[i]]
        new.pop[,id] <- rtruncnorm(nrow(new.pop), min.val, max.val, pop.pos[[1]][[i]], pop.pos[[2]][[i]])
    }
    for(i in unused.columns)
    {
        new.pop[,i] <- pop[1,i]
    }
    
    return(new.pop)
}


list.to.matrix <- function(list.to.convert)
{
    m <- matrix(nrow=length(list.to.convert), ncol=length(list.to.convert[[1]]))
    for(i in 1:nrow(m))
    {
        for(j in 1:ncol(m))
        {
            m[i,j] <- list.to.convert[[i]][[j]]
        }
    }
    return(m)
}


generate.pattern.codes.from.positions <- function(m.lists, 
                                             param.thresholds = 0.05) #Thresholds : error %
{
    pattern.list <- list()
    if(length(m.lists)>0)
    {
        m.matrix <- list.to.matrix(m.lists)
        nmb.pop <- nrow(m.matrix)
        nmb.params <- ncol(m.matrix)
        tmp.pattern.list <- matrix(nrow = nmb.pop, ncol = nmb.params)
        for(i in 1:nmb.params) #parcourt les parametres
        {
            values.order <- order(as.numeric(m.matrix[,i]))
            sorted.values <- matrix(ncol=2, nrow=nmb.pop)
            sorted.values[,1] <- sort(as.numeric(m.matrix[,i]))
            for(j in 1:nmb.pop)
            {
                val <- 1
                if(j>1)
                {
                    val <- sorted.values[j-1,2]
                    if(abs(sorted.values[j,1]/sorted.values[j-1,1]-1) > param.thresholds)
                    {
                        val <- val+1
                    }
                }
                sorted.values[j,2] <- val
            }
            tmp.pattern.list[as.numeric(unlist(values.order)),i] <- sorted.values[,2]
        }
        
        pattern.list <- lapply(1:nrow(tmp.pattern.list), function(i)
        {
            return(as.vector(as.numeric(tmp.pattern.list[i,])))
        })
    }
    
    return(pattern.list)
}

generate.fcs <- function(nmb.events, nmb.dim, nmb.populations, nmb.rare.populations, rare.pop.freq, 
                         min.freq.list, max.freq.list, min.val=-0.5, max.val=4.5)
{
    fcs <- NULL
    if(rare.pop.freq<100)
    {
        if(sum(min.freq.list)<=100 && sum(max.freq.list)>=100)
        {
            pop.frequencies <- min.freq.list
            freq.rep.order <- rev(order(max.freq.list-pop.frequencies))
            nmb.free.pop <- nmb.populations
            #==POP FREQUENCIES GENERATED HERE
            while((100-sum(pop.frequencies))>1/1000000)
            {
                remaining.percentage <- trunc(1000000*(100-sum(pop.frequencies)))/1000000
                diff.freq <- max.freq.list-pop.frequencies
                
                for(i in 1:nmb.populations)
                {
                    if(pop.frequencies[freq.rep.order[i]] < max.freq.list[[freq.rep.order[i]]] &&
                       remaining.percentage>0)
                    {
                        attributed.percentage <- min(remaining.percentage/nmb.free.pop, 
                                                     diff.freq[freq.rep.order[i]])
                        
                        pop.frequencies[freq.rep.order[i]] <- pop.frequencies[freq.rep.order[i]] + 
                            attributed.percentage
                        
                        remaining.percentage <- remaining.percentage - attributed.percentage
                        diff.freq <- max.freq.list-pop.frequencies
                    }
                }
            }
            #==
            pop.patterns.list <- generate.pattern.codes(nmb.populations, nmb.dim, 2)
            pop.position.list <- lapply(1:nmb.populations, function(i)
            {
                return(generate.position.from.pattern(pop.patterns.list[[i]], min.val, max.val, 2))
            })
            
            file.mat <- NULL
            event.id <- 0
            for(current.pop in 1:nmb.populations)
            {
                pop.events <- max(1,as.integer(pop.frequencies[current.pop]*nmb.events/100))
                event.id <- event.id+pop.events
                if(event.id<nmb.events && current.pop==nmb.populations)
                {
                    pop.events <- pop.events+(nmb.events-event.id)
                }
                pop.mat <- generate.pop.from.position(pop.position.list[[current.pop]][[1]],
                                                      pop.position.list[[current.pop]][[2]],
                                                      pop.events,lower.value = min.val, upper.value = max.val)
                pop.mat <- cbind(pop.mat, 0)
                pop.mat <- cbind(pop.mat, current.pop)
                if(is.null(file.mat))
                {
                    file.mat <- pop.mat
                }
                else
                {
                    file.mat <- rbind(file.mat, pop.mat)
                }
            }
            
            colnames(file.mat) <- c(rep("NA",times = nmb.dim),"time","Population")
            for (current_dim in 1:nmb.dim)
            {
                colnames(file.mat)[current_dim] <- paste("PARAM_",current_dim, sep="")
            }
            
            fcs <- flowFrame(file.mat)
            descR <- description(fcs)
            lapply(c(1:dim(file.mat)[2]),function(x)
            {
                descR[[paste0("$P",x,"R")]] <<- 262144
            })
            
            nmb.grp <- min(nmb.events,1000)
            descR[["TIMESTEP"]] <- 1/nmb.grp
            for(e in 1:nmb.events)
            {
                file.mat[e,nmb.dim+1] <- descR[["TIMESTEP"]] *  as.integer(e / (nmb.events/nmb.grp))
                
            }
            fcs <- flowFrame(file.mat, description = descR)
        }
    }
    
    return(fcs)
}

create.pop.table.from.fcs.matrix <- function(fcs.mat, populations.list, used.markers, cluster.column, pop.max=200)
{
    tmp.df <- NULL
    tmp.mat <- NULL
    if(!is.null(fcs.mat) && length(populations.list) > 0 && length(populations.list) < pop.max && !is.null(cluster.column))
    {
        pop.ids <- populations.list
        col.id <- cluster.column
        
        tmp.mat <- matrix(ncol=length(used.markers), nrow=length(pop.ids))
        colnames(tmp.mat) <- c(colnames(fcs.mat)[used.markers])
        lapply(1:length(pop.ids), function(j)
        {
            if(length(pop.ids[[j]])>0)
            {
                curr.pop.id <- unique(fcs.mat[pop.ids[[j]],col.id])
                curr.pop <- fcs.mat[fcs.mat[,col.id]==curr.pop.id, used.markers]
                if(length(pop.ids[[j]])==1)
                {
                    curr.pop <- t(fcs.mat[fcs.mat[,col.id]==curr.pop.id, used.markers])
                }
                tmp.val <- extract.position.from.pop(curr.pop)
                tmp.mat[j,] <<- as.numeric(unlist(unlist(tmp.val[[1]])))
            }
        })
        rownames(tmp.mat) <- names(populations.list)
        
        tmp.df <- expand.grid(Populations=names(pop.ids), Markers=colnames(fcs.mat)[used.markers])
        tmp.df$value <- as.numeric(unlist(as.list(tmp.mat)))
    }
    return(list(tmp.df, tmp.mat))
}