
sim_ghost = function(years, error, phi, p){
  
  phi = rep(phi, years)

  lambda = -log(1-p)
  
  marked = rep(n.ind, years-1)
  n.occasions = years
  
  CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
  CH.n = CH
  
  mark.occ <- rep(1:length(marked), marked[1:length(marked)])
  
  for (i in 1:sum(marked)){
    CH[i, mark.occ[i]] <- 1
    CH.n[i, mark.occ[i]] <- 1
    if (mark.occ[i]==n.occasions) next
    for (t in (mark.occ[i]+1):n.occasions){
      sur <- rbinom(1, 1, phi[t])
      if (sur==0) break		
      rp <- rbinom(1, 1, p)
      if (rp==1){
        CH[i,t] <- 1
        
        # if seen, how many times?
        nobs = max(rpois(1, lambda), 1)
        CH.n[i,t] <- nobs
      }
    } 
  } 
  
  # loop through occasions, introduce misreads by creating ghosts
  # if misID occurs, new CH created an appended
  # ghosts can only be seen once
  
  CH.ghost = matrix(0, nrow = 1, ncol = ncol(CH))
  
  for(i in 1:ncol(CH)){
    ch <- CH.n[,i]
    
    # detected individuals - each detection gives potential misread
    x = which(ch > 0)
    misreads = rbinom(rep(1, length(x)), ch[x], error)
    
    # if the number of misreads is >= number of obs, switch that ind to not detected
    nd = which(misreads >= ch[x])
    CH[nd,i] <- 0
    
    n.misreads = sum(misreads)
    
    # create ghost capture histories
    new.ghosts = matrix(0, nrow = n.misreads, ncol = ncol(CH))
    new.ghosts[,i] <- 1
    
    CH.ghost = rbind(CH.ghost, new.ghosts)
  }
  
  # remove first all-zero row of ghost ch
  CH.ghost = CH.ghost[-1,]
  
  # bind ghosts to real CH
  CH = rbind(CH, CH.ghost)

  # check for and remove inds never seen
  ns = which(rowSums(CH) == 0)
  if(length(ns > 0)){
    CH = CH[-ns,]
  }

  # # remove first observations
  CH.clean = CH.n
  first = apply(CH.clean, 1, function(x) min(which(x != 0)))
  for(i in 1:nrow(CH.clean)){
    CH.clean[i,first[i]] <- CH.clean[i,first[i]] - 1
  }
  CH.clean[CH.clean > 1] <- 1
  
  # check for and remove inds never seen
  ns = which(rowSums(CH.clean) == 0)
  if(length(ns > 0)){
    CH.clean = CH.clean[-ns,]
  }
  
  # format for marked - no data cleaning
  dat = data.frame(ch = apply(CH, 1, function(x) paste(x, collapse = "")))
  dat$ch = as.character(dat$ch)

  # with data cleaning
  dat_clean = data.frame(ch = apply(CH.clean, 1, function(x) paste(x, collapse = "")))
  dat_clean$ch = as.character(dat_clean$ch)

  # fit both data sets
  Phi.t = list(formula = ~time)
  p.dot = list(formula = ~1)
  
  results1 = mark(data = dat, ddl = NULL,
                  model.parameters = list(Phi = Phi.t, p = p.dot))
  results2 <- mark(data = dat_clean, ddl = NULL,
                  model.parameters = list(Phi = Phi.t, p = p.dot))
  return(list(results1, results2))
}
