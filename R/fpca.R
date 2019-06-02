
require(data.table)
require(ggplot2)
require(dplyr)
require(latex2exp)

# compute fpca scores of compensators, setting a threshold on explained variance
compute_fpca_scores = function(cumulative_hazard,
                               min_explained_variance = 0.99,
                               plot_eigenfunctions = FALSE, 
                               verbose = FALSE,
                               const_plot = 10,
                               cex=0.8,
                               y.intersp=2.5){
  
  cumulative_hazard = data.table(cumulative_hazard)
  times = sort(unique(cumulative_hazard$time))
  
  # reformat cumulative Hazard in a matrix format
  cumulative_hazard_matrix = .reformat_cumulative_hazard_to_matrix(cumulative_hazard,verbose)
  
  # Compute eigenvalues and eigenfunctions
  h <- (times[length(times)]-times[1])/(length(times)-1)
  pca <- prcomp(t(cumulative_hazard_matrix))
  efuncs <- pca$rotation*sqrt(1/h)
  evalues <- pca$sdev^2*h
  
  # Choose the minimum K such that the sum of eigenvalues (variances) is
  # at least min_explained_variance of total
  percent.var <- evalues/sum(evalues)
  K <- min(which(cumsum(percent.var)>=min_explained_variance))
  if(verbose){
    cat('\n')
    print(paste('To explain',min_explained_variance*100,'% of the variability',
                K,'components are needed.'))
  }
  
  # Fix eigenfunctions: sign chosen to have positive integral, for ease of interpretation
  for(i in 1:K) {
    tempfun <- efuncs[,i]
    tempsign <- sum(tempfun)
    efuncs[,i] <- ifelse(tempsign<0, -1,1) * tempfun
  }
    
  # plot eigenfunction
  mean_hazard = cumulative_hazard[,list('cumhaz' = mean(cumhaz)), by = 'time']
  if(plot_eigenfunctions){
    .plot_eigenfunctions(cumulative_hazard_matrix,
                         mean_hazard,
                         times,percent.var,
                         K,
                         efuncs,
                         const_plot,
                         cex,
                         y.intersp)
    }
  
  # Compute scores on the first K principal components
  score <- NULL
  for(i in 1:K){
    score <- cbind(score,
                   (as.matrix(t(cumulative_hazard_matrix-mean_hazard$cumhaz)) %*% cbind(efuncs[,i]))*h)
  }
  # reformat scores fpca
  rownames(score) = NULL
  column_names = NULL
  for(i in 1:K){column_names = c(column_names,paste('PC',i,sep = ''))}
  colnames(score) = column_names
  score = data.table('id'=colnames(cumulative_hazard_matrix) ,score)
  return(score)
}


# reformat cumulative Hazard in a matrix format
.reformat_cumulative_hazard_to_matrix = function(cumulative_hazard,verbose = FALSE){
  cumulative_hazard_matrix = list()
  individual_ids = unique(cumulative_hazard$id)
  if(verbose){pb <- txtProgressBar(min = 0, max = length(individual_ids), style = 3)}
  for(i in 1:length(individual_ids)){
    if(verbose){setTxtProgressBar(pb, i)}
    individual_id = individual_ids[i]
    cumulative_hazard_matrix[[individual_id]] = cumulative_hazard[id == individual_id,'cumhaz']
  }
  cumulative_hazard_matrix = bind_cols(cumulative_hazard_matrix)
  colnames(cumulative_hazard_matrix) = as.character(individual_ids)
  
  return(cumulative_hazard_matrix)
}


# Plot first K eigenfunctions
.plot_eigenfunctions = function(cumulative_hazard_matrix,
                                mean_hazard,
                                times,
                                percent.var,
                                K,
                                efuncs,
                                const_plot = 10,
                                cex=0.8,
                                y.intersp=2.5){
  
    m <- dim(cumulative_hazard_matrix)[1]
    mu <- apply(cumulative_hazard_matrix,MARGIN=1,"mean")
    par(mfcol=c(2,K),mar=c(2.5,2,2.5,1))
    ylim <- range(efuncs[,1:K])
    ylim2 <- c(0,1.3*max(mean_hazard$cumhaz))

    for(i in 1:K) {
      v1 <- efuncs[,i]
      plot(v1, type="l", ylab="", xlab="time", ylim=ylim, lwd=2, 
           main=paste("Component ", i,", ",100*round(percent.var[i],3), "%", sep=""))

      v2 <- v1 * const_plot
      widep <- ((1:m)%%30)==0
      plot(mu, type="l", ylim=ylim2, lwd=3)
      points(times[widep],mean_hazard$cumhaz[widep] + v2[widep], lwd=3, pch="+")
      points(times[widep],mean_hazard$cumhaz[widep] - v2[widep], lwd=3, pch="-")

      legend(min(mu),ylim2[2], # where
             legend=c(TeX('$\\mu(t)$'),
                      TeX(paste('$\\mu(t) + ',const_plot,'*PC_',i,'(t)$',sep = '')),
                      TeX(paste('$\\mu(t) - ',const_plot,'*PC_',i,'(t)$',sep = ''))), # text
             lty=c(1,NA,NA), lwd = c(2,1,1), pch = c(NA,'+','-'), # symbols
             cex=cex,y.intersp=y.intersp,box.lty=0 # other options
            )
    }
}