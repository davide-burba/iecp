
require(ggplot2)
require(latex2exp)
require(data.table)

# compute residuals
get_residuals = function(input,cumulative_hazard,times = NA,verbose = FALSE){
    
    if(all(is.na(times))){times = seq(min(input$stop), max(input$stop),by = 1)}
    realizations_long_format = .compute_realizations_long_format(input,times,verbose)
    residuals = cbind(cumulative_hazard, Nt = realizations_long_format$Nt)
    residuals['residuals'] = residuals$cumhaz - residuals$Nt
    residuals = data.table(residuals)
    return(residuals)
}


# compute dataframe of realizations in long format 
.compute_realizations_long_format = function(input,times = NA, verbose = FALSE){
  individual_ids = unique(input$id)
  if(all(is.na(times))){times = seq(min(input$stop), max(input$stop),by = 1)}
    
  daily_realizations = NULL
  if(verbose){pb <- txtProgressBar(min = 0, max = length(individual_ids), style = 3)}
  # could be parallelised on individuals
  for(i in 1:length(individual_ids)){ 
    if(verbose){setTxtProgressBar(pb, i)}
    individual_id = individual_ids[i]
    tmp = input[id==individual_id]
    Nt = c()
    for(t in times){
      nt = tmp[start<=t,Nm]
      nt = nt[length(nt)]
      Nt = c(Nt,nt)
    }
    daily_realizations = rbind(daily_realizations,data.frame(individual_id,time = times, Nt))
  }
  if(verbose){cat("\n")}
  return(daily_realizations)
}  


# Plot martingale Residuals
plot_martingale_residuals = function(residuals,sample_size = NA){
    # compute time-varying mean residuals
    residuals = data.table(residuals)
    mean_residuals = residuals[,list('residuals' = mean(residuals)), by = 'time']
      
    # select a sample to plot
    if(is.na(sample_size)){sample_size = min(500,length(unique(residuals$id)))}
    sample_individuals = sample(unique(residuals$id),sample_size)
    residuals = data.frame(residuals)
    sample = residuals[residuals$id %in% sample_individuals,]
    
    # plot
    ggplot(data = mean_residuals, aes(x= time, y=residuals)) +
        geom_line() +
        geom_line(data = sample, 
                  aes(x= time, y=residuals, group = factor(id), color=factor(id),alpha = 0.1)) +
        geom_line(data = mean_residuals, aes(x= time, y=residuals)) +
        xlab('time') + 
        ylab(TeX("$\\hat{\\M}_i(t) = \\hat{\\Lambda}_i(t) - N_i(t)")) + 
        theme(legend.position="none",plot.title = element_text(hjust = 0.5))
}

