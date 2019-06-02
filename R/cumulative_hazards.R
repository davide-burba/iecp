require(data.table)
require(splines)
require(ggplot2)
require(latex2exp)
require(foreach)
require(doParallel)

#' Evaluate cumulative Hazard in a grid of points; returns a dataframe in long format
compute_cumulative_hazard = function(model,
                                     input,
                                     smoothed_baseline,
                                     times = NA,
                                     verbose = FALSE){
  input = data.table(input)
  # add missing columns, corresponding to operations between present columns
  input = .add_missing_columns(input,model)
  
  # build grid of time, if not given
  if(all(is.na(times))){times = seq(min(input$stop), max(input$stop),by = 1)}
  
  # Compute constant at times coefficients
  individual_coefficients = .compute_coefficients_ck(input,model, verbose)
    
  # Compute interval deltas of cumulative hazard and sum it up
  cumulative_hazard = .compute_cumulative_hazard(individual_coefficients,
                                                 times,
                                                 smoothed_baseline,
                                                 verbose)
  return(cumulative_hazard)
}


# Add missing columns, corresponding to operations between present columns
.add_missing_columns = function(input,model){
    name_coefficients = names(model$coefficients)
    # find operations between columns
    cols = name_coefficients[!name_coefficients %in% colnames(input)]
    # fix product operation (called ":" in formulas)
    cols_tmp_names = gsub(":", "*", cols)
    # if presents, add relative columns
    if (length(cols)>0){
        for(i in 1:length(cols)){
            col_name = cols[i]
            col_tmp_name = cols_tmp_names[i]
            val = eval(parse(text = paste0("values = input[ ,",  col_tmp_name, "]")))
            input[,(col_name) := val]
        }
    }
    return(input)
}


#' Compute coefficient ck = exp{beta*x_i(t_k)} at each interval for each individuals 
.compute_coefficients_ck = function(input,model, verbose = FALSE){
  if(verbose){print('Computing coefficients ck = exp{beta*x_i(t_k)}')}
  # take fitted coefficients beta of Cox model
  name_coefficients = names(model$coefficients)
  coefficients = model$coefficients
  
  individual_coefficients = NULL
  individuals_ids = unique(input$id)
  if(verbose){pb <- txtProgressBar(min = 0, max = length(individuals_ids), style = 3)}
  
  # compute coefficients for each individual 
  for(i in 1:length(individuals_ids)){
    if(verbose){setTxtProgressBar(pb, i)}
    individuals_id = individuals_ids[i]
    individuals_df = input[id == individuals_id,]
    ck = c()
    
    # compute coefficients for each time interval (start,stop]
    for(k in individuals_df$Nm){
      beta_times_xik = sum(individuals_df[, ..name_coefficients][k+1]*coefficients)
      ck = c(ck,exp(beta_times_xik))
    }
    individual_coefficients = rbind(individual_coefficients,
                                     data.frame(id = individuals_id,
                                                k = individuals_df$Nm,
                                                start = individuals_df$start,
                                                stop = individuals_df$stop,
                                                ck))
  }
  if(verbose){cat("\n")}
  return(data.table(individual_coefficients))
}


#' Compute daily deltas of cumulative hazard and sum it up; returns a dataframe in long format
.compute_cumulative_hazard = function(individual_coefficients,
                                      times,
                                      smoothed_baseline,
                                      verbose = FALSE){
    if(verbose){print('Computing cumulative Hazard on the grid')}
    
    # evaluate hazard at jump times
    Lambda0_starts = .Lambda0_fun(individual_coefficients$start,smoothed_baseline)
    Lambda0_stops = .Lambda0_fun(individual_coefficients$stop,smoothed_baseline)
    deltas = (Lambda0_stops - Lambda0_starts)*individual_coefficients$ck
    ids = unique(individual_coefficients$id)
    hazards = NULL
    for(id in ids){
        hazards = c(hazards, cumsum(deltas[individual_coefficients$id == id]))
    }
    individual_coefficients = data.table(individual_coefficients)#tmp
    individual_coefficients[,hazard :=hazards]
    individual_coefficients[,Lambda0_starts :=Lambda0_starts]
    Lamda0s = .Lambda0_fun(times,smoothed_baseline)

    # compute cumulative_hazard on a grid
    if(verbose){pb <- txtProgressBar(min = 0, max = length(ids), style = 3)}
    cumulative_hazard = NULL
    for( i in 1:length(ids)){
        individual = ids[i]
        if(verbose){setTxtProgressBar(pb, i)}
        id_df = individual_coefficients[id == individual]
        hazard_id = NULL

        for(i in 1:length(times)){
            t = times[i]
            Lambda0_t = Lamda0s[i]

            id_df_before = id_df[stop<=t] 
            id_df_after = id_df[stop>t]

            if(dim(id_df_before)[1]==0){
                starting_hazard = 0
            }else{
                starting_hazard = id_df_before$hazard[dim(id_df_before)[1]]
            }

            if(dim(id_df_after)[1]>0){
                lambda0_previous_start = id_df_after$Lambda0_starts[1]
                ck = id_df_after$ck[1]
                hazard_id_t = starting_hazard + ck*(Lambda0_t-lambda0_previous_start)
            }else{
                hazard_id_t = starting_hazard 
            }
            hazard_id = c(hazard_id,hazard_id_t)
        }
        cumulative_hazard = rbind(cumulative_hazard,
                                  data.frame(id = individual,time = times,cumhaz = hazard_id))
    }
    if(verbose){cat("\n")}
    return(cumulative_hazard)
}


# Lambda_0
#      t  = evaluation points
.Lambda0_fun <- function(t,smoothed_baseline){
  return(.basis(t, smoothed_baseline$knots) %*% smoothed_baseline$coef)
}


########################################################################################
### basis: returns the evaluations of a b-spline basis on a vector of points x
###
### Inputs:
###     x:  vector of abscissa values for evaluation
### knots:  knots of the basis
###   deg:  polynomial degree of the basis
###
### Outputs:
###         a matrix with the evaluations in values of x (rows) of the B-splines basis
###         functions (columns)
########################################################################################
.basis <- function(x, knots, deg=2){
  return(bs(x,
            knots=knots[2:(length(knots)-1)],
            degree=deg,
            Boundary.knots=c(knots[1],knots[length(knots)]),
            intercept=TRUE))
}


# Plot compensators (sample)
plot_cumulative_hazard = function(cumulative_hazard, sample_size = NA){
    
    if(is.na(sample_size)){sample_size = min(500,length(unique(cumulative_hazard$id)))}
    
    # select a sample to plot
    sample_individualss = sample(unique(cumulative_hazard$id),sample_size)
    sample = cumulative_hazard[cumulative_hazard$id %in% sample_individualss,]
    
    # plot
        ggplot(sample, aes(x= time, y=cumhaz, group = factor(id), color=factor(id),alpha = 0.1)) +
            geom_line() +
            xlab('time') + 
            ylab(TeX("$\\hat{\\Lambda}_i(t)")) + 
            theme(legend.position="none",plot.title = element_text(hjust = 0.5))
}