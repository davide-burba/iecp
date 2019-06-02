

#' reformat_dataset: add the following new columns
#' - start,stop: interval times as requested by survival::coxph
#' - Nm: number of events before the considered one for each individual N_i(t-)
#' - y: sum of the past marks for each individual
reformat_dataset = function(input,start_time = 0, verbose = FALSE){
  if(verbose){
    print('Adding new columns: start,stop, Nm [N_i(t-)], y [sum of past marks, if mark column is present]')
    }
  
  start = stop = Nm = y = c()
  input = data.table(input)
  individual_ids = unique(input$id)
  if(verbose){pb <- txtProgressBar(min = 0, max = length(individual_ids), style = 3)}
  
  for(i in 1:length(individual_ids)){
    if(verbose){setTxtProgressBar(pb, i)}
    
    individual_id = individual_ids[i]
    tmp = input[input$id == individual_id,]
    # compute times (start,stop]
    events_time = tmp[,'time_event']
    events_time = as.numeric(unlist(events_time))
    stop = c(stop,events_time)
    start = c(start,c(start_time,events_time[-length(events_time)]))
    # compute NM [N_i(t-)]
    Nm = c(Nm,0:(length(events_time)-1))
    # compute sum of past marks "y", if marks are present
    if ('mark' %in% colnames(input)){
        y_individual = tmp[,'mark']
        sum_y_individual = cumsum(as.numeric(unlist(y_individual)))
        y_individual = c(0,sum_y_individual[-length(sum_y_individual)])
        y = c(y, y_individual)
    }
  }
  if(verbose){cat("\n")}
  input[,'start' := start]
  input[,'stop' := stop]
  input[,'Nm' := Nm]
  if ('mark' %in% colnames(input)){input[,'y' := y]}
    
  # Keep only relevant features (discard "time_event" and "mark" to avoid confusion, if present)
  const_features = colnames(input)[!colnames(input) %in% c('id',
                                                           'time_event',
                                                           'status',
                                                           'mark',
                                                           'start',
                                                           'stop',
                                                           'Nm',
                                                           'y')]

  if ('mark' %in% colnames(input)){
      features = c('id','start','stop','status',const_features,'Nm','y')
  }else{
      features = c('id','start','stop','status',const_features,'Nm')
  }
  input = input[,..features]
  
  return(input)
}