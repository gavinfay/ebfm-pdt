### MSPROD equation with system cap
## Solves the multsipecies operating model dynamics for a single time step given parameters, a set of harvest rates, and the current biomass
## Given a maximum catch on the system
dNbydt <- function(t,X=c(1,0),parms=list(r=rep(0.4,1),
                                         KGuild=rep(2,1),
                                         Ktot=10,
                                         alpha=matrix(0,nrow=1,ncol=1),
                                         Guildmembership=1,
                                         BetweenGuildComp=matrix(0,nrow=1,ncol=1),
                                         WithinGuildComp=matrix(0,nrow=1,ncol=1),
                                         hrate=0,maxcat=0)) {
  Nsp <- length(X)/2
  N <- X[1:Nsp]
  testcat <- sum(parms$hrate*N,na.rm=TRUE)
  frac <- 1
  if (testcat>parms$maxcat) frac <- parms$maxcat/testcat
  hrate <- parms$hrate*frac
  NG <- aggregate(N,by=list(parms$Guildmembership),sum,na.rm=TRUE)
  NG <- t(parms$BetweenGuildComp)%*%NG$x
  dN <- parms$r*N*(1-(N/parms$KGuild[parms$Guildmembership])-(t(parms$WithinGuildComp)%*%N)/parms$KGuild[parms$Guildmembership]-NG[parms$Guildmembership]/(parms$Ktot-parms$KGuild[parms$Guildmembership]))- parms$alpha%*%N*N-hrate*N
  #dN <- pmax(rep(0.0001,length(N)),r*N*(1-(N/KGuild[Guildmembership])-(t(WithinGuildComp)%*%N)/KGuild[Guildmembership]-NG[Guildmembership]/(Ktot-KGuild[Guildmembership]))- alpha%*%N*N-hrate*N)
  cat <- hrate*N
  predloss <-  parms$alpha%*%N*N
  betweenloss <- parms$r*N*NG[parms$Guildmembership]/(parms$Ktot-parms$KGuild[parms$Guildmembership])
  withinloss <- parms$r*N*(parms$WithinGuildComp%*%N)/parms$KGuild[parms$Guildmembership]
  results <- list(deriv=c(dN,cat),predloss=predloss,withinloss=withinloss,betweenloss=betweenloss)
  return(results)
}

get_om_pars <- function() {
  pars = NULL
  pars$Nsp <- 4
  pars$Nyr <- 50
  pars$Fuse <- c(0.01	,
            0.04	,
            0.07	,
            0.1	,
            0.13	,
            0.16	,
            0.19	,
            0.22	,
            0.25	,
            0.28	,
            0.31	,
            0.34	,
            0.37	,
            0.4	,
            0.43	,
            0.46	,
            0.49	,
            0.52	,
            0.55	,
            0.58	,
            0.61	,
            0.64	,
            0.67	,
            0.7	,
            0.7	,
            0.7	,
            0.7	,
            0.7	,
            0.7	,
            0.7	,
            0.7	,
            0.7	,
            0.7	,
            0.7	,
            0.7	,
            0.7	,
            0.7	,
            0.7	,
            0.7	,
            0.7	,
            0.65	,
            0.6	,
            0.55	,
            0.5	,
            0.45	,
            0.4	,
            0.35	,
            0.3	,
            0.25,
            0.2)
  pars$r <- c(0.7, 0.4, 0.25, 0.5)
  pars$K <- rep(1000, pars$Nsp)
  pars$q <- c(0.8,0.5,0.2,0.65)
  pars$InitBio <- rep(1000, pars$Nsp)
  pars$complex <- c(1,2,2,1)
  return(pars)
}

get_om <- function(pars) {
  iyr = 1
  Nsp <- length(pars$r)
  parms = list(
    r = pars$r,
    KGuild = pars$K,
    Ktot = sum(pars$K, na.rm = TRUE),
    alpha = matrix(0, nrow = pars$Nsp, ncol = pars$Nsp),
    Guildmembership = 1:pars$Nsp,
    BetweenGuildComp = matrix(0, nrow =
                                pars$Nsp, ncol = pars$Nsp),
    WithinGuildComp = matrix(0, nrow =
                               pars$Nsp, ncol = pars$Nsp),
    q = pars$q,
    hrate = pars$q * 0,
    maxcat = 1e+06
  )
  init_vec <- c(pars$InitBio,rep(0,length(pars$InitBio)))
  names(init_vec) <- c(paste0("biomass_",1:pars$Nsp),paste0("catch_",1:pars$Nsp))
  XX <- purrr::accumulate(pars$Fuse, om_update, .init = init_vec, parms = parms)
  return(XX)
}

om_update <- function(X, Fstar, parms) {
  iyr <- 1
  Nsp <- length(parms$r)
  parms$hrate <- Fstar*parms$q
  #X <- c(N[1:Nsp],rep(0,Nsp))
  x <- deSolve::ode(X,seq(iyr-1,(iyr),1),dNbydt,parms=parms,method="rk4")
  om_vec <- x[2,2:(1+2*Nsp)]
  om_vec[(1+Nsp):(2*Nsp)] <- om_vec[(1+Nsp):(2*Nsp)] - x[1,(2+Nsp):(1+2*Nsp)]
  names(om_vec) <- c(paste0("biomass_",1:Nsp),paste0("catch_",1:Nsp))
  return(om_vec)
}


run_om <- function(input) {
 
  #pars <- get_om_pars() 
  om <- get_om(input)
  names(om) <- 1:length(om)
  #om_results <- bind_cols(om)
  #names(om_results) <- c(paste0("biomass_",1:pars$Nsp),paste0("catch_",1:pars$Nsp)) 
  om_long <- do.call(bind_rows, om) %>% 
    rowid_to_column() %>% 
    pivot_longer(values_to = "value", names_to = "variable", -rowid) %>% 
    separate(col = "variable", into = c("type", "isp")) %>% 
    rename("t" = rowid) %>% 
    mutate(complex = ifelse(isp %in% c(1,4), 1, 2))
  return(om_long)
}

gen_data <- function(om_long) {
  set.seed(365)
  sigma_b <- 0.35
  sigma_c <- 0.15
  em_data <- om_long %>% 
    mutate(data_sd = ifelse(type == "biomass",sigma_b,sigma_c),
           err = rnorm(nrow(om_long),0,data_sd),
           data = value*exp(err-0.5*data_sd*data_sd)) %>% 
    mutate(t = ifelse(type=="catch",t-1,t)) %>% 
    filter(t != 0) %>% 
    mutate(isp = as.numeric(isp)) %>% 
    I()
  #em_data
  ss_data <- em_data %>%
    select(t, isp, type, data) %>% 
    ungroup() %>% 
    I()
  complex_data <- em_data %>%
    select(t, isp, type, data, complex) %>%
    group_by(t, complex, type) %>% 
    summarize(data = sum(data)) %>% 
    ungroup() %>% 
    mutate(isp = complex+4) %>% 
    select(-complex) %>% 
    I()
  system_data <- em_data %>% 
    group_by(t, type) %>% 
    summarize(data = sum(data)) %>% 
    ungroup() %>% 
    mutate(isp = rep(7,nrow(.))) %>% 
    I()
  gendata <- bind_rows(ss_data, complex_data, system_data)
#  gendata <- list(ss_data = ss_data,
 #                  complex_data = complex_data,
  #                 system_data = system_data)
  return(gendata)
}

do_assess <- function(data) {
  ini.parms <- log(c(5000, 0.6)) #log(1200), log(0.1), 0.3)
  # Fit the logistic model to data:  
  biomass <- as_vector(data$biomass)
  catch <- as_vector(data$catch)
  results <- assess(catch,biomass,calc.vcov=FALSE,ini.parms)
  # results <- optim(c(0.6, -0.001), get_preds, 
  #       biomass = as_vector(data$biomass),
  #       catch = as_vector(data$catch),
  #       method = "BFGS")
  return(results)
}

get_preds <- function(beta, biomass, catch) {
  pred_asp <- beta[1]*biomass + beta[2]*(biomass^2)
  obs_asp <- catch + diff(biomass)
  rss <- sum((obs_asp - pred_asp)^2,na.rm=TRUE)
  return(rss)
}


run_assessments <- function(om_long) {
  emdata <- gen_data(om_long) %>% 
    group_by(isp) %>% 
    pivot_wider(names_from = type,
                values_from = data) %>% 
    arrange(t) %>% 
    nest() %>% 
    I()
  
  
  results <- emdata %>% 
    mutate(results = map(data, do_assess),
           pars = map(results, "pars"),
           bmsy = exp(map_dbl(pars, 1))/2,
           msy = exp(map_dbl(pars, 2))*bmsy/2,
           fmsy = exp(map_dbl(pars, 2))/2,
           q = map_dbl(results, "q")) %>% 
    I()
#  print(results)

  # # Extract the maximum likelihood and parameter estimates  
  # biomass.mle <- redfish$biomass
  # print(biomass.mle)
  # pars.mle <- redfish$pars
  # print(pars.mle)
  
  
  
  
}

  
hline_plot <- function(point=FALSE) {
  list(if(point)
    geom_hline(aes(yintercept = q*bmsy), col = "blue", lty =2))
}  
floor_line <- function(point=FALSE) {
  list(if(point)
    geom_hline(aes(yintercept = floor), col = "red"))
}  
floor_vline <- function(point=FALSE) {
  list(if(point)
    geom_vline(aes(xintercept = floor/bmsy), col = "red"))
}  


surveyplot <- function(assess_results, settings, input) {
                       
    # settings = list(
    # #showTimeSeries = input$showTimeSeries,
    # useCeiling = "Yes",
    # assessType = "single species")
    # #targetF = input$targetF,
    # #floorB = input$floorB,
    # #floorOption = input$floorOption
    # 
  
    complex_id <- tibble(isp = 1:input$Nsp,
                       complex = input$complex)
  emdata_plot <- assess_results %>% 
    ungroup() %>% 
    filter(isp <= input$Nsp) %>% 
    select(isp, data, bmsy, q) %>% 
    unnest(cols = c(data)) %>% 
    pivot_longer(cols = c(biomass,catch), names_to = "type", values_to = "value") %>% 
    left_join(complex_id) %>% 
    I()
  floors <- emdata_plot %>% 
    group_by(isp, type) %>% 
    summarize(floor = settings$floorB*max(value, na.rm = TRUE)) %>% 
    ungroup() %>% 
    I()
  emdata_plot <- left_join(emdata_plot, floors)
  emdata_plot
  
  p1 <- emdata_plot %>% 
    filter(type == "biomass") %>% 
    ggplot() +
    aes(x = t, y = value) +
    geom_point() +
    hline_plot(settings$assessType == "single species") +
    floor_line(settings$assessType == "stock complex") +
    facet_wrap(~fct_reorder(factor(isp), complex), scales = "free") +
    labs(y = "survey index",
         x = "year",
#         title = "New plot <b style='color:#009E73'>title</b>", 
#         subtitle = "A <b style='color:#D55E00'>subtitle</b>") +
         title = "survey index: species",
         subtitle = "with estimated <b style='color:blue;'>BMSY</b> & <b style='color:red;'>species floors</b>") +
    theme(plot.title = element_markdown(lineheight = 1.1),
          plot.subtitle = element_markdown(lineheight = 1.1)) +
    NULL

  if (settings$showTimeSeries == "kobe plots") {
  kobe_dat <- assess_results %>% 
    ungroup() %>% 
    mutate(bio = map(results, "biomass"),
           bio = map(bio, ~.[-length(.)]),
           cat = map(data, "catch"),
           f = map2(cat,bio,function(x,y) x/y),
           bbmsy = map(bio, ~(./bmsy)),
           ffmsy = map(f,~(./fmsy)),
           obs = map(data, "biomass"),
           t = map(data, "t"),
           obs = map(obs, ~(./q/bmsy))) %>% 
    left_join(filter(floors, type == "biomass")) %>% 
    select(isp, t, obs, bbmsy, ffmsy, bmsy, floor) %>% 
    unnest(cols = -c(isp, bmsy, floor)) %>% 
    I()
  #kobe_dat
  
  p1 <- kobe_dat %>% 
    filter(isp <= input$Nsp) %>%
    left_join(complex_id) %>% 
    ggplot() +
    aes(x = obs, y = ffmsy) +
    geom_point() +
    geom_line(aes(x = bbmsy, y = ffmsy),col="darkolivegreen") +
    geom_hline(yintercept = 1, lty = 2) +
    geom_vline(xintercept = 1, lty = 2, col = "blue") +
    floor_vline(settings$assessType == "stock complex") +
    #hline_plot(settings$assessType == "single species") +
    #floor_line(settings$assessType == "stock complex") +
    facet_wrap(~fct_reorder(factor(isp), complex), scales = "free") +
    xlim(0,2) +
    ylim(0,2) +
    labs(y = "F/FMSY",
         x = "B/BMSY",
         #         title = "New plot <b style='color:#009E73'>title</b>", 
         #         subtitle = "A <b style='color:#D55E00'>subtitle</b>") +
         title = "Assessment results: species") +
         #subtitle = "with estimated <b style='color:blue;'>BMSY</b> & <b style='color:red;'>species floors</b>") +
    theme(plot.title = element_markdown(lineheight = 1.1)) + #,
          #plot.subtitle = element_markdown(lineheight = 1.1)) +
    NULL
  }
  
  syspp <- max(assess_results$isp)
  emdata_plot <- assess_results %>% 
    ungroup() %>% 
    filter(isp > input$Nsp, isp < syspp) %>% 
    select(isp, data, bmsy, q) %>% 
    unnest(cols = c(data)) %>% 
    ungroup() %>% 
    mutate(isp = isp - input$Nsp) %>% 
    pivot_longer(cols = c(biomass,catch), names_to = "type", values_to = "value") %>% 
    #left_join(complex_id) %>% 
    I()
  p2 <- emdata_plot %>% 
    filter(type == "biomass") %>% 
    ggplot() +
    aes(x = t, y = value) +
    geom_point() +
    hline_plot(settings$assessType == "stock complex") +
    facet_wrap(~factor(isp), ncol = 1, scales = "free") +
    labs(x = "year",
         title = "survey index: complexes") +
    NULL
  emdata_plot <- assess_results %>% 
    ungroup() %>% 
    filter(isp == syspp) %>% 
    select(isp, data, bmsy, q) %>% 
    unnest(cols = c(data)) %>% 
    ungroup() %>% 
    pivot_longer(cols = c(biomass,catch), names_to = "type", values_to = "value") %>% 
    #left_join(complex_id) %>% 
    I()
  p3 <- emdata_plot %>% 
    filter(type == "biomass") %>% 
    ggplot() +
    aes(x = t, y = value) +
    geom_point() +
    hline_plot(settings$useCeiling == "Yes") +
    #facet_wrap(~factor(isp), ncol = 1, scales = "free") +
    labs(x = "year",
         title = "survey index: system") +
    NULL  
  p1 + p2 + (p3/plot_spacer()) +
    plot_layout(widths = c(2, 1, 1))
}    
  
# biomass <- NULL
# biomass[[1]] <- em_data %>% 
#   filter(isp == 1, type == "biomass") %>% 
#   select(data) %>% 
#   as_vector()
# catch <- NULL
# catch[[1]] <- em_data %>% 
#   filter(isp == 1, type == "catch") %>% 
#   select(data) %>% 
#   as_vector() 
# xx <- map2(biomass, catch, ~do_assess(biomass = .x, catch = .y))
# 






schaefer <- function(B,C,K,r) {
  #function schaefer takes the current biomass, a catch, 
  #and the model parameters to compute next year's biomass
  res <- B + B * r * (1 - B/K) - C
  return(max(0.001,res))  # we add a constraint to prevent negative biomass
}

# Now a function to do the biomass projection:  
dynamics <- function(pars,C,yrs) {
  # dynamics takes the model parameters, the time series of catch, 
  # & the yrs to do the projection over
  
  # first extract the parameters from the pars vector (we estimate K in log-space)
  K <- exp(pars[1])
  r <- exp(pars[2])
  
  # find the total number of years
  nyr <- length(C) + 1
  
  # if the vector of years was not supplied we create 
  # a default to stop the program crashing
  if (missing(yrs)) yrs <- 1:nyr
  
  #set up the biomass vector
  B <- numeric(nyr)
  
  #intialize biomass at carrying capacity
  B[1] <- K
  # project the model forward using the schaefer model
  for (y in 2:nyr) {
    B[y] <- schaefer(B[y-1],C[y-1],K,r)
  }
  
  #return the time series of biomass
  return(B[yrs])
  
  #end function dynamics
}  


# We are going to condition the operating model by estimating the parameters based on the historical biomass index data.  
# 
# To do this we make a function that shows how well the current parameters fit the data, we assume that the observation errors around the true biomass are log-normally distributed.  

# function to calculate the negative log-likelihood
nll <- function(pars,C,U,getq=FALSE) {  #this function takes the parameters, the catches, and the index data
  #sigma <- exp(pars[3])  # additional parameter, the standard deviation of the observation error
  B <- dynamics(pars,C)  #run the biomass dynamics for this set of parameters
  B <- B[-length(B)]
  #Uhat <- exp(pars[4])*B   #calculate the predicted biomass index - here we assume an unbiased absolute biomass estimate
  #print(U)
  #print(B)
  q = exp(sum(log(U/B))/length(U))
  Uhat = q*B
  sigma = sqrt(sum((log(U)-log(q)-log(B))^2, na.rm = TRUE)/length(U))
  
 # obj_fun += ncpue*log(f_sigma) + 0.5*ncpue;
  output <- length(U)*log(sigma) + length(U)/2
#  output <- -sum(dnorm(log(U),log(Uhat),sigma,log=TRUE),na.rm=TRUE)   #calculate the negative log-likelihood
  ifelse(getq==TRUE,return(q),return(output))
  #end function nll
}


# Function to perform the assessment and estimate the operating model parameters  
# (i.e. to fit the logistic model to abundance data)
assess <- function(catch,index,calc.vcov=FALSE,pars.init) {
  # assess takes catch and index data, initial values for the parameters,
  # and a flag saying whether to compute uncertainty estimates for the model parameters
  
  #fit model
  # optim runs the function nll() repeatedly with differnt values for the parameters,
  # to find the values that give the best fit to the index data
  res <- optim(pars.init,nll,C=catch,U=index,getq = FALSE, method = "BFGS") #,hessian=TRUE)
  
  # store the output from the model fit
  output <- list()
  output$pars <- res$par
  output$q <- nll(res$par, C=catch, U=index, getq = TRUE)
  output$biomass <- dynamics(res$par,catch)
  output$convergence <- res$convergence
  output$nll <- res$value
  if (calc.vcov)
    output$vcov <- solve(res$hessian)
  
  return(output)
  #end function assess
}



ref_table <- function(settings, assess_results, input) {
  sysspp <- input$Nsp + max(input$complex, na.rm=TRUE) + 1
  refs <- assess_results %>% 
    ungroup() %>% 
    mutate(estbio = map(results, "biomass"),
           estbio = map(estbio, ~.[-length(.)]),
           blast = map_dbl(estbio, ~.[length(.)]),
           cfmsy = blast * fmsy,
           ffmsy = map_dbl(data, ~.$catch[length(.$catch)-1]),
           ffmsy = ffmsy/map_dbl(estbio, ~.[length(.)-1])/fmsy,
           bbmsy = blast/bmsy,
           ddmax = map_dbl(data, ~last(.$biomass))/map_dbl(data, ~max(.$biomass, na.rm = TRUE)),
           bfloor = ddmax/settings$floorB,
           bfloor = ifelse(bfloor>1,1,bfloor),
           ftarg = settings$targetF*fmsy,
           fuse = ftarg) %>% 
    I()
    if(settings$assessType == "stock complex")
      refs <- mutate(refs, fuse = ftarg*bfloor)
    refs <- mutate(refs, cfuse = blast*fuse) %>% 
    I()
  ceiling_val <- refs$cfmsy[sysspp]
  refs <- refs  %>% mutate(ceiling = rep(ceiling_val, nrow(.))) %>%
    I()
  #refs
  complex_res <- refs %>% 
    ungroup() %>% 
    filter(isp <= input$Nsp) %>% 
    mutate(complex = input$complex) %>% 
    group_by(complex) %>% 
    summarize(floor_mean = mean(bfloor),
              floor_min = min(bfloor)) %>% 
      # input$floorOption == "avg status" ~ mean(bfloor),
      # input$floorOption == "min status" ~ min(bfloor))) %>%  #,
      #TRUE ~ NA)) %>% 
   #ifelse(input$floorOption == "avg status", 
    I()
  #complex_res

  if (settings$assessType == "single species") {
    ss_maxcat <- refs %>% 
      filter(isp <= input$Nsp) %>% 
      summarize(totcat = sum(cfuse),
                ceiling = mean(ceiling)) %>% 
      mutate(ceiling_mult = ifelse(totcat/ceiling>1,ceiling/totcat,1)) %>% 
      I()
    if (settings$useCeiling == "No") ss_maxcat$ceiling_mult <- 1
    out_table <- refs %>% 
      mutate(advice = ss_maxcat$ceiling_mult*cfuse) %>% 
      select(isp, fmsy, bmsy, blast, ffmsy, bbmsy, cfmsy, fuse, cfuse, ceiling, advice) %>% 
      #select(isp, fmsy, bmsy, blast, ffmsy, bbmsy, cfmsy) %>% 
      filter(isp <= input$Nsp)
    if(settings$useCeiling == "No")
       out_table <- select(out_table, -ceiling)
  }

  
  if (settings$assessType == "stock complex") {
    #refs$bfloor
    sc_refs <- refs %>%
      ungroup() %>% 
      filter(isp > input$Nsp & isp < sysspp) %>% 
      mutate(bfloor = case_when(
        settings$floorOption == "avg status" ~ complex_res$floor_mean,
        settings$floorOption == "min status" ~ complex_res$floor_min),
        fuse = ftarg*bfloor,
        cfuse = blast*fuse) %>% 
      I()
    #sc_refs$bfloor
    
    sc_maxcat <- sc_refs %>% 
      summarize(totcat = sum(cfuse),
                ceiling = mean(ceiling)) %>% 
      mutate(ceiling_mult = ifelse(totcat/ceiling>1,ceiling/totcat,1)) %>% 
      I()
    if (settings$useCeiling == "No") {sc_maxcat$ceiling_mult <- 1}
    out_table <- sc_refs %>% 
      mutate(advice = sc_maxcat$ceiling_mult*cfuse,
             complex = isp - input$Nsp) %>% 
      select(complex, fmsy, bmsy, blast, ffmsy, bbmsy, cfmsy, bfloor, fuse, cfuse, ceiling, advice)
      #select(isp, fmsy, bmsy, blast, ffmsy, bbmsy, cfmsy) %>% 
    if(settings$useCeiling == "No")
      out_table <- select(out_table, -ceiling)
  }
  
  return(out_table)  
  
}



##################
### Shiny App Main tab functions
#################
do_ts_plot <- function(settings) {
  input <- get_om_pars()
  #input$Fuse <- rep(0.25,length(input$Fuse))
  #input$K <- rep(2000,4)
  om_long <- run_om(input)
  assess_results <- run_assessments(om_long)
  surveyplot(assess_results, settings, input)
}


report_table <- function(settings) {
  input <- get_om_pars()
  om_long <- run_om(input)
  assess_results <- run_assessments(om_long)
  out_table <- ref_table(settings, assess_results, input) %>% 
    #remove_rownames() %>%
    #ungroup() %>% 
    as_tibble()
  return(out_table)
}


settings = list(
  showTimeSeries = "No",
  useCeiling = "Yes",
  assessType = "stock complex",
  targetF = 0.75,
  floorB = 0.25,
  floorOption = "avg status")
# do_ts_plot(settings)
# # print(report_table(settings))
# # dimnames(report_table(settings))
