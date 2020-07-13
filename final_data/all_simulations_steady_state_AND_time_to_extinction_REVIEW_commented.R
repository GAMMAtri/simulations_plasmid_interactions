# NOTES: ----
# Simulations were run using R version 3.4.4 (2018-03-15) -- "Someone to Lean On", on a Platform: x86_64-pc-linux-gnu (64-bit) 
# with processor: 32-core processor × 64; memory: 125,8 GiB; os: Ubuntu 18.04.4 LTS
# WARNING: the notation for the equation parameters below differs from that in the manuscript

# libraries ----
library(deSolve); library(reshape2); library(ggplot2); library(gridExtra); library(grid); library(openxlsx); library(data.table)
library(parallel); library(foreach); library(doParallel); library(rootSolve)

# simulations concerning two plasmids (includes steady state and time to extinction) ----
# parameters ----
gama = c(1e-13,1e-12,1e-11) # conjugation efficiency
intra = c(1e-3,1e0,1e1) # intracellular effect on conjugation
inter = c(1e-2,1e0,1e1) # intercellular effect on conjugation
loss = 10^(seq(from = -8, to = -4, by = 2)) # loss rate
cost = c(seq(from = .85, to = .95, by = .05),.975) # relative fitness
loss2 = c(1,10) # effect on loss
epi = seq(from = -.05, to = .1, by = .05) # epistasis

# making all combinations of the variables
tb = expand.grid(gama,gama,intra,intra,inter,inter,cost,cost,epi,loss,loss2)

# adding ID for each plasmid and combination
# name of each plasmid p1 and p2 (ie, X and Y) is in the format: 
# "conjugation efficiency"_"intracellular_conjugation_effect"_"intercellular_conjugation_effect"_"fitness"_"epistastic_effect"_"loss"_"loss_interaction_effect"
# name of the combination has the format "plasmid1":"plasmid2" (ie, X:Y)
tmp = tb; for (i in 1:ncol(tmp)) {tmp[,i] = as.factor(tmp[,i])}
tb['p1'] = apply(X = tmp, MARGIN = 1, FUN = function(x) { paste(x[c("Var1","Var3","Var5","Var7","Var9","Var10","Var11")], collapse = '_' )})
tb['p2'] = apply(X = tmp, MARGIN = 1, FUN = function(x) { paste(x[c("Var2","Var4","Var6","Var8","Var9","Var10","Var11")], collapse = '_' )})
tb['comb'] = apply(X = tb[,c('p1','p2')], MARGIN = 1, FUN = function(x){ paste(x, collapse = ':')})
colnames(tb) = c('alone1','alone2','intra1fx','intra2fx','inter1fx','inter2fx','fit_alone1','fit_alone2','epis','loss','loss2','p1','p2','comb')

# precalculating conjugation, fitness and loss in function of the interaction
tb['intra1'] = tb$alone1 * tb$intra1fx; tb['intra2'] = tb$alone2 * tb$intra2fx # intracellular conjugation efficiency of plasmid 1 (X); and plasmid 2 (Y)
tb['inter1'] = tb$alone1 * tb$inter1fx; tb['inter2'] = tb$alone2 * tb$inter2fx # intercellular conjugation efficiency of plasmid 1 (X); and plasmid 2 (Y)
tb['fitness'] = (tb$fit_alone1 * tb$fit_alone2) + tb$epis # fitness of doubles XY
tb['cotransf'] = apply(X = tb[,c('intra1','intra2')], MARGIN = 1, FUN = min) # contransfer of X and Y from doubles
tb['exx1fx'] = ifelse(test = tb$intra1fx > 1, yes = tb$alone1, no = tb$intra1) # entry exclusion effect for plasmid 1 (X)
tb['exx2fx'] = ifelse(test = tb$intra2fx > 1, yes = tb$alone2, no = tb$intra2) # entry exclusion effect for plasmid 2 (Y)
tb['lossXY'] = tb$loss * tb$loss2 # combined loss
rownames(tb) = tb$comb; View(tb)

# subsetting (and keeping only) cases where the two plasmids do not interact, or there is only one type of interaction
tmp = tb[(tb$intra1fx + tb$intra2fx + tb$inter1fx + tb$inter2fx + tb$loss2) == 5 & tb$epis == 0,]; paste(unique(tmp$intra1fx),unique(tmp$intra2fx),unique(tmp$inter1fx),unique(tmp$inter2fx),unique(tmp$loss2),unique(tmp$epis)) # no interactions
tmp2 = tb[(tb$intra1fx + tb$intra2fx + tb$inter1fx + tb$inter2fx + tb$loss2) == 5 & tb$epis != 0,]; paste(unique(tmp2$intra1fx),unique(tmp2$intra2fx),unique(tmp2$inter1fx),unique(tmp2$inter2fx),unique(tmp2$loss2),unique(tmp2$epis)) # epistasis
tmp3 = tb[(tb$intra2fx + tb$inter1fx + tb$inter2fx + tb$loss2) == 4 & tb$epis == 0 & tb$intra1fx != 1,]; paste(unique(tmp3$intra1fx),unique(tmp3$intra2fx),unique(tmp3$inter1fx),unique(tmp3$inter2fx),unique(tmp3$loss2),unique(tmp3$epis)) # intrafx1
tmp4 = tb[(tb$intra1fx + tb$inter1fx + tb$inter2fx + tb$loss2) == 4 & tb$epis == 0 & tb$intra2fx != 1,]; paste(unique(tmp4$intra1fx),unique(tmp4$intra2fx),unique(tmp4$inter1fx),unique(tmp4$inter2fx),unique(tmp4$loss2),unique(tmp4$epis)) # intrafx2
tmp5 = tb[(tb$intra1fx + tb$intra2fx + tb$inter2fx + tb$loss2) == 4 & tb$epis == 0 & tb$inter1fx != 1,]; paste(unique(tmp5$intra1fx),unique(tmp5$intra2fx),unique(tmp5$inter1fx),unique(tmp5$inter2fx),unique(tmp5$loss2),unique(tmp5$epis)) # interfx1
tmp6 = tb[(tb$intra1fx + tb$intra2fx + tb$inter1fx + tb$loss2) == 4 & tb$epis == 0 & tb$inter2fx != 1,]; paste(unique(tmp6$intra1fx),unique(tmp6$intra2fx),unique(tmp6$inter1fx),unique(tmp6$inter2fx),unique(tmp6$loss2),unique(tmp6$epis)) # interfx2
tmp7 = tb[(tb$intra1fx + tb$intra2fx + tb$inter1fx + tb$inter2fx) == 4 & tb$epis == 0 & tb$loss2 != 1,]; paste(unique(tmp7$intra1fx),unique(tmp7$intra2fx),unique(tmp7$inter1fx),unique(tmp7$inter2fx),unique(tmp7$loss2),unique(tmp7$epis)) # loss2
tb0 = tb
tb = rbind(tmp,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7)

# steady state ----
# running the simulation in parallel/multiple cores for each row of conditions in the table tb generated above
registerDoParallel(detectCores()*3/4)
system.time({
  tab = foreach (i = 1:nrow(tb), .inorder = F) %do% { # saving the results of the simulation for steady state in the (dummy) list "tab"
    ini = c(O = 0, X = 5e5, Y = 5e5, XY = 0, R = 100) # initial conditions: 5x10^5 cells X & 5 x10^5 cells Y only (no Ø nor XY cells) and initial resources = 100
    
    parms = c(Ro = 100, D = 1e-6, W = .05, Q = 5, G = 3 * log(2), # pre-set parameters:
              # resource concentration, resource required per cell division, turnover, Monod, growth rate
              
              s = tb[i,'loss'], s2 = tb[i,'lossXY'], # segregation rates for each plasmid X and Y respectivelly
              hX = tb[i,'alone1'], hY = tb[i,'alone2'], # conjugation rates alone for each plasmid X and Y respectivelly
              hXy = tb[i,'intra1'], hYx = tb[i,'intra2'], # intracellular conjugation rates for each plasmid X and Y respectivelly
              hX_Y = tb[i,'inter1'], hY_X = tb[i,'inter2'], # intercellular conjugation rates for each plasmid X and Y respectivelly
              hXy2 = tb[i,'exx1fx'], hYx2 = tb[i,'exx2fx'], # exclusion on conjugation rates for each plasmid X and Y respectivelly
              hXY = tb[i,'cotransf'], # co-transfer rate
              FX = tb[i,'fit_alone1'], FY = tb[i,'fit_alone2'], FXY = tb[i,'fitness']) # fitness of cells X, Y and XY
    
    # defining the model
    model <- function (time, y, parms) {
      with(as.list(c(y, parms)), {
        
        dO = O*G*(R/(R+Q)) - W*O + s*(X+Y) - hX*O*X - hY*O*Y - (hXY + (hXy-hXY) + (hYx-hXY))*O*XY # Ø cells
        
        dX =  X*G*FX*(R/(R+Q)) - W*X + s2*XY - s*X + hX*O*X + (hXy-hXY)*O*XY - hY_X*X*Y - hYx2*X*XY # X cells  
        
        dY =  Y*G*FY*(R/(R+Q)) - W*Y + s2*XY - s*Y + hY*O*Y + (hYx-hXY)*O*XY - hX_Y*X*Y - hXy2*Y*XY # Y cells 
        
        dXY =  XY*G*FXY*(R/(R+Q)) - W*XY - 2*s2*XY + hXY*O*XY + (hX_Y+hY_X)*Y*X + hYx2*X*XY + hXy2*Y*XY # XY cells
        
        dR =  W*(Ro-R) - D*G*(R/(R+Q))* (O + X*FX + Y*FY + XY*FXY) # resources
        
        list(c(dO,dX,dY,dXY,dR))
      })
    }
    
    tmp = runsteady(y = ini, fun = model, parms = parms, times = c(0, 1e10), maxsteps = 1e9) # running the model for steady state
    
    # appeding precision, steady, time and number os steps required to reach teh steady state, combination
    cbind(data.frame(t(tmp$y),t(unlist(attributes(tmp)[c('precis','steady','time','steps')]))), data.frame(ID = rownames(tb[i,]))) 
  } 
})
stopImplicitCluster()

tab_cmplt = rbindlist(tab) # reshaping the results list into a single dataframe (tab_cmplt)

unique(tab_cmplt$steady) # if 1, all simulations reached the stead state
nrow(tab_cmplt[tab_cmplt$X > 1 | tab_cmplt$XY > 1,]); View(tab_cmplt[tab_cmplt$X > 1 | tab_cmplt$XY > 1,]) # cases where plasmid X is stably maintained
# checking maximum time, steps and precision of the siulations
max(as.numeric(as.character(tab_cmplt$time))); max(as.numeric(as.character(tab_cmplt$steps))); max(as.numeric(as.character(tab_cmplt$precis))) 

write.table(x = tab_cmplt, file = 'Downloads/gama_steady_state_singles_expanded_all_REVIEW.txt') # saving the results of the steady state simulations as a tab separated file

# saving the table parameters used as a tab separated file
tab_cmplt$ID = as.character(tab_cmplt$ID)
tb['time'] = apply(X = as.matrix(tb$comb), MARGIN = 1, FUN = function(x) { ceiling(unlist(tab_cmplt[tab_cmplt$ID == x, 'time'])) } )
write.table(x = tb, file = 'Downloads/gama_parameters_singles_expanded_all_REVIEW.txt')




# timewise simulation for extinctions ----
# subsetting the parameter table to keep only the cases of extinction to be able to calculate their time to extinction
tb0 = tb
tb = tb[tb$comb %in% unlist(tab_cmplt[tab_cmplt$X < 1 & tab_cmplt$XY < 1, "ID"]), ]
fivenum(tb$time) # distribution of times required for steady state

# running the model to get times to extinction: the model and conditions are the same used above
# running the simulation in parallel/multiple cores for each row of conditions in the table tb generated above
registerDoParallel(detectCores()*3/4)
system.time({
  tab2 =  # saving the results of the simulation for times of extinction in the (dummy) list "tab2"
    foreach (i = 1:nrow(tb), .inorder = F) %do% {
      ini = c(O = 0, X = 5e5, Y = 5e5, XY = 0, R = 100)  # initial conditions: 5x10^5 cells X & 5 x10^5 cells Y only (no Ø nor XY cells) and initial resources = 100
      
      parms = c(Ro = 100, D = 1e-6, W = .05, Q = 5, G = 3 * log(2), # pre-set parameters:
                # resource concentration, resource required per cell division, turnover, Monod, growth rate
                
                s = tb[i,'loss'], s2 = tb[i,'lossXY'],# segregation rates for each plasmid X and Y respectivelly
                hX = tb[i,'alone1'], hY = tb[i,'alone2'], # conjugation rates alone for each plasmid X and Y respectivelly
                hXy = tb[i,'intra1'], hYx = tb[i,'intra2'], # intracellular conjugation rates for each plasmid X and Y respectivelly
                hX_Y = tb[i,'inter1'], hY_X = tb[i,'inter2'], # intercellular conjugation rates for each plasmid X and Y respectivelly
                hXy2 = tb[i,'exx1fx'], hYx2 = tb[i,'exx2fx'], # exclusion on conjugation rates for each plasmid X and Y respectivelly
                hXY = tb[i,'cotransf'], # co-transfer rate
                FX = tb[i,'fit_alone1'], FY = tb[i,'fit_alone2'], FXY = tb[i,'fitness']) # fitness of cells X, Y and XY
      
      times <- seq(0, tb[i,'time'], 1) # cycling times hourly from 0 up to the time of steady state specific for each combination
      
      # defining the model
      model <- function (time, y, parms) {
        with(as.list(c(y, parms)), {
          
          dO = O*G*(R/(R+Q)) - W*O + s*(X+Y) - hX*O*X - hY*O*Y - (hXY + (hXy-hXY) + (hYx-hXY))*O*XY # Ø cells
          
          dX =  X*G*FX*(R/(R+Q)) - W*X + s2*XY - s*X + hX*O*X + (hXy-hXY)*O*XY - hY_X*X*Y - hYx2*X*XY # X cells 
          
          dY =  Y*G*FY*(R/(R+Q)) - W*Y + s2*XY - s*Y + hY*O*Y + (hYx-hXY)*O*XY - hX_Y*X*Y - hXy2*Y*XY # Y cells
          
          dXY =  XY*G*FXY*(R/(R+Q)) - W*XY - 2*s2*XY + hXY*O*XY + (hX_Y+hY_X)*Y*X + hYx2*X*XY + hXy2*Y*XY # XY cells
          
          dR =  W*(Ro-R) - D*G*(R/(R+Q))* (O + X*FX + Y*FY + XY*FXY) # resources
          
          list(c(dO,dX,dY,dXY,dR))
        })
      }
      
      tmp = data.frame(ode(ini, times, model, parms)) # running the model for the times specified above (in the vector "times")
      
      # keeping the results as the time of extinction for X, XY, X or XY (X_XY), Y and the ID of the parameters combination
      # time of extinction is the first time where the entity amount is < 1
      cbind(data.frame(X = min(tmp[tmp$X < 1, 'time']), XY = min(tmp[tmp$XY < 1 & tmp$time > min(tmp[tmp$XY > 1, 'time']), 'time']),
                       X_XY = min(tmp[tmp$X < 1 & tmp$XY < 1, 'time']), Y = min(tmp[tmp$Y < 1, 'time']), ID = rownames(tb[i,])))

    } 
})
stopImplicitCluster()

tab_cmplt2 = rbindlist(tab2) # reshaping the results list into a single dataframe (tab_cmplt2)

# checking if there is any non natural value in the final table (ie, any error)
any(is.na(tab_cmplt2$X),is.na(tab_cmplt2$XY),is.na(tab_cmplt2$X_XY),is.na(tab_cmplt2$Y))
any(is.infinite(tab_cmplt2$X),is.infinite(tab_cmplt2$XY),is.infinite(tab_cmplt2$X_XY),is.infinite(tab_cmplt2$Y))

write.table(x = tab_cmplt2, file = 'Downloads/gama_time_to_extinction_singles_expanded_all_REVIEW.txt') # saving the results of the time of extinction simulations as a tab separated file











# simulations concerning two plasmids (includes steady state and time to extinction) ----
# there is a fixed competitor (Y) and multiple interactions are allowed ----
# parameters ----
gama = c(1e-13,1e-12,1e-11) # conjugation efficiency of X
gama2 = 1e-12 # conjugation efficiency of Y
intra = c(1e-3,1e0,1e1) # intracellular effect on conjugation
inter = c(1e-2,1e0,1e1) # intercellular effect on conjugation
loss = 10^(seq(from = -8, to = -4, by = 2)) # loss rate
cost = c(seq(from = .85, to = .95, by = .05),.975) # relative fitness of X
cost2 = .9 # relative fitness of Y
epi = seq(from = -.05, to = .1, by = .05) # epistasis
loss2 = c(1,10) # effect on loss

# making all combinations of the variables
tb = expand.grid(gama,gama2,intra,intra,inter,inter,cost,cost2,epi,loss,loss2)

# adding ID for each plasmid and combination
# name of each plasmid p1 and p2 (ie, X and Y) is in the format: 
# "conjugation efficiency"_"intracellular_conjugation_effect"_"intercellular_conjugation_effect"_"fitness"_"epistastic_effect"_"loss"_"loss_interaction_effect"
# name of the combination has the format "plasmid1":"plasmid2" (ie, X:Y)
tmp = tb; for (i in 1:ncol(tmp)) {tmp[,i] = as.factor(tmp[,i])}
tb['p1'] = apply(X = tmp, MARGIN = 1, FUN = function(x) { paste(x[c("Var1","Var3","Var5","Var7","Var9","Var10","Var11")], collapse = '_' )})
tb['p2'] = apply(X = tmp, MARGIN = 1, FUN = function(x) { paste(x[c("Var2","Var4","Var6","Var8","Var9","Var10","Var11")], collapse = '_' )})
tb['comb'] = apply(X = tb[,c('p1','p2')], MARGIN = 1, FUN = function(x){ paste(x, collapse = ':')})
colnames(tb) = c('alone1','alone2','intra1fx','intra2fx','inter1fx','inter2fx','fit_alone1','fit_alone2','epis','loss','loss2','p1','p2','comb')

# precalculating conjugation, fitness and loss in function of the interaction
tb['intra1'] = tb$alone1 * tb$intra1fx; tb['intra2'] = tb$alone2 * tb$intra2fx # intracellular conjugation efficiency of plasmid 1 (X); and plasmid 2 (Y)
tb['inter1'] = tb$alone1 * tb$inter1fx; tb['inter2'] = tb$alone2 * tb$inter2fx # intercellular conjugation efficiency of plasmid 1 (X); and plasmid 2 (Y)
tb['fitness'] = (tb$fit_alone1 * tb$fit_alone2) + tb$epis # fitness of doubles XY
tb['cotransf'] = apply(X = tb[,c('intra1','intra2')], MARGIN = 1, FUN = min) # contransfer of X and Y from doubles
tb['exx1fx'] = ifelse(test = tb$intra1fx > 1, yes = tb$alone1, no = tb$intra1)  # entry exclusion effect for plasmid 1 (X)
tb['exx2fx'] = ifelse(test = tb$intra2fx > 1, yes = tb$alone2, no = tb$intra2)  # entry exclusion effect for plasmid 2 (Y)
tb['lossXY'] = tb$loss * tb$loss2 # combined loss
rownames(tb) = tb$comb; View(tb)


# steady state ----
# running the simulation in parallel/multiple cores for each row of conditions in the table tb generated above; the model and conditions are the same used above
registerDoParallel(detectCores()*3/4)
system.time({
  tab = foreach (i = 1:nrow(tb), .inorder = F) %do% { # saving the results of the simulation for steady state in the (dummy) list "tab"
    ini = c(O = 0, X = 5e5, Y = 5e5, XY = 0, R = 100) # initial conditions: 5x10^5 cells X & 5 x10^5 cells Y only (no Ø nor XY cells) and initial resources = 100
    
    parms = c(Ro = 100, D = 1e-6, W = .05, Q = 5, G = 3 * log(2), # pre-set parameters:
              # resource concentration, resource required per cell division, turnover, Monod, growth rate
              
              s = tb[i,'loss'], s2 = tb[i,'lossXY'],# # segregation rates for each plasmid X and Y respectivelly
              hX = tb[i,'alone1'], hY = tb[i,'alone2'], # conjugation rates alone for each plasmid X and Y respectivelly
              hXy = tb[i,'intra1'], hYx = tb[i,'intra2'], # intracellular conjugation rates for each plasmid X and Y respectivelly
              hX_Y = tb[i,'inter1'], hY_X = tb[i,'inter2'], # intercellular conjugation rates for each plasmid X and Y respectivelly
              hXy2 = tb[i,'exx1fx'], hYx2 = tb[i,'exx2fx'], # exclusion on conjugation rates for each plasmid X and Y respectivelly
              hXY = tb[i,'cotransf'], # co-transfer rate
              FX = tb[i,'fit_alone1'], FY = tb[i,'fit_alone2'], FXY = tb[i,'fitness']) # fitness of cells X, Y and XY
    
    # defining the model
    model <- function (time, y, parms) {
      with(as.list(c(y, parms)), {
        
        dO = O*G*(R/(R+Q)) - W*O + s*(X+Y) - hX*O*X - hY*O*Y - (hXY + (hXy-hXY) + (hYx-hXY))*O*XY # Ø cells
        
        dX =  X*G*FX*(R/(R+Q)) - W*X + s2*XY - s*X + hX*O*X + (hXy-hXY)*O*XY - hY_X*X*Y - hYx2*X*XY # X cells
        
        dY =  Y*G*FY*(R/(R+Q)) - W*Y + s2*XY - s*Y + hY*O*Y + (hYx-hXY)*O*XY - hX_Y*X*Y - hXy2*Y*XY # Y cells
        
        dXY =  XY*G*FXY*(R/(R+Q)) - W*XY - 2*s2*XY + hXY*O*XY + (hX_Y+hY_X)*Y*X + hYx2*X*XY + hXy2*Y*XY # XY cells
        
        dR =  W*(Ro-R) - D*G*(R/(R+Q))* (O + X*FX + Y*FY + XY*FXY) # resources
        
        list(c(dO,dX,dY,dXY,dR))
      })
    }
    
    tmp = runsteady(y = ini, fun = model, parms = parms, times = c(0, 1e10), maxsteps = 1e9) # running the model for steady state
    
    # appeding precision, steady, time and number os steps required to reach teh steady state, combination
    cbind(data.frame(t(tmp$y),t(unlist(attributes(tmp)[c('precis','steady','time','steps')]))), data.frame(ID = rownames(tb[i,])))
  } 
})
stopImplicitCluster()

tab_cmplt = rbindlist(tab) # reshaping the results list into a single dataframe (tab_cmplt)

unique(tab_cmplt$steady) # if 1, all simulations reached the stead state
nrow(tab_cmplt[tab_cmplt$X > 1 | tab_cmplt$XY > 1,]); View(tab_cmplt[tab_cmplt$X > 1 | tab_cmplt$XY > 1,]) # cases where plasmid X is stably maintained
# checking maximum time, steps and precision of the siulations
max(as.numeric(as.character(tab_cmplt$time))); max(as.numeric(as.character(tab_cmplt$steps))); max(as.numeric(as.character(tab_cmplt$precis)))

write.table(x = tab_cmplt, file = 'Downloads/gama_steady_state_singles_fixed_competitor_REVIEW.txt') # saving the results of the steady state simulations as a tab separated file

# saving the table parameters used as a tab separated file
tab_cmplt$ID = as.character(tab_cmplt$ID)
tb['time'] = apply(X = as.matrix(tb$comb), MARGIN = 1, FUN = function(x) { ceiling(unlist(tab_cmplt[tab_cmplt$ID == x, 'time'])) } )
write.table(x = tb, file = 'Downloads/gama_parameters_singles_fixed_competitor_REVIEW.txt')




# timewise simulation for extinctions ----
# subsetting the parameter table to keep only the cases of extinction to be able to calculate their time to extinction
tb0 = tb
tb = tb[tb$comb %in% unlist(tab_cmplt[tab_cmplt$X < 1 & tab_cmplt$XY < 1, "ID"]), ]
fivenum(tb$time) # distribution of times required for steady state

# running the model to get times to extinction: the model and conditions are the same used above
# running the simulation in parallel/multiple cores for each row of conditions in the table tb generated above
registerDoParallel(detectCores()*3/4)
system.time({
  tab2 = # saving the results of the simulation for times of extinction in the (dummy) list "tab2" 
    foreach (i = 1:nrow(tb), .inorder = F) %do% {
      ini = c(O = 0, X = 5e5, Y = 5e5, XY = 0, R = 100) # initial conditions: 5x10^5 cells X & 5 x10^5 cells Y only (no Ø nor XY cells) and initial resources = 100
      
      parms = c(Ro = 100, D = 1e-6, W = .05, Q = 5, G = 3 * log(2), # pre-set parameters:
                # resource concentration, resource required per cell division, turnover, Monod, growth rate
                
                s = tb[i,'loss'], s2 = tb[i,'lossXY'], # segregation rates for each plasmid X and Y respectivelly
                hX = tb[i,'alone1'], hY = tb[i,'alone2'], # conjugation rates alone for each plasmid X and Y respectivelly
                hXy = tb[i,'intra1'], hYx = tb[i,'intra2'], # intracellular conjugation rates for each plasmid X and Y respectivelly
                hX_Y = tb[i,'inter1'], hY_X = tb[i,'inter2'], # intercellular conjugation rates for each plasmid X and Y respectivelly
                hXy2 = tb[i,'exx1fx'], hYx2 = tb[i,'exx2fx'], # exclusion on conjugation rates for each plasmid X and Y respectivelly
                hXY = tb[i,'cotransf'], # co-transfer rate
                FX = tb[i,'fit_alone1'], FY = tb[i,'fit_alone2'], FXY = tb[i,'fitness']) # fitness of cells X, Y and XY
      
      times <- seq(0, tb[i,'time'], 1) # cycling times hourly from 0 up to the time of steady state specific for each combination
      
      # defining the model
      model <- function (time, y, parms) {
        with(as.list(c(y, parms)), {
          
          dO = O*G*(R/(R+Q)) - W*O + s*(X+Y) - hX*O*X - hY*O*Y - (hXY + (hXy-hXY) + (hYx-hXY))*O*XY # Ø cells
          
          dX =  X*G*FX*(R/(R+Q)) - W*X + s2*XY - s*X + hX*O*X + (hXy-hXY)*O*XY - hY_X*X*Y - hYx2*X*XY # X cells
          
          dY =  Y*G*FY*(R/(R+Q)) - W*Y + s2*XY - s*Y + hY*O*Y + (hYx-hXY)*O*XY - hX_Y*X*Y - hXy2*Y*XY # Y cells 
          
          dXY =  XY*G*FXY*(R/(R+Q)) - W*XY - 2*s2*XY + hXY*O*XY + (hX_Y+hY_X)*Y*X + hYx2*X*XY + hXy2*Y*XY # XY cells
          
          dR =  W*(Ro-R) - D*G*(R/(R+Q))* (O + X*FX + Y*FY + XY*FXY) # resources
          
          list(c(dO,dX,dY,dXY,dR))
        })
      }
      
      tmp = data.frame(ode(ini, times, model, parms)) # running the model for the times specified above (in the vector "times")
      
      # keeping the results as the time of extinction for X, XY, X or XY (X_XY), Y and the ID of the parameters combination
      # time of extinction is the first time where the entity amount is < 1
      cbind(data.frame(X = min(tmp[tmp$X < 1, 'time']), XY = min(tmp[tmp$XY < 1 & tmp$time > min(tmp[tmp$XY > 1, 'time']), 'time']),
                       X_XY = min(tmp[tmp$X < 1 & tmp$XY < 1, 'time']), Y = min(tmp[tmp$Y < 1, 'time']), ID = rownames(tb[i,])))
      
    } 
})
stopImplicitCluster()

tab_cmplt2 = rbindlist(tab2) # reshaping the results list into a single dataframe (tab_cmplt2)

# checking if there is any non natural value in the final table (ie, any error)
any(is.na(tab_cmplt2$X),is.na(tab_cmplt2$XY),is.na(tab_cmplt2$X_XY),is.na(tab_cmplt2$Y))
any(is.infinite(tab_cmplt2$X),is.infinite(tab_cmplt2$XY),is.infinite(tab_cmplt2$X_XY),is.infinite(tab_cmplt2$Y))

write.table(x = tab_cmplt2, file = 'Downloads/gama_time_to_extinction_singles_fixed_competitor_REVIEW.txt') # saving the results of the time of extinction simulations as a tab separated file









# simulations concerning only a single plasmid (includes steady state and time to extinction) ----
# parameters ----
# interactions are set as 1, therefore there is no effect from plasmid Y
gama = c(1e-13,1e-12,1e-11,1e-10) # conjugation efficiency of X
gama2 = 0 # conjugation efficiency of Y
intra = 1 # intracellular effect on conjugation
inter = 1 # intercellular effect on conjugation
loss = 10^(seq(from = -8, to = -4, by = 2)) # loss rate
cost = c(seq(from = .85, to = .95, by = .05),.975,1,1.05) # relative fitness of X
cost2 = 0 # relative fitness of Y
epi = 0 # epistasis
loss2 = 1 # effect on loss

# making all combinations of the variables
tb = expand.grid(gama,gama2,intra,intra,inter,inter,cost,cost2,epi,loss,loss2)

# adding ID for each plasmid and combination
# name of each plasmid p1 and p2 (ie, X and Y) is in the format: 
# "conjugation efficiency"_"intracellular_conjugation_effect"_"intercellular_conjugation_effect"_"fitness"_"epistastic_effect"_"loss"_"loss_interaction_effect"
# name of the combination has the format "plasmid1":"plasmid2" (ie, X:Y)
tmp = tb; for (i in 1:ncol(tmp)) {tmp[,i] = as.factor(tmp[,i])}
tb['p1'] = apply(X = tmp, MARGIN = 1, FUN = function(x) { paste(x[c("Var1","Var3","Var5","Var7","Var9","Var10","Var11")], collapse = '_' )})
tb['p2'] = apply(X = tmp, MARGIN = 1, FUN = function(x) { paste(x[c("Var2","Var4","Var6","Var8","Var9","Var10","Var11")], collapse = '_' )})
tb['comb'] = apply(X = tb[,c('p1','p2')], MARGIN = 1, FUN = function(x){ paste(x, collapse = ':')})
colnames(tb) = c('alone1','alone2','intra1fx','intra2fx','inter1fx','inter2fx','fit_alone1','fit_alone2','epis','loss','loss2','p1','p2','comb')

# precalculating conjugation, fitness and loss in function of the interaction
tb['intra1'] = tb$alone1 * tb$intra1fx; tb['intra2'] = tb$alone2 * tb$intra2fx # intracellular conjugation efficiency of plasmid 1 (X); and plasmid 2 (Y)
tb['inter1'] = tb$alone1 * tb$inter1fx; tb['inter2'] = tb$alone2 * tb$inter2fx # intercellular conjugation efficiency of plasmid 1 (X); and plasmid 2 (Y)
tb['fitness'] = (tb$fit_alone1 * tb$fit_alone2) + tb$epis # fitness of doubles XY
tb['cotransf'] = apply(X = tb[,c('intra1','intra2')], MARGIN = 1, FUN = min) # contransfer of X and Y from doubles
tb['exx1fx'] = ifelse(test = tb$intra1fx > 1, yes = tb$alone1, no = tb$intra1)  # entry exclusion effect for plasmid 1 (X)
tb['exx2fx'] = ifelse(test = tb$intra2fx > 1, yes = tb$alone2, no = tb$intra2)  # entry exclusion effect for plasmid 2 (Y)
tb['lossXY'] = tb$loss * tb$loss2 # combined loss
rownames(tb) = tb$comb; View(tb)

tb = tb[!(tb$alone1 == 1e-10 & tb$fit_alone1 >= 1), ]

# steady state ----
# running the simulation in parallel/multiple cores for each row of conditions in the table tb generated above; the model and conditions are the same used above
registerDoParallel(detectCores()*3/4)
system.time({
  tab = foreach (i = 1:nrow(tb), .inorder = F) %do% { # saving the results of the simulation for steady state in the (dummy) list "tab"
    ini = c(O = 0, X = 5e5, Y =0, XY = 0, R = 100) # initial conditions: 5x10^5 cells X only (no Ø, Y nor XY cells) and initial resources = 100
    
    parms = c(Ro = 100, D = 1e-6, W = .05, Q = 5, G = 3 * log(2), # pre-set parameters:
              # resource concentration, resource required per cell division, turnover, Monod, growth rate
              
              s = tb[i,'loss'], s2 = tb[i,'lossXY'],# # segregation rates for each plasmid X and Y respectivelly
              hX = tb[i,'alone1'], hY = tb[i,'alone2'], # conjugation rates alone for each plasmid X and Y respectivelly
              hXy = tb[i,'intra1'], hYx = tb[i,'intra2'], # intracellular conjugation rates for each plasmid X and Y respectivelly
              hX_Y = tb[i,'inter1'], hY_X = tb[i,'inter2'], # intercellular conjugation rates for each plasmid X and Y respectivelly
              hXy2 = tb[i,'exx1fx'], hYx2 = tb[i,'exx2fx'], # exclusion on conjugation rates for each plasmid X and Y respectivelly
              hXY = tb[i,'cotransf'], # co-transfer rate
              FX = tb[i,'fit_alone1'], FY = tb[i,'fit_alone2'], FXY = tb[i,'fitness']) # fitness of cells X, Y and XY
    
    # defining the model
    model <- function (time, y, parms) {
      with(as.list(c(y, parms)), {
        
        dO = O*G*(R/(R+Q)) - W*O + s*(X+Y) - hX*O*X - hY*O*Y - (hXY + (hXy-hXY) + (hYx-hXY))*O*XY # Ø cells
        
        dX =  X*G*FX*(R/(R+Q)) - W*X + s2*XY - s*X + hX*O*X + (hXy-hXY)*O*XY - hY_X*X*Y - hYx2*X*XY # X cells
        
        dY =  Y*G*FY*(R/(R+Q)) - W*Y + s2*XY - s*Y + hY*O*Y + (hYx-hXY)*O*XY - hX_Y*X*Y - hXy2*Y*XY # Y cells
        
        dXY =  XY*G*FXY*(R/(R+Q)) - W*XY - 2*s2*XY + hXY*O*XY + (hX_Y+hY_X)*Y*X + hYx2*X*XY + hXy2*Y*XY # XY cells
        
        dR =  W*(Ro-R) - D*G*(R/(R+Q))* (O + X*FX + Y*FY + XY*FXY) # resources
        
        list(c(dO,dX,dY,dXY,dR))
      })
    }
    
    tmp = runsteady(y = ini, fun = model, parms = parms, times = c(0, 1e10), maxsteps = 1e9) # running the model for steady state
    
    # appeding precision, steady, time and number os steps required to reach teh steady state, combination
    cbind(data.frame(t(tmp$y),t(unlist(attributes(tmp)[c('precis','steady','time','steps')]))), data.frame(ID = rownames(tb[i,])))
  } 
})
stopImplicitCluster()

tab_cmplt = rbindlist(tab) # reshaping the results list into a single dataframe (tab_cmplt)
unique(tab_cmplt$steady) # if 1, all simulations reached the stead state
nrow(tab_cmplt[tab_cmplt$X > 1,]); View(tab_cmplt[tab_cmplt$X > 1,]) # cases where plasmid X is stably maintained
# checking maximum time, steps and precision of the siulations
max(as.numeric(as.character(tab_cmplt$time))); max(as.numeric(as.character(tab_cmplt$steps))); max(as.numeric(as.character(tab_cmplt$precis)))

write.table(x = tab_cmplt, file = 'Downloads/gama_steady_state_solo_REVIEW.txt') # saving the results of the steady state simulations as a tab separated file

# saving the table parameters used as a tab separated file
tab_cmplt$ID = as.character(tab_cmplt$ID)
tb['time'] = apply(X = as.matrix(tb$comb), MARGIN = 1, FUN = function(x) { ceiling(unlist(tab_cmplt[tab_cmplt$ID == x, 'time'])) } )
write.table(x = tb, file = 'Downloads/gama_parameters_solo_REVIEW.txt')




# timewise simulation for extinctions ----
# subsetting the parameter table to keep only the cases of extinction to be able to calculate their time to extinction
tb0 = tb
tb = tb[tb$comb %in% unlist(tab_cmplt[tab_cmplt$X < 1 & tab_cmplt$XY < 1, "ID"]), ]
fivenum(tb$time) # distribution of times required for steady state

# running the model to get times to extinction: the model and conditions are the same used above
# running the simulation in parallel/multiple cores for each row of conditions in the table tb generated above
registerDoParallel(detectCores()*3/4)
system.time({
  tab2 = # saving the results of the simulation for times of extinction in the (dummy) list "tab2" 
    foreach (i = 1:nrow(tb), .inorder = F) %do% {
      ini = c(O = 0, X = 5e5, Y =0, XY = 0, R = 100) # initial conditions: 5x10^5 cells X only (no Ø, Y nor XY cells) and initial resources = 100
      
      parms = c(Ro = 100, D = 1e-6, W = .05, Q = 5, G = 3 * log(2), # pre-set parameters:
                # resource concentration, resource required per cell division, turnover, Monod, growth rate
                
                s = tb[i,'loss'], s2 = tb[i,'lossXY'], # segregation rates for each plasmid X and Y respectivelly
                hX = tb[i,'alone1'], hY = tb[i,'alone2'], # conjugation rates alone for each plasmid X and Y respectivelly
                hXy = tb[i,'intra1'], hYx = tb[i,'intra2'], # intracellular conjugation rates for each plasmid X and Y respectivelly
                hX_Y = tb[i,'inter1'], hY_X = tb[i,'inter2'], # intercellular conjugation rates for each plasmid X and Y respectivelly
                hXy2 = tb[i,'exx1fx'], hYx2 = tb[i,'exx2fx'], # exclusion on conjugation rates for each plasmid X and Y respectivelly
                hXY = tb[i,'cotransf'], # co-transfer rate
                FX = tb[i,'fit_alone1'], FY = tb[i,'fit_alone2'], FXY = tb[i,'fitness']) # fitness of cells X, Y and XY
      
      times <- seq(0, tb[i,'time'], 1) # cycling times hourly from 0 up to the time of steady state specific for each combination
      
      # defining the model
      model <- function (time, y, parms) {
        with(as.list(c(y, parms)), {
          
          dO = O*G*(R/(R+Q)) - W*O + s*(X+Y) - hX*O*X - hY*O*Y - (hXY + (hXy-hXY) + (hYx-hXY))*O*XY # Ø cells
          
          dX =  X*G*FX*(R/(R+Q)) - W*X + s2*XY - s*X + hX*O*X + (hXy-hXY)*O*XY - hY_X*X*Y - hYx2*X*XY # X cells
          
          dY =  Y*G*FY*(R/(R+Q)) - W*Y + s2*XY - s*Y + hY*O*Y + (hYx-hXY)*O*XY - hX_Y*X*Y - hXy2*Y*XY # Y cells 
          
          dXY =  XY*G*FXY*(R/(R+Q)) - W*XY - 2*s2*XY + hXY*O*XY + (hX_Y+hY_X)*Y*X + hYx2*X*XY + hXy2*Y*XY # XY cells
          
          dR =  W*(Ro-R) - D*G*(R/(R+Q))* (O + X*FX + Y*FY + XY*FXY) # resources
          
          list(c(dO,dX,dY,dXY,dR))
        })
      }
      
      tmp = data.frame(ode(ini, times, model, parms)) # running the model for the times specified above (in the vector "times")
      
      # keeping the results as the time of extinction for X, and the ID of the parameters combination
      # time of extinction is the first time where the entity amount is < 1
      cbind(data.frame(X = min(tmp[tmp$X < 1, 'time']), ID = rownames(tb[i,])))
    } 
})
stopImplicitCluster()

tab_cmplt2 = rbindlist(tab2) # reshaping the results list into a single dataframe (tab_cmplt2)

# checking if there is any non natural value in the final table (ie, any error)
any(is.na(tab_cmplt2$X))
any(is.infinite(tab_cmplt2$X))

write.table(x = tab_cmplt2, file = 'Downloads/gama_time_to_extinction_solo_REVIEW.txt') # saving the results of the time of extinction simulations as a tab separated file








# simulations concerning only a single plasmid (includes steady state and time to extinction) ----
# this is a control for the initial density
# parameters ----
# interactions are set as 1, therefore there is no effect from plasmid Y
gama = c(1e-13,1e-12,1e-11) # conjugation efficiency of X
gama2 = 0 # conjugation efficiency of Y
intra = 1 # intracellular effect on conjugation
inter = 1 # intercellular effect on conjugation
loss = 10^(seq(from = -8, to = -4, by = 2)) # loss rate
cost = c(seq(from = .85, to = .95, by = .05),.975) # relative fitness of X
cost2 = 0 # relative fitness of Y
epi = 0 # epistasis
loss2 = 1 # effect on loss

# making all combinations of the variables
tb = expand.grid(gama,gama2,intra,intra,inter,inter,cost,cost2,epi,loss,loss2)

# adding ID for each plasmid and combination
# name of each plasmid p1 and p2 (ie, X and Y) is in the format: 
# "conjugation efficiency"_"intracellular_conjugation_effect"_"intercellular_conjugation_effect"_"fitness"_"epistastic_effect"_"loss"_"loss_interaction_effect"
# name of the combination has the format "plasmid1":"plasmid2" (ie, X:Y)
tmp = tb; for (i in 1:ncol(tmp)) {tmp[,i] = as.factor(tmp[,i])}
tb['p1'] = apply(X = tmp, MARGIN = 1, FUN = function(x) { paste(x[c("Var1","Var3","Var5","Var7","Var9","Var10","Var11")], collapse = '_' )})
tb['p2'] = apply(X = tmp, MARGIN = 1, FUN = function(x) { paste(x[c("Var2","Var4","Var6","Var8","Var9","Var10","Var11")], collapse = '_' )})
tb['comb'] = apply(X = tb[,c('p1','p2')], MARGIN = 1, FUN = function(x){ paste(x, collapse = ':')})
colnames(tb) = c('alone1','alone2','intra1fx','intra2fx','inter1fx','inter2fx','fit_alone1','fit_alone2','epis','loss','loss2','p1','p2','comb')

# precalculating conjugation, fitness and loss in function of the interaction
tb['intra1'] = tb$alone1 * tb$intra1fx; tb['intra2'] = tb$alone2 * tb$intra2fx # intracellular conjugation efficiency of plasmid 1 (X); and plasmid 2 (Y)
tb['inter1'] = tb$alone1 * tb$inter1fx; tb['inter2'] = tb$alone2 * tb$inter2fx # intercellular conjugation efficiency of plasmid 1 (X); and plasmid 2 (Y)
tb['fitness'] = (tb$fit_alone1 * tb$fit_alone2) + tb$epis # fitness of doubles XY
tb['cotransf'] = apply(X = tb[,c('intra1','intra2')], MARGIN = 1, FUN = min) # contransfer of X and Y from doubles
tb['exx1fx'] = ifelse(test = tb$intra1fx > 1, yes = tb$alone1, no = tb$intra1)  # entry exclusion effect for plasmid 1 (X)
tb['exx2fx'] = ifelse(test = tb$intra2fx > 1, yes = tb$alone2, no = tb$intra2)  # entry exclusion effect for plasmid 2 (Y)
tb['lossXY'] = tb$loss * tb$loss2 # combined loss
rownames(tb) = tb$comb; View(tb)


# steady state ----
# running the simulation in parallel/multiple cores for each row of conditions in the table tb generated above; the model and conditions are the same used above
registerDoParallel(detectCores()*3/4)
system.time({
  tab = foreach (i = 1:nrow(tb), .inorder = F) %do% { # saving the results of the simulation for steady state in the (dummy) list "tab"
    ini = c(O = 0, X = 1e6, Y =0, XY = 0, R = 100) # initial conditions: 1x10^6 cells X only (no Ø, Y nor XY cells) and initial resources = 100
    
    parms = c(Ro = 100, D = 1e-6, W = .05, Q = 5, G = 3 * log(2), # pre-set parameters:
              # resource concentration, resource required per cell division, turnover, Monod, growth rate
              
              s = tb[i,'loss'], s2 = tb[i,'lossXY'],# # segregation rates for each plasmid X and Y respectivelly
              hX = tb[i,'alone1'], hY = tb[i,'alone2'], # conjugation rates alone for each plasmid X and Y respectivelly
              hXy = tb[i,'intra1'], hYx = tb[i,'intra2'], # intracellular conjugation rates for each plasmid X and Y respectivelly
              hX_Y = tb[i,'inter1'], hY_X = tb[i,'inter2'], # intercellular conjugation rates for each plasmid X and Y respectivelly
              hXy2 = tb[i,'exx1fx'], hYx2 = tb[i,'exx2fx'], # exclusion on conjugation rates for each plasmid X and Y respectivelly
              hXY = tb[i,'cotransf'], # co-transfer rate
              FX = tb[i,'fit_alone1'], FY = tb[i,'fit_alone2'], FXY = tb[i,'fitness']) # fitness of cells X, Y and XY
    
    # defining the model
    model <- function (time, y, parms) {
      with(as.list(c(y, parms)), {
        
        dO = O*G*(R/(R+Q)) - W*O + s*(X+Y) - hX*O*X - hY*O*Y - (hXY + (hXy-hXY) + (hYx-hXY))*O*XY # Ø cells
        
        dX =  X*G*FX*(R/(R+Q)) - W*X + s2*XY - s*X + hX*O*X + (hXy-hXY)*O*XY - hY_X*X*Y - hYx2*X*XY # X cells
        
        dY =  Y*G*FY*(R/(R+Q)) - W*Y + s2*XY - s*Y + hY*O*Y + (hYx-hXY)*O*XY - hX_Y*X*Y - hXy2*Y*XY # Y cells
        
        dXY =  XY*G*FXY*(R/(R+Q)) - W*XY - 2*s2*XY + hXY*O*XY + (hX_Y+hY_X)*Y*X + hYx2*X*XY + hXy2*Y*XY # XY cells
        
        dR =  W*(Ro-R) - D*G*(R/(R+Q))* (O + X*FX + Y*FY + XY*FXY) # resources
        
        list(c(dO,dX,dY,dXY,dR))
      })
    }
    
    tmp = runsteady(y = ini, fun = model, parms = parms, times = c(0, 1e10), maxsteps = 1e9) # running the model for steady state
    
    # appeding precision, steady, time and number os steps required to reach teh steady state, combination
    cbind(data.frame(t(tmp$y),t(unlist(attributes(tmp)[c('precis','steady','time','steps')]))), data.frame(ID = rownames(tb[i,])))
  } 
})
stopImplicitCluster()

tab_cmplt = rbindlist(tab) # reshaping the results list into a single dataframe (tab_cmplt)
unique(tab_cmplt$steady) # if 1, all simulations reached the stead state
nrow(tab_cmplt[tab_cmplt$X > 1,]); View(tab_cmplt[tab_cmplt$X > 1,]) # cases where plasmid X is stably maintained
# checking maximum time, steps and precision of the siulations
max(as.numeric(as.character(tab_cmplt$time))); max(as.numeric(as.character(tab_cmplt$steps))); max(as.numeric(as.character(tab_cmplt$precis)))

write.table(x = tab_cmplt, file = 'Downloads/gama_steady_state_solo_density_REVIEW.txt') # saving the results of the steady state simulations as a tab separated file

# saving the table parameters used as a tab separated file
tab_cmplt$ID = as.character(tab_cmplt$ID)
tb['time'] = apply(X = as.matrix(tb$comb), MARGIN = 1, FUN = function(x) { ceiling(unlist(tab_cmplt[tab_cmplt$ID == x, 'time'])) } )
write.table(x = tb, file = 'Downloads/gama_parameters_solo_density_REVIEW.txt')




# timewise simulation for extinctions ----
# subsetting the parameter table to keep only the cases of extinction to be able to calculate their time to extinction
tb0 = tb
tb = tb[tb$comb %in% unlist(tab_cmplt[tab_cmplt$X < 1 & tab_cmplt$XY < 1, "ID"]), ]
fivenum(tb$time) # distribution of times required for steady state

# running the model to get times to extinction: the model and conditions are the same used above
# running the simulation in parallel/multiple cores for each row of conditions in the table tb generated above
registerDoParallel(detectCores()*3/4)
system.time({
  tab2 = # saving the results of the simulation for times of extinction in the (dummy) list "tab2" 
    foreach (i = 1:nrow(tb), .inorder = F) %do% {
      ini = c(O = 0, X = 1e6, Y = 0, XY = 0, R = 100) # initial conditions: 1x10^6 cells X only (no Ø, Y nor XY cells) and initial resources = 100
      
      parms = c(Ro = 100, D = 1e-6, W = .05, Q = 5, G = 3 * log(2), # pre-set parameters:
                # resource concentration, resource required per cell division, turnover, Monod, growth rate
                
                s = tb[i,'loss'], s2 = tb[i,'lossXY'], # segregation rates for each plasmid X and Y respectivelly
                hX = tb[i,'alone1'], hY = tb[i,'alone2'], # conjugation rates alone for each plasmid X and Y respectivelly
                hXy = tb[i,'intra1'], hYx = tb[i,'intra2'], # intracellular conjugation rates for each plasmid X and Y respectivelly
                hX_Y = tb[i,'inter1'], hY_X = tb[i,'inter2'], # intercellular conjugation rates for each plasmid X and Y respectivelly
                hXy2 = tb[i,'exx1fx'], hYx2 = tb[i,'exx2fx'], # exclusion on conjugation rates for each plasmid X and Y respectivelly
                hXY = tb[i,'cotransf'], # co-transfer rate
                FX = tb[i,'fit_alone1'], FY = tb[i,'fit_alone2'], FXY = tb[i,'fitness']) # fitness of cells X, Y and XY
      
      times <- seq(0, tb[i,'time'], 1) # cycling times hourly from 0 up to the time of steady state specific for each combination
      
      # defining the model
      model <- function (time, y, parms) {
        with(as.list(c(y, parms)), {
          
          dO = O*G*(R/(R+Q)) - W*O + s*(X+Y) - hX*O*X - hY*O*Y - (hXY + (hXy-hXY) + (hYx-hXY))*O*XY # Ø cells
          
          dX =  X*G*FX*(R/(R+Q)) - W*X + s2*XY - s*X + hX*O*X + (hXy-hXY)*O*XY - hY_X*X*Y - hYx2*X*XY # X cells
          
          dY =  Y*G*FY*(R/(R+Q)) - W*Y + s2*XY - s*Y + hY*O*Y + (hYx-hXY)*O*XY - hX_Y*X*Y - hXy2*Y*XY # Y cells 
          
          dXY =  XY*G*FXY*(R/(R+Q)) - W*XY - 2*s2*XY + hXY*O*XY + (hX_Y+hY_X)*Y*X + hYx2*X*XY + hXy2*Y*XY # XY cells
          
          dR =  W*(Ro-R) - D*G*(R/(R+Q))* (O + X*FX + Y*FY + XY*FXY) # resources
          
          list(c(dO,dX,dY,dXY,dR))
        })
      }
      
      tmp = data.frame(ode(ini, times, model, parms)) # running the model for the times specified above (in the vector "times")
      
      # keeping the results as the time of extinction for X, and the ID of the parameters combination
      # time of extinction is the first time where the entity amount is < 1
      cbind(data.frame(X = min(tmp[tmp$X < 1, 'time']), ID = rownames(tb[i,])))
    } 
})
stopImplicitCluster()

tab_cmplt2 = rbindlist(tab2) # reshaping the results list into a single dataframe (tab_cmplt2)

# checking if there is any non natural value in the final table (ie, any error)
any(is.na(tab_cmplt2$X))
any(is.infinite(tab_cmplt2$X))

write.table(x = tab_cmplt2, file = 'Downloads/gama_time_to_extinction_solo_density_REVIEW.txt') # saving the results of the time of extinction simulations as a tab separated file

