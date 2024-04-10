#Modelling with Sicegar####
# Data tracking
require("sicegar")
require("tidyverse")

# Graphing
require("ggplot2")
require("cowplot")

# Text manipulation
require("stringr")
require("rlist")
require("rstatix")
library("eply")

#Testing model input parameters ####
#models a single selected time series functionally. 
#useful for testing changes to categorize() parameters
temp <- plot_data %>% 
  dplyr:::filter(unique_ID=="20230829_10_3_10") %>% 
  ungroup() %>% 
  dplyr:::select(HPI, normalised_sum_red) %>% 
  dplyr:::rename("time"=HPI, "intensity"=normalised_sum_red)

temp <- temp_mean %>% 
  ungroup() %>% 
  dplyr:::select(HPI, meanGreen) %>% 
  dplyr:::rename("time"=HPI, "intensity"=meanGreen)


temp1<- sicegar::normalizeData(temp)

sigmoidalModel <- sicegar::multipleFitFunction(dataInput = temp1,
                                               model = "sigmoidal",
                                               n_runs_min = 20,
                                               n_runs_max = 500,
                                               showDetails = FALSE)
sigmoidalModel <- sicegar::parameterCalculation(sigmoidalModel)

doubleSigmoidalModel <- multipleFitFunction(dataInput = temp1,
                                            model = "doublesigmoidal",
                                            n_runs_min = 20,
                                            n_runs_max = 500,
                                            showDetails = FALSE)
doubleSigmoidalModel <- parameterCalculation(doubleSigmoidalModel)

doubleSigmoidalModel$start_intensity

fig01 <- sicegar::figureModelCurves(dataInput = temp1,
                                    sigmoidalFitVector = sigmoidalModel,
                                    showParameterRelatedLines = FALSE)
fig02 <- figureModelCurves(dataInput = temp1,
                           doubleSigmoidalFitVector = doubleSigmoidalModel,
                           showParameterRelatedLines = FALSE)
  


plot_grid(fig01,fig02)

decisionProcess <- categorize(threshold_minimum_for_intensity_maximum = 10,
                              threshold_intensity_range = 4,
                              threshold_t0_max_int = 100,                  #always need to assess raw data to determine this param
                              threshold_sm_tmax_IntensityRatio = 0.65,
                              threshold_dsm_tmax_IntensityRatio = 0.65,
                              threshold_AIC = -10,
                              parameterVectorSigmoidal = sigmoidalModel,
                              parameterVectorDoubleSigmoidal = doubleSigmoidalModel) 

#Batch curve fitting ####
seedNo=14159   #maintains the same random seed for iterative curve fitting
set.seed(seedNo)

min_max_val <- 10 #set the threshold_minimum_for_intensity_maximum parameter for categorize(). Needs to be considered on a per dataset basis

data_list_red<- plot_data %>%  
  filter(HPI<=24) %>% 
  group_by(unique_ID) %>%  
  dplyr:::group_split()%>% #Split df into a list of dfs, one per cell
  stats:::setNames(sort(unique(plot_data$unique_ID)))

data_list_green <- plot_data %>%
  filter(HPI<=24) %>% 
  group_by(unique_ID) %>%  
  dplyr:::group_split()%>% #Split df into a list of dfs, one per cell
  stats:::setNames(sort(unique(plot_data$unique_ID)))

#mCh model fitting ####
sicegar_fitFunction_red<-function(data)     #if threshold_minimum_for_intensity_maximum is set too low this can return an error $model argument is of length 0, or similar
 {
    dataInput1<-data.frame(time=data$HPI, intensity=data$normalised_sum_red)
    dataInputName<-unique(data$unique_ID)
    dataInputCondition<-unique(data$Condition)
    #print(dataInputName)
    fitObj <- sicegar::fitAndCategorize(dataInput = dataInput1,
                              dataInputName = dataInputName,
                              threshold_minimum_for_intensity_maximum = min_max_val,
                              threshold_intensity_range = 4,
                              threshold_t0_max_int = 100,                  #always need to assess raw data to determine this param
                              threshold_sm_tmax_IntensityRatio = 0.65,
                              threshold_dsm_tmax_IntensityRatio = 0.65,
                              threshold_AIC = -10,
                              parameterVectorSigmoidal = sigmoidalModel,
                              parameterVectorDoubleSigmoidal = doubleSigmoidalModel)
    
    # Data Scaling Parameters
    dataScalingParametersDF <- data.frame(t(fitObj$normalizedInput$dataScalingParameters))
    colnames(dataScalingParametersDF) <- paste0("SCALING_", colnames(dataScalingParametersDF))
    
    # Sigmoidal
    if(typeof(fitObj$sigmoidalModel)!="logical")
    {
      sigmoidalDF <- data.frame(fitObj$sigmoidalModel)
      sigmoidalDF[grep(pattern = "dataScalingParameters", colnames(sigmoidalDF))]<-NULL
      colnames(sigmoidalDF) <- paste0("SM_", colnames(sigmoidalDF))
    }
    
    # Double Sigmoidal
    if(typeof(fitObj$doubleSigmoidalModel)!="logical")
    {
      doublesigmoidalDF <- data.frame(fitObj$doubleSigmoidalModel)
      doublesigmoidalDF[grep(pattern = "dataScalingParameters", colnames(doublesigmoidalDF))]<-NULL
      colnames(doublesigmoidalDF) <- paste0("DSM_", colnames(doublesigmoidalDF))
    }
    
    #Decision Process
    if(typeof(fitObj$decisionProcess)!="logical")
    {
      decisionProcessDF <- as.data.frame(t(fitObj$decisionProcess))
      colnames(decisionProcessDF) <- paste0("DEC_", colnames(decisionProcessDF))
    }
    
    # Summary File
    summaryDF <- data.frame(t(fitObj$summaryVector))
    summaryDF[,c(1,2)] <- NULL
    if(decisionProcessDF["DEC_decision"]=="sigmoidal")
    {
      summaryDF %>%
        dplyr::rename(midPoint1_x = midPoint_x,
                      midPoint1_y = midPoint_y,
                      slope1 = slope)->summaryDF
    }
    
    if(ncol(summaryDF)!=0)
    {
      colnames(summaryDF) <- paste0("COMB_", colnames(summaryDF)) 
    }
    
    if(exists("doublesigmoidalDF") & exists("sigmoidalDF"))
    {
      output <- dplyr::bind_cols(dataScalingParametersDF, 
                                 sigmoidalDF, 
                                 doublesigmoidalDF, 
                                 decisionProcessDF,
                                 summaryDF)
    }
    else
    {
      output <- dplyr::bind_cols(dataScalingParametersDF, 
                                 decisionProcessDF,
                                 summaryDF)
    }
    
    names(output) <- sub("\\.", "_",names(output))
    output2 <- purrr::flatten(as.vector(output))
    names(output2) <- names(output)
    output2$unique_ID <- dataInputName  
    output2$Condition <- dataInputCondition
    return(output2)
}

sicegar_fitFunction_flexible_red <- possibly(.f = sicegar_fitFunction_red, otherwise = NULL) #sicegar:::parameterCalculation throws out random length zero errors when evaluating. 

 #loops the sicegar_fitFunction over the list of dataframes
sicegar_results_red<- map(data_list_red, sicegar_fitFunction_flexible_red) 

sicegar_results_filter_red<- sicegar_results_red %>% 
  compact()   #removes NULL elements from the list

data_filter<- function(.data_list) #a set of criteria for filtering out relevant elements of sicegar_results
{ 
  as.data.frame(.data_list) %>% 
  dplyr:::select(contains("COMB"), unique_ID, Condition, DEC_decision)
}

temp_red<- map(sicegar_results_filter_red, data_filter)

sicegar_results_df_red <- do.call("rbind.fill", temp_red) #binds together the filtered sicegar_results, filling any missing columns with NAs

#write.csv(sicegar_results_df_red, "sicegar_results_red_MOI.csv")

model_frequencies_red <- sicegar_results_df_red %>% 
  group_by(DEC_decision) %>% 
  dplyr::summarise(count=n())

#merge sicegar results for both pEL-mCh and A3-GFP
sicegar_results_df_red_temp<-sicegar_results_df_red %>% 
  setNames(paste0(names(.),'_red')) %>%
  dplyr:::rename("unique_ID"=unique_ID_red)
  
#YFP model fitting ####

min_max_val <- 20 #set the threshold_minimum_for_intensity_maximum parameter for categorize(). Needs to be considered on a per dataset basis

sicegar_fitFunction_green<-function(data)     #if threshold_minimum_for_intensity_maximum is set too low this can return an error $model argument is of length 0, or similar
{
  dataInput1<-data.frame(time=data$HPI, intensity=data$normalised_sum_green)
  dataInputName<-unique(data$unique_ID)
  dataInputCondition<-unique(data$Condition)
  #print(dataInputName)
  fitObj <- sicegar::fitAndCategorize(dataInput = dataInput1,
                                      dataInputName = dataInputName,
                                      threshold_minimum_for_intensity_maximum = min_max_val,
                                      threshold_intensity_range = 4,
                                      threshold_t0_max_int = 30,                  #always need to assess raw data to determine this param
                                      threshold_sm_tmax_IntensityRatio = 0.65,
                                      threshold_dsm_tmax_IntensityRatio = 0.65,
                                      threshold_AIC = -10,
                                      parameterVectorSigmoidal = sigmoidalModel,
                                      parameterVectorDoubleSigmoidal = doubleSigmoidalModel)
  
  # Data Scaling Parameters
  dataScalingParametersDF <- data.frame(t(fitObj$normalizedInput$dataScalingParameters))
  colnames(dataScalingParametersDF) <- paste0("SCALING_", colnames(dataScalingParametersDF))
  
  # Sigmoidal
  if(typeof(fitObj$sigmoidalModel)!="logical")
  {
    sigmoidalDF <- data.frame(fitObj$sigmoidalModel)
    sigmoidalDF[grep(pattern = "dataScalingParameters", colnames(sigmoidalDF))]<-NULL
    colnames(sigmoidalDF) <- paste0("SM_", colnames(sigmoidalDF))
  }
  
  # Double Sigmoidal
  if(typeof(fitObj$doubleSigmoidalModel)!="logical")
  {
    doublesigmoidalDF <- data.frame(fitObj$doubleSigmoidalModel)
    doublesigmoidalDF[grep(pattern = "dataScalingParameters", colnames(doublesigmoidalDF))]<-NULL
    colnames(doublesigmoidalDF) <- paste0("DSM_", colnames(doublesigmoidalDF))
  }
  
  #Decision Process
  if(typeof(fitObj$decisionProcess)!="logical")
  {
    decisionProcessDF <- as.data.frame(t(fitObj$decisionProcess))
    colnames(decisionProcessDF) <- paste0("DEC_", colnames(decisionProcessDF))
  }
  
  # Summary File
  summaryDF <- data.frame(t(fitObj$summaryVector))
  summaryDF[,c(1,2)] <- NULL
  if(decisionProcessDF["DEC_decision"]=="sigmoidal")
  {
    summaryDF %>%
      dplyr::rename(midPoint1_x = midPoint_x,
                    midPoint1_y = midPoint_y,
                    slope1 = slope)->summaryDF
  }
  
  if(ncol(summaryDF)!=0)
  {
    colnames(summaryDF) <- paste0("COMB_", colnames(summaryDF)) 
  }
  
  if(exists("doublesigmoidalDF") & exists("sigmoidalDF"))
  {
    output <- dplyr::bind_cols(dataScalingParametersDF, 
                               sigmoidalDF, 
                               doublesigmoidalDF, 
                               decisionProcessDF,
                               summaryDF)
  }
  else
  {
    output <- dplyr::bind_cols(dataScalingParametersDF, 
                               decisionProcessDF,
                               summaryDF)
  }
  
  names(output) <- sub("\\.", "_",names(output))
  output2 <- purrr::flatten(as.vector(output))
  names(output2) <- names(output)
  output2$unique_ID <- dataInputName  
  output2$Condition <- dataInputCondition
  return(output2)
}

sicegar_fitFunction_flexible_green <- possibly(.f = sicegar_fitFunction_green, otherwise = NULL) #sicegar:::parameterCalculation throws out random length zero errors when evaluating. 

#loops the sicegar_fitFunction over the list of dataframes
sicegar_results_green<- map(data_list_green, sicegar_fitFunction_flexible_green) 

sicegar_results_filter_green<- sicegar_results_green %>% 
  compact()   #removes NULL elements from the list

data_filter<- function(.data_list) #a set of criteria for filtering out relevant elements of sicegar_results
{ 
  as.data.frame(.data_list) %>% 
    dplyr:::select(contains("COMB"), unique_ID, Condition, DEC_decision)
}

#try making a function that pulls out the NULL results and reruns them through sicegar_fitFunction_flexible

temp_green<- map(sicegar_results_filter_green, data_filter)

sicegar_results_df_green <- do.call("rbind.fill", temp_green) #binds together the filtered sicegar_results, filling any missing columns with NAs

#write.csv(sicegar_results_df_green, "sicegar_results_green_MOI.csv")

model_frequencies_green <- sicegar_results_df_green %>% 
  group_by(DEC_decision) %>% 
  dplyr::summarise(count=n())

sicegar_results_df_merged_24HPI <- sicegar_results_df_red_temp %>% 
  left_join(., sicegar_results_df_green, by="unique_ID")

write.csv(sicegar_results_df_merged_24HPI, "sicegar_results_MOI_24HPI.csv")

plot_data_with_sicegar_2_temp<- plot_data %>% 
  left_join(., sicegar_results_df_merged, by="unique_ID")
  
plot_data_with_sicegar_2_temp<-plot_data_with_sicegar_2_temp %>% 
  group_by(Condition_red) %>% 
  dplyr::mutate(COMB_startPoint_x_red2=(start_time*10)/60,
                COMB_startPoint_x2=(start_time_green*10)/60) %>% 
  group_by(unique_ID) %>% 
  filter(((start_time - 1) > min(Metadata_Frame, na.rm=TRUE) | min(Metadata_Frame, na.rm=TRUE) == 0 | min(Metadata_Frame, na.rm=TRUE) < 20)) 
  

model_frequencies <- sicegar_results_df_merged %>% 
  filter(DEC_decision_red=="double_sigmoidal"|DEC_decision_red=="sigmoidal") %>% 
  group_by(DEC_decision, DEC_decision_red, Condition_red) %>% 
  dplyr::summarise(count=n())

plot_data <- read.csv("plot_data.csv")
sicegar_results_df_merged <- read.csv("sicegar_results_MOI.csv")

write.csv(plot_data_with_sicegar_2_temp, "plot_data_with_sicegar.csv")
