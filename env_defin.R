library(bigreadr)
library(dplyr)
library(plyr)

# Set parameters
PROJ_PATH      <- "/public/home/gw_hychu/swj/03_trait/proteohubProject"
TRAIT          <- "/public/home/gw_hychu/swj/03_trait/participant.csv.gz"

# Load sample id
all_id <- fread2(paste0(PROJ_PATH, "/02_data/sample_id/all_eid"))[, 1]
male_id <- fread2(paste0(PROJ_PATH, "/02_data/sample_id/male_eid"))[, 1]
female_id <- fread2(paste0(PROJ_PATH, "/02_data/sample_id/female_eid"))[, 1]
id_name <- c("all", "female", "male")
trait_exposure_dat <- trait_exposure_dat[match(all_id, trait_exposure_dat$eid), ]
trait_exposure_dat$sex <- ifelse(all_id%in%male_id, "male", "female") 


# Income, Education and Employment groups #
## Categorize income groups
incm_cnd <- recode(trait_exposure_dat[["p738_i0"]],
                   "Less than 18,000" = "1",
                   "18,000 to 30,999" = "2",
                   "31,000 to 51,999" = "3",
                   "52,000 to 100,000" = "4",
                   "Greater than 100,000" = "5") %>% as.integer()
incm_out <- recode(trait_exposure_dat[["p738_i0"]],
                   "Less than 18,000" = "<18,000",
                   "18,000 to 30,999" = "18,000~30,999",
                   "31,000 to 51,999" = "31,000~51,999",
                   "52,000 to 100,000" = "52,000~100,000",
                   "Greater than 100,000" = ">1000,000", 
                   "Do not know" = "NA", 
                   "Prefer not to answer" = "NA") 
incm_out[which(incm_out == "")] <- "NA"
incm_out[which(incm_out == "NA")] <- NA
incm_out <-  factor(incm_out, levels = c("<18,000", "18,000~30,999", 
                                         "31,000~51,999", "52,000~100,000", 
                                         ">1000,000"))

## Categorize education groups (p6138_i0)
### set 6138_i0 as the highest education status
### 1-2: College or above
### 3-6: High school or equivalent
### -7: Less than high school (33853828)
edu_list <- str_split_i(trait_exposure_dat[["p6138_i0"]], "\\|", 1)
edu_cnd <- recode(edu_list,
                  "College or University degree" = "3",
                  "A levels/AS levels or equivalent" = "3",
                  "O levels/GCSEs or equivalent" = "2",
                  "CSEs or equivalent" = "2",
                  "NVQ or HND or HNC or equivalent" = "2",
                  "Other professional qualifications eg: nursing, teaching" = "2",
                  "None of the above" = "1") %>% as.integer()
edu_out <- recode(edu_list,
                  "College or University degree" = "College or above",
                  "A levels/AS levels or equivalent" = "College or above",
                  "O levels/GCSEs or equivalent" = "High school or equivalent",
                  "CSEs or equivalent" = "High school or equivalent",
                  "NVQ or HND or HNC or equivalent" = "High school or equivalent",
                  "Other professional qualifications eg: nursing, teaching" = "High school or equivalent",
                  "None of the above" = "Less than high school",
                  "Prefer not to answer" = "NA") 
edu_out[which(edu_out == "")] <- "NA"
edu_out[which(edu_out == "NA")] <- NA
edu_out <-  factor(edu_out, levels = c("College or above", 
                                       "High school or equivalent", 
                                       "Less than high school"))

################## Employment ##################
# Employed Status including those in paid employment or self-employed, retired, 
# doing unpaid or voluntary work, or being full or part time students)
emp_list <- str_split_i(trait_exposure_dat[["p6142_i0"]], "\\|", 1)
emp_cnd <- recode(emp_list,
                  "In paid employment or self-employed" = "2",
                  "Retired" = "2",
                  "Doing unpaid or voluntary work" = "2",
                  "Full or part-time student" = "2",
                  "Looking after home and/or family" = "1",
                  "Unable to work because of sickness or disability" = "1",
                  "None of the above" = "NA") %>% as.integer()
emp_out <- recode(emp_list,
                  "In paid employment or self-employed" = "Employed",
                  "Retired" = "Employed",
                  "Doing unpaid or voluntary work" = "Employed",
                  "Full or part-time student" = "Employed",
                  "Unemployed" = "Unemployed",
                  "Looking after home and/or family" = "Unemployed",
                  "Unable to work because of sickness or disability" = "Unemployed",
                  "None of the above" = "Unemployed", 
                  "Prefer not to answer" = "NA")
emp_out[which(emp_out == "")] <- "NA"
emp_out[which(emp_out == "NA")] <- NA
emp_out <-  factor(emp_out, levels = c("Employed", "Unemployed"))

################## Estimate SES ##################
ses_df <- data.frame(Income = incm_cnd, 
                     Education = edu_cnd, 
                     Employment = emp_cnd,
                     row.names = trait_exposure_dat$eid)
SES <- rep(NA, nrow(ses_df))
nna_idx <- which(rowSums(is.na(ses_df)) == 0)
set.seed(20250207)
SES_LCA3 <- poLCA(cbind(Income, Education, Employment) ~ 1, 
                  data = ses_df[nna_idx,],
                  nclass = 3, 
                  maxiter = 10000, 
                  graphs = F, 
                  tol = 1e-6, 
                  na.rm = T, 
                  probs.start = NULL, 
                  nrep = 1, 
                  verbose = T, 
                  calc.se = T)
SES[nna_idx] <- SES_LCA3$predclass
high_idx <- which.max(SES_LCA3$probs$Income[,5] + SES_LCA3$probs$Income[,4])
low_idx <- which.max(SES_LCA3$probs$Income[,1])
medium_idx <- setdiff(1:3, c(high_idx, low_idx))
SES_out <- factor(SES, levels = c(low_idx, medium_idx, high_idx), 
                  labels = c("Low", "Medium", "High"))


################## no current smoke ##################
# "p1239_i0", 	 ## 02.01 Current tobacco smoking
# "p1249_i0", 	 ## 02.02 Past tobacco smoking
# "p2644_i0", 	 ## 02.03 Light smokers, at least 100 smokes in lifetime
crt_smk_cnd <- recode(trait_exposure_dat[["p1239_i0"]],
                      "Yes, on most or all days" = "1",
                      "Only occasionally" = "0",
                      "No" = "0") %>% as.integer()
pst_smk_cnd <- recode(trait_exposure_dat[["p1249_i0"]],
                      "Smoked on most or all days" = "1",
                      "Smoked occasionally" = "0",
                      "Just tried once or twice" = "0",
                      "I have never smoked" = "0") %>% as.integer()
npst_smk_cnd <- ifelse(trait_exposure_dat[["p2644_i0"]] == "Yes", 1, 0)
nosmoke_cnd <- rep(NA, nrow(trait_exposure_dat))
nosmoke_cnd[crt_smk_cnd | pst_smk_cnd | npst_smk_cnd] <- 0
nosmoke_cnd[!crt_smk_cnd & !pst_smk_cnd & !npst_smk_cnd] <- 1
nosmoke_out <- ifelse(nosmoke_cnd == 0, "No", "Yes") %>% 
  factor(., levels = c("No", "Yes")) 
table(nosmoke_out)

################### activity ##################
# calculate metabolic equivalent times for each participant by adding 
# time spent on each activity weighted by its metabolic equivalent score.
# metabolic equivalent scores:
# walking for pleasure: 3.3
# strenuous sports: 8
# other exercises: 4
# "p971_i0", 	   ## 03.01.01 Frequency of walking for pleasure in last 4 weeks
# "p991_i0", 	   ## 03.01.02 Frequency of strenuous sports in last 4 weeks
# "p3637_i0", 	 ## 03.01.03 Frequency of other exercises in last 4 weeks
# "p981_i0", 	   ## 03.02.01 Duration walking for pleasure
# "p1001_i0", 	 ## 03.02.02 Duration of strenuous sports
# "p3647_i0", 	 ## 03.02.03 Duration of other exercises
act_metab_weight <- c(3.3, 8, 4)
act_freq_mat <- lapply(c("p971_i0", "p991_i0", "p3637_i0"), function(x){
  
  recode(trait_exposure_dat[[x]],
         "Once in the last 4 weeks" = "0.25",
         "2-3 times in the last 4 weeks" = "0.625",
         "Once a week" = "1",
         "2-3 times a week" = "2.5",
         "4-5 times a week" = "4.5") %>% as.numeric()
}) %>% Reduce("cbind", .)
act_dura_mat <- lapply(c("p981_i0", "p1001_i0", "p3647_i0"), function(x){
  
  recode(trait_exposure_dat[[x]],
         "Less than 15 minutes" = "0.25",
         "Between 15 and 30 minutes" = "0.375",
         "Between 30 minutes and 1 hour" = "0.75",
         "Between 1 and 1.5 hours" = "1.25",
         "Between 1.5 and 2 hours" = "1.75",
         "Between 2 and 3 hours" = "2.5",
         "Over 3 hours" = "3") %>% as.numeric()
}) %>% Reduce("cbind", .)
act_freq_mat[is.na(act_freq_mat)] <- 0
act_dura_mat[is.na(act_dura_mat)] <- 0
act_metab_score <- (act_freq_mat[,1] * act_dura_mat[,1] * act_metab_weight[1]) +
  (act_freq_mat[,2] * act_dura_mat[,2] * act_metab_weight[2]) + 
  (act_freq_mat[,3] * act_dura_mat[,3] * act_metab_weight[3])
act_metab_score[act_metab_score == 0] <- NA
act_cnd <- act_metab_score > quantile(act_metab_score, 1/3, na.rm = T)
act_out <- ifelse(act_cnd, "High", "Low") %>%
  factor(., levels = c("Low", "High"))

################## healthy diet ##################
### Fruits: 3 servings/day
# "p1309_i0", 	 ## 04.01.01 Fresh fruit intake per DAY
# "p1319_i0", 	 ## 04.01.02 Dried fruit intake per DAY
# Amount per serving: 1309 ?C 1 piece OR 1319 ?C 5 pieces 
fruit_mat <- cbind(as.numeric(trait_exposure_dat[["p1309_i0"]]), 
                   as.numeric(trait_exposure_dat[["p1319_i0"]])/5) 
fruit_cnd <- rowSums(fruit_mat) >= 3
fruit_cnd[rowSums(fruit_mat, na.rm = T) >= 3] <- T
### vegetable: 3 servings/day
# "p1289_i0", 	 ## 04.02.01 Cooked vegetable intake per DAY
# "p1299_i0", 	 ## 04.02.02 Salad / raw vegetable intake per DAY
# Amount per serving: 3 heaped tablespoons 
veg_mat <- cbind(as.numeric(trait_exposure_dat[["p1289_i0"]])/3, 
                 as.numeric(trait_exposure_dat[["p1299_i0"]])/3)
veg_cnd <- rowSums(veg_mat) >= 3
veg_cnd[rowSums(veg_mat,na.rm = T) >= 3] <- T
### fish: 2 servings/week
# "p1329_i0", 	 ## 04.03.01 Oily fish intake
# "p1339_i0", 	 ## 04.03.02 Non-oily fish intake
# Amount per serving: Once/week 
fish_mat <- lapply(c("p1329_i0", "p1339_i0"), function(x){
  recode(trait_exposure_dat[[x]],
         "Never" = "0", 
         "Less than once a week" = "0.5", 
         "Once a week" = "1",
         "2-4 times a week" = "3", 
         "5-6 times a week" = "5.5", 
         "Once or more daily" = "7") %>% as.integer()
}) %>% Reduce("cbind", .)
fish_cnd <- rowSums(fish_mat) >= 2
fish_cnd[rowSums(fish_mat,na.rm = T) >= 2] <- T
### Diary: 2 servings/day 
# "p1408_i0",    ## 04.04.01 Cheese intake
# "p1418_i0",    ## 04.04.02 Milk type used
# Amount per serving: 1408 ?C 1 piece/day OR
# 1418 ?C 1 glass/day if consumption of any type of milk 
chz_intake <- ifelse(trait_exposure_dat[["p1408_i0"]] %in% c("", "Do not know", "Prefer not to answer"), NA,
                     ifelse(trait_exposure_dat[["p1408_i0"]] == "Once or more daily", T, F))
milk_intake <- ifelse(trait_exposure_dat[["p1408_i0"]] %in% c("", "Do not know", "Prefer not to answer"), NA,
                      ifelse(trait_exposure_dat[["p1408_i0"]] == "Never/rarely have milk", F, T))
diary_cnd <- chz_intake & milk_intake
### Vegetable oils: 2 servings/day 
# "p2654_i0", 	 ## 04.04.02 Non-butter spread type details
# "p1438_i0", 	 ## 04.04_05 Bread intake each WEEK
# Amount per serving: 1 serving/day if in combination with 
# eating at least TWO slices of bread (ID 1438) 
veg_oil_na_mat <- lapply(c("p1438_i0", "p2654_i0"), function(x){
  trait_exposure_dat[[x]] %in% c("", "Do not know", "Prefer not to answer") |
    is.na(trait_exposure_dat[[x]])
}) %>% Reduce("cbind", .)
veg_oil_intake <- rep(0, nrow(trait_exposure_dat))
veg_oil_idx <- which(trait_exposure_dat[["p2654_i0"]] %in% 
                       c("Flora Pro-Active/Benecol", 
                         "Soft (tub) margarine", 
                         "Olive oil based spread (eg: Bertolli)",
                         "Polyunsaturated/sunflower oil based spread (eg: Flora)",
                         "Other low or reduced fat spread"))
veg_oil_intake[veg_oil_idx] <- trait_exposure_dat[["p1438_i0"]][veg_oil_idx] %>% as.numeric()
veg_oil_cnd <- veg_oil_intake/7 >= 2
veg_oil_cnd[rowMeans(veg_oil_na_mat) == 1] <- NA
### Whole grains: 3 servings/day
# "p1438_i0", 	 ## 04.04_05 Bread intake each WEEK
# "p1448_i0", 	 ## 04.05 Bread type
# "p1458_i0", 	 ## 04.06.01 Cereal intake each WEEK
# "p1468_i0", 	 ## 04.06.02 Cereal type
# Amount per serving: 1438/1448 ?C 1 slice/day OR 1458/1468 ?C 1 bowl/day 
# treat "Less than one" as 0
grain_na_mat <- lapply(c("p1438_i0", "p1448_i0", "p1458_i0", "p1468_i0"), function(x){
  trait_exposure_dat[[x]] %in% c("", "Do not know", "Prefer not to answer") |
    is.na(trait_exposure_dat[[x]])
}) %>% Reduce("cbind", .)
whole_grain_mat <- matrix(0, nrow(trait_exposure_dat), 2)
# wholemeal/wholegrain bread
whole_grain_idx <- which(trait_exposure_dat[["p1448_i0"]] == "Wholemeal or wholegrain")
whole_grain_mat[whole_grain_idx, 1] <- trait_exposure_dat[["p1438_i0"]][whole_grain_idx]
# bran/oat/muesli cereal 
cereal_idx <- which(trait_exposure_dat[["p1468_i0"]] %in% c("Bran cereal (e.g. All Bran, Branflakes)", 
                                                            "Oat cereal (e.g. Ready Brek, porridge)", 
                                                            "Muesli"))
whole_grain_mat[cereal_idx, 2] <- trait_exposure_dat[["p1458_i0"]][cereal_idx]
whole_grain_mat[whole_grain_mat == "Less than one"] <- "0.5"
whole_grain_mat <- apply(whole_grain_mat, 2, function(x){as.numeric(x)})
whole_grain_cnd <- rowSums(whole_grain_mat, na.rm = T)/7 >= 3
whole_grain_cnd[rowMeans(grain_na_mat) == 1] <- NA
### Refine grains: 2 servings/day
# "p1438_i0", 	 ## 04.04_05 Bread intake each WEEK
# "p1448_i0", 	 ## 04.05 Bread type
# "p1458_i0", 	 ## 04.06.01 Cereal intake each WEEK
# "p1468_i0", 	 ## 04.06.02 Cereal type
# Amount per serving: 1438/1448 ?C 1 slice/day OR 1458/1468 ?C 1 bowl/day 
refine_grain_mat <- matrix(0, nrow(trait_exposure_dat), 2)
# wholemeal/wholegrain bread
refine_grain_idx <- which(trait_exposure_dat[["p1448_i0"]] %in% c("White", "Brown", "Other type of bread"))
refine_grain_mat[refine_grain_idx, 1] <- trait_exposure_dat[["p1438_i0"]][refine_grain_idx] %>% as.numeric()
# bran/oat/muesli cereal 
refine_cereal_idx <- which(trait_exposure_dat[["p1468_i0"]] %in% c("Biscuit cereal (e.g. Weetabix)", 
                                                                   "Other (e.g. Cornflakes, Frosties)"))
refine_grain_mat[refine_cereal_idx, 2] <- trait_exposure_dat[["p1458_i0"]][refine_cereal_idx] %>% as.numeric()
refine_grain_mat[refine_grain_mat == "Less than one"] <- "0.5"
refine_grain_mat <- apply(refine_grain_mat, 2, function(x){as.numeric(x)})
refine_grain_cnd <- rowSums(refine_grain_mat, na.rm = T)/7 <= 2
refine_grain_cnd[rowMeans(grain_na_mat) == 1] <- NA
### processed red meat: 1 times/week
# "p1349_i0", 	 ## 04.07.01 Processed meat intake
# "p3680_i0",    ## 04.08_09 Processed meat intake
# Amount per serving: 1349 ?C 1 piece/day
pmeat_cnd <- ifelse(trait_exposure_dat[["p1349_i0"]] %in% 
                      c("Never", "Less than once a week", "Once a week"), T,
                    ifelse(trait_exposure_dat[["p1349_i0"]] %in% 
                             c("2-4 times a week", "5-6 times a week", "Once or more daily"), F, NA))
pmeat_cnd[which(trait_exposure_dat[["p3680_i0"]] == "0")] <- T
### non-processed red meat: ?? 2 times/week
# "p1359_i0", 	 ## 04.08.01 Poultry intake
# "p1369_i0", 	 ## 04.08.02 Beef intake
# "p1379_i0", 	 ## 04.08.03 Lamb/mutton intake
# "p1389_i0", 	 ## 04.08.04 Pork intake
# "p3680_i0",    ## 04.08_09 Processed meat intake
# Amount per serving: 1359-1389 ?C once/week 
npmeat_mat <- lapply(c("p1359_i0", "p1369_i0", "p1379_i0", "p1389_i0"), function(x){
  recode(trait_exposure_dat[[x]],
         "Never" = "0", 
         "Less than once a week" = "0.5", 
         "Once a week" = "1",
         "2-4 times a week" = "3", 
         "5-6 times a week" = "5.5", 
         "Once or more daily" = "7") %>% as.integer()
}) %>% Reduce("cbind", .)
npmeat_cnd <- rowSums(npmeat_mat) <= 2
npmeat_cnd[rowSums(npmeat_mat, na.rm = T) > 2] <- F
npmeat_cnd[which(trait_exposure_dat[["p3680_i0"]] == "0")] <- T
### Sugar-sweetened beverages : No consumption 
# "p6144_i0",    ## Never eat eggs, dairy, wheat, sugar
sugar_cnd <- rep(F, nrow(trait_exposure_dat))
sugar_cnd[grep("Sugar or foods/drinks containing sugar", trait_exposure_dat[["p6144_i0"]])] <- T
sugar_cnd[sugar_cnd %in% c("", "Do not know", "Prefer not to answer")] <- NA
### Healthy diet: At least 5 of the 10 food groups
diet_mat <- data.frame(fruit_cnd, 
                       veg_cnd, 
                       fish_cnd, 
                       diary_cnd,
                       veg_oil_cnd,
                       whole_grain_cnd,
                       refine_grain_cnd,
                       pmeat_cnd, 
                       npmeat_cnd, 
                       sugar_cnd)
diet_cnd <- rowSums(diet_mat) >= 5
# Can satisfy even when missing = F
diet_cnd[rowSums(diet_mat,na.rm = T) >= 5] <- T
# Never satisfy even when missing = T
diet_cnd[rowSums(is.na(diet_mat)) + rowSums(diet_mat,na.rm = T) < 5] <- F
diet_out <- ifelse(diet_cnd, "Yes", "No") %>%
  factor(., levels = c("No", "Yes"))

################## alcohol ##################
# women ??1 drink/day; men ??2 drinks/day
# "p1558_i0", 	 ## 05.01 Alcohol intake frequency. 
# "p4407_i0", 	 ## 05.02 Average monthly red wine intake
# "p4418_i0", 	 ## 05.03 Average monthly champagne plus white wine intake
# "p4429_i0", 	 ## 05.04 Average monthly beer plus cider intake
# "p4440_i0", 	 ## 05.05 Average monthly spirits intake
# "p4451_i0", 	 ## 05.06 Average monthly fortified wine intake
# "p4462_i0", 	 ## 05.07 Average monthly intake of other alcoholic drinks
# "p1568_i0", 	 ## 05.08 Average weekly red wine intake
# "p1578_i0", 	 ## 05.09 Average weekly champagne plus white wine intake
# "p1588_i0", 	 ## 05.10 Average weekly beer plus cider intake
# "p1598_i0", 	 ## 05.11 Average weekly spirits intake
# "p1608_i0", 	 ## 05.12 Average weekly fortified wine intake
# "p5364_i0" 	   ## 05.13 Average weekly intake of other alcoholic drinks
# never alcohol use
noalcohol_cnd <- ifelse(trait_exposure_dat[["p1558_i0"]] == "Never", T,
                        ifelse(trait_exposure_dat[["p1558_i0"]] %in% c("", "Do not know", "Prefer not to answer"), NA, F))
#
alcohol_intake_month_mat <- lapply(c("p4407_i0", "p4418_i0", "p4429_i0", "p4440_i0", "p4451_i0", "p4462_i0"), function(x){
  as.numeric(trait_exposure_dat[[x]])
}) %>% Reduce("cbind", .)
alcohol_intake_week_mat <- lapply(c("p1568_i0", "p1578_i0", "p1588_i0", "p1598_i0", "p1608_i0", "p5364_i0"), function(x){
  as.numeric(trait_exposure_dat[[x]])
}) %>% Reduce("cbind", .)
alcohol_intake_day_mean_month <- rowSums(alcohol_intake_month_mat, na.rm = T) / 30
alcohol_intake_day_mean_month[rowMeans(is.na(alcohol_intake_month_mat)) == 1] <- NA
alcohol_intake_day_mean_week <- rowSums(alcohol_intake_week_mat, na.rm = T) / 7
alcohol_intake_day_mean_week[rowMeans(is.na(alcohol_intake_week_mat)) == 1] <- NA
alcohol_intake_day <- data.frame(alcohol_intake_day_mean_month,
                                 alcohol_intake_day_mean_week) %>%
  rowMeans(., na.rm = T)
alcohol_intake_day[is.na(alcohol_intake_day_mean_month) & is.na(alcohol_intake_day_mean_week)] <- NA

alcohol_intake_day_weight <- ifelse(trait_exposure_dat[["sex"]] == "female",
                                    alcohol_intake_day * 2, alcohol_intake_day)
noalcohol_cnd[which(alcohol_intake_day_weight <= 2)] <- T
noalcohol_out <- ifelse(diet_cnd, "Yes", "No") %>%
  factor(., levels = c("No", "Yes"))

################## life factors mat format ##################
lf_df <- data.frame(Nosmoke = nosmoke_cnd, 
                    Activity = act_cnd, 
                    Diet = diet_cnd, 
                    Noalcohol = noalcohol_cnd)
lifescore <- rowSums(lf_df) >= 3
lifescore[rowSums(lf_df, na.rm = T) >= 3] <- T
lifescore[rowSums(is.na(lf_df)) + rowSums(lf_df,na.rm = T) < 3] <- F
lifescore_out <- ifelse(diet_cnd, "High", "Low") %>%
  factor(., levels = c("Low", "High"))


#######################################
env_code <- c("p24003",            # 1.1 Nitrogen dioxide air pollution; 2010
              "p24004",            # 1.2 Nitrogen oxides air pollution; 2010
              "p24005",            # 1.3 Particulate matter air pollution (pm10); 2010
              "p24006",            # 1.4 Particulate matter air pollution (pm2.5); 2010
              "p24007",            # 1.5 Particulate matter air pollution (pm2.5) absorbance; 2010
              "p24008",            # 1.6 Particulate matter air pollution 2.5-10um; 2010
              "p24009",            # 1.7 Traffic intensity on the nearest road
              "p24010",            # 1.8 Inverse distance to the nearest road
              "p24011",            # 1.9 Traffic intensity on the nearest major road
              "p24012",            # 1.10 Inverse distance to the nearest major road
              "p24013",            # 1.11 Total traffic load on major roads
              "p24014",            # 1.12 Close to major road
              "p24015",            # 1.13 Sum of road length of major roads within 100m
              "p24016",            # 1.14 Nitrogen dioxide air pollution; 2005
              "p24017",            # 1.15 Nitrogen dioxide air pollution; 2006
              "p24018",            # 1.16 Nitrogen dioxide air pollution; 2007
              "p24019",            # 1.17 Particulate matter air pollution (pm10); 2007
              "p24020",            # 1.18 Average daytime sound level of noise pollution
              "p24021",            # 1.19 Average evening sound level of noise pollution
              "p24022",            # 1.20 Average night-time sound level of noise pollution
              "p24023",            # 1.21 Average 16-hour sound level of noise pollution
              "p24024",            # 1.22 Average 24-hour sound level of noise pollution
              "p24500_i0",         # 2.1 Greenspace percentage, buffer 1000m 
              "p24501_i0",         # 2.2 Domestic garden percentage, buffer 1000m
              "p24502_i0",         # 2.3 Water percentage, buffer 1000m
              "p24503_i0",         # 2.4 Greenspace percentage, buffer 300m 
              "p24504_i0",         # 2.5 Domestic garden percentage, buffer 300m
              "p24505_i0",         # 2.6 Water percentage, buffer 300m
              "p24506_i0",         # 2.7 Natural environment percentage, buffer 1000m
              "p24507_i0",         # 2.8 Natural environment percentage, buffer 300m
              "p24508_i0"          # 2.9 Distance (Euclidean) to coast
)
trait_env_dat <- fread2(TRAIT, select = c("eid", env_code))
trait_env_dat <- trait_env_dat[match(all_id, trait_env_dat$eid), ]
trait_env_dat$eid <- NULL
colnames(trait_env_dat) <- c(paste0("PHE0100", c(1: 9)), paste0("PHE010", c(10: 22)), 
                             paste0("PHE0200", c(1: 9)))

################### life data frame ########################
exposure_dat <- cbind.data.frame(trait_env_dat, 
                                 PHE03001 = incm_cnd, 
                                 PHE03002 = edu_cnd, 
                                 PHE03003 = emp_cnd,
                                 PHE03004 = SES_out, 
                                 PHE03005 = trait_exposure_dat[["p22189"]],
                                 PHE04001 = crt_smk_cnd, 
                                 PHE04002 = pst_smk_cnd,
                                 PHE04003 = npst_smk_cnd, 
                                 PHE04004 = nosmoke_cnd, 
                                 PHE04005 = act_cnd, 
                                 PHE04006 = fruit_cnd, 
                                 PHE04007 = veg_cnd, 
                                 PHE04008 = fish_cnd, 
                                 PHE04009 = diary_cnd,
                                 PHE04010 = veg_oil_cnd,
                                 PHE04011 = whole_grain_cnd,
                                 PHE04012 = refine_grain_cnd,
                                 PHE04013 = pmeat_cnd, 
                                 PHE04014 = npmeat_cnd, 
                                 PHE04015 =sugar_cnd, 
                                 PHE04016 = diet_cnd, 
                                 PHE04017 = noalcohol_cnd, 
                                 PHE04018 = lifescore_out)
saveRDS(exposure_dat, file = paste0(PROJ_PATH, "/02_data/out_pheno/exposure_factors_all.rds"))
saveRDS(exposure_dat[all_id%in%male_id, ], 
        file = paste0(PROJ_PATH, "/02_data/out_pheno/exposure_factors_male.rds"))
saveRDS(exposure_dat[all_id%in%female_id, ], 
        file = paste0(PROJ_PATH, "/02_data/out_pheno/exposure_factors_female.rds"))