#------------------------------------------------------------------------------#
#                                                                              #
# Nonparametric analysis of illness onset data for 2019-nCoV in Wuhan, China   #
# from Liu et al. (2020) (n=6 transmission pairs)                              #
# Oswaldo Gressani. All rights reserved.                                       #
#                                                                              #
#------------------------------------------------------------------------------# 

rm(list = ls())

# Source: Li, Q., Guan, X., et al. (2020). Early transmission dynamics in Wuhan, 
# China, of novel coronavirus–infected pneumonia. 
# New England journal of medicine, 382(13), 1199-1207.

infector_dates <- data.frame(Case_id_infector = c("Case 1.1","Case 1.1", "Case 2.1", 
                                            "Case 3.1", "Case 4.1", "Case 5.1"),
                             Calendar_date = as.Date(c("2019/12/20", "2019/12/20",
                                                       "2019/12/27", "2019/12/12",
                                                       "2019/12/21", "2020/01/04")))

infectee_dates <- data.frame(Case_id_infector = c("Case 1.2","Case 1.3", "Case 2.2", 
                                                  "Case 3.3", "Case 4.3", "Case 5.2"),
                             Calendar_date = as.Date(c("2019/12/25", "2019/12/29",
                                                       "2020/01/03", "2019/12/19",
                                                       "2019/12/24", "2020/01/11")))

# siobs <- infectee_dates$Calendar_date-infector_dates$Calendar_date

# Making data doubly interval-censored
infectee_left  <- as.numeric(infectee_dates$Calendar_date)
infectee_right <- infectee_left + 1
infector_left  <- as.numeric(infector_dates$Calendar_date)
infector_right <- infector_left + 1
sl <- infectee_left - infector_right
sr <- infectee_right - infector_left

SIdatnCoV2019 <- data.frame(Pair_index = seq_len(6), sl = sl, sr = sr)

# Export data to local file in txt format
write.table(SIdatnCoV2019, file = "Cov19-Li-2020.txt")