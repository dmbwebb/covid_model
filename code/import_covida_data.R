
# Import COVIDA data ------------------------------------------------------


#Clean the workspace
# rm(list=ls())
cat("\014")
local({r <- getOption("repos"); r["CRAN"] <- "http://cran.r-project.org"; options(repos=r)}) #set repo


#Load Packages
pkg<-list("dplyr","haven","lubridate",'openxlsx')
lapply(pkg, require, character.only=T)
rm(pkg)







# load Data ------------------------------------------------------------------
dta<-read_dta("data/Datos_Salesforce_treated_mar03_slim.dta")

# dta %>% print_names
# 
# dta <- dta %>% 
#   select(
#     fechatomamuestra,
#     fecharecepciónmuestralab,
#     fechahoradeaperturaalone1,
#     contact_COVID,
#     contact,
#     symptom,
#     stratum,
#     positive,
#     ocup_cat,
#     test_month,
#     weight_ocup,
#     númerodepersonasconquienesvive
#   ) %>% 
#   write_dta("data/Datos_Salesforce_treated_mar03_slim.dta")

#According to my encoding fecharecepci? always comes as  fecharecepci?nmuestralab. Change in CPu if not your case
dta_covida <- dta %>%
  mutate(test_day=as.Date(fechatomamuestra, "%d/%m/%Y"),
         #lab_reception=as.Date(fecharecepci?nmuestralab, "%d/%m/%Y"),#, origin="1960-01-01"), #for Camilo PC
         #lab_reception=as.Date(fecharecepciónmuestralab, "%d/%m/%Y"),#, origin="1960-01-01"), #for Camilo PC
         lab_reception=as.Date(fecharecepciónmuestralab,"%d/%m/%Y"), #, origin="1960-01-01"), #for Camilo PC
         fechahoradeaperturaalone1=as.Date(fechahoradeaperturaalone1, "%d/%m/%Y"),#, origin="1960-01-01"),
         test_day= ifelse(is.na(test_day), lab_reception,test_day), #fix test_day with fecharecepciónmuestralab
         test_day= ifelse(is.na(test_day), fechahoradeaperturaalone1,test_day), #fix test_day with fecharecepciónmuestralab
         test_day= as_date(test_day),
         exclude=ifelse( symptom==1 | contact_COVID==1 | contact==1,1,0),
         mes=month(test_day),
         year=year(test_day),
         semana=week(test_day),
         stratum = case_when(stratum %in% c(1, 2) ~ 1,
                             stratum %in% c(3, 4) ~ stratum - 1,
                             stratum %in% c(5, 6) ~ 4),
         #date_m=dmy(paste("01",mes,year,sep="-")),
         date_m=floor_date(test_day, "month"),
         semana=ifelse(semana==53,52,semana),
         date_week=paste(year,semana,"1",sep="-"),
         date_week=as.Date(date_week, "%Y-%U-%u")
  ) %>%
  filter_track(date_m>as.Date("2020-04-01")) %>%  #drops 3 obs in may, we start in june
  filter_track(!is.na(positive)) %>%  #Missing test result (18453 obs, we end up with)57165 obs
  filter(exclude==0) %>%  #Exclude those with symtoms and contacts (40292 obs)
  filter_track(!(ocup_cat=="militares y fuerza publica" &  test_day==as.Date("2020-07-02"))) %>% 
  filter()

# We need this for lal the figures and tables....
# dta_covida<- dta_covida %>% 
#   mutate(ocup_cat=ifelse(grepl("aseo",ocupacion_desagregada),"personal limpieza",ocup_cat),
#          poblacion_desagregada=ifelse(grepl("aseo",ocupacion_desagregada),108800,poblacion_desagregada),
#          ocup_cat=ifelse(grepl("gruesa",ocupacionasis),	"obreros de construccion",ocup_cat	),
#          poblacion_desagregada=ifelse(grepl("gruesa",ocupacion_desagregada),306000,poblacion_desagregada),
#          ocupacion_desagregada=ifelse(grepl("gruesa",ocupacionasis),"obreros de construccion",ocupacion_desagregada	),
#          ocup_cat=ifelse(ocupacionasis=="constructores de casas",	"obreros de construccion",ocup_cat	),
#          ocupacion_desagregada=ifelse(ocupacionasis=="constructores de casas",	"obreros de construccion",ocupacion_desagregada	)
#   )


# write_dta(dta_covida, "Data/Datos_Salesforce_treated_feb19_clean.dta")







# Fig 2b ------------------------------------------------------------------





#Load Packages
pkg<-list("dplyr","ggplot2","stringr","haven",'tidyr', "lubridate","ggsci")
lapply(pkg, require, character.only=T)
rm(pkg)





# Parameters --------------------------------------------------------------
set.seed(101010)
name<-"analytic"
#days_oct<-as.numeric(dmy("30-11-2020")-dmy("01-06-2020"))
days_oct<-30*5
#days_fin<-as.numeric(dmy("03-03-2021")-dmy("01-06-2020"))
days_fin<-30*9

# covida ------------------------------------------------------------------


# dta_covida<-read_dta("Data/Datos_Salesforce_treated_feb19_clean.dta") 


# poblacion<-read_dta("Data/pob_strat.dta")
# poblacion <- poblacion %>% 
#   rename(poblacion_agregada=pob_stratum)


# REPLACE THIS WITH STRATUM POPS (from import_data.R)
sds_jobs <- read_dta(str_glue("data/casos_SDS_poblaciones_20Feb2021_slim.dta"))
# sds_jobs %>% print_names
# Clean up
sds_jobs_clean <- sds_jobs %>% 
  # sample_n(100) %>% 
  rename(
    any_of(c(
      "case_id" = "caso",
      "date_symptoms" = "fechainiciosintomas",
      "date_consultation" = "fecha_consulta",
      "date_results" = "fechadiagnostresultlaboratorio",
      "stratum" = "estratosocioeconomico",
      "stratum_2" = "estrato",
      "stratum_pop" = "poblacion_estrato",
      "recovered" = "recupsaluddatabog",
      "date_death" = "fecharecuperadosaluddata"
    ))
  ) %>% 
  mutate(stratum = coalesce(as.numeric(stratum), stratum_2)) %>% 
  select(-stratum_2) %>% 
  count_prop(stratum) %>% 
  mutate_track(across(stratum, ~ if_else(as.numeric(.x) %in% 1:6, as.numeric(.x), NA_real_))) %>% 
  mutate(across(starts_with("date"), dmy)) %>% 
  count_prop(stratum, stratum_pop) %>% 
  count_prop(year(date_symptoms), year(date_results), year(date_consultation))


stratum_pops <- sds_jobs_clean %>% 
  # count_prop(stratum, stratum_pop)
  mutate(stratum = as.integer(stratum)) %>% 
  mutate(stratum = case_when(stratum %in% c(1, 2) ~ 1L,
                             stratum %in% c(3, 4) ~ stratum - 1L,
                             stratum %in% c(5, 6) ~ 4L)) %>% 
  group_by(stratum, stratum_pop) %>% 
  summarise() %>% 
  filter(!is.na(stratum), !is.na(stratum_pop)) %>% 
  summarise(stratum_pop = sum(stratum_pop))  



poblacion <- stratum_pops %>% 
  rename(poblacion_agregada = stratum_pop)


#db<-dta_covida



# Calculate rates by strata --------------------------------------------
#June November
rates_oct <-broom::tidy(lm(positive~as.factor(stratum)-1,dta_covida %>%   filter(mes>4 & mes<12) ,weights = weight_ocup), conf.int = TRUE)


rates_oct<- rates_oct %>% 
  mutate(term=str_remove_all(term,"as.factor\\(stratum\\)"),
         term=as.numeric(term)) %>% 
  rename(stratum=term,
         rate_pos=estimate,
         q025=conf.low,
         q975=conf.high)  %>% 
  select(stratum,rate_pos,q025,q975) %>%
  # left_join(stratum_pops) %>% 
  # rename(poblacion_agregada = stratum_pop) %>%
  left_join(.,poblacion) %>%
  mutate(tot_day_cases_covida=((rate_pos*poblacion_agregada)/17)*days_oct,
         q025_tot_day_cases_covida=((q025*poblacion_agregada)/17)*days_oct,
         q975_tot_day_cases_covida=((q975*poblacion_agregada)/17)*days_oct) %>%
  mutate(acumm_covid_covida=tot_day_cases_covida/poblacion_agregada,
         q025_acumm_covid_covida=q025_tot_day_cases_covida/poblacion_agregada,
         q975_acumm_covid_covida=q975_tot_day_cases_covida/poblacion_agregada,
         grp=1) %>%
  na.omit()





#June March
rates_jan <-broom::tidy(lm(positive~as.factor(stratum)-1,dta_covida ,weights = weight_ocup), conf.int = TRUE)

rates_jan<- rates_jan %>% 
  mutate(term=str_remove_all(term,"as.factor\\(stratum\\)"),
         term=as.numeric(term)) %>% 
  rename(stratum=term,
         rate_pos=estimate,
         q025=conf.low,
         q975=conf.high)  %>% 
  select(stratum,rate_pos,q025,q975) %>%
  left_join(.,poblacion)%>%
  mutate(tot_day_cases_covida=((rate_pos*poblacion_agregada)/17)*days_fin,
         q025_tot_day_cases_covida=((q025*poblacion_agregada)/17)*days_fin,
         q975_tot_day_cases_covida=((q975*poblacion_agregada)/17)*days_fin) %>%
  mutate(acumm_covid_covida=tot_day_cases_covida/poblacion_agregada,
         q025_acumm_covid_covida=q025_tot_day_cases_covida/poblacion_agregada,
         q975_acumm_covid_covida=q975_tot_day_cases_covida/poblacion_agregada,
         grp=2) %>%
  na.omit()


covida_rates <-bind_rows(rates_oct,rates_jan)





# Positivity rates by hh --------------------------------------------------

# PLOT 
dta_covida %>% select(matches("hogar")) %>% print_names


covida_hhsize <- dta %>% 
  mutate(
    stratum = case_when(stratum %in% c(1, 2) ~ 1,
                        stratum %in% c(3, 4) ~ stratum - 1,
                        stratum %in% c(5, 6) ~ 4)
  ) %>% 
  count_prop(númerodepersonasconquienesvive) %>% 
  mutate(hhsize = as.integer(str_replace_all(númerodepersonasconquienesvive, "Más de 10", "10"))) %>% 
  count_prop(hhsize)


broom::tidy(lm(positive ~ factor(stratum) + factor(stratum)*hhsize + factor(test_month) - 1, data = covida_hhsize))

covida_hhsize %>% 
  filter(!is.na(stratum)) %>% 
  count_prop(stratum) %>% 
  mutate(i_group = recode_i_group(stratum)) %>% 
  mutate(stratum_merge = case_when(stratum %in% 1:2 ~ 1, 
                                   stratum %in% 3:4 ~ 2)) %>% 
  ggplot(aes(x = hhsize, y = positive, colour = i_group, fill = i_group)) + 
  geom_smooth(method = "lm", se = FALSE, size = 2) + 
  geom_ribbon(method = "lm", linetype = "dotted", stat = "smooth", alpha = .05) + 
  theme_classic() + 
  scale_colour_viridis(option = "viridis",
                       begin = 0.3, 
                       discrete = TRUE) + 
  scale_fill_viridis(option = "viridis",
                     begin = 0.3, 
                     discrete = TRUE) +
  scale_y_continuous(labels = scales::percent) + 
  scale_x_continuous(breaks = 1:10) + 
  coord_cartesian(ylim = c(0, 0.13), xlim = c(1, 10)) + 
  labs(x = "Household size", y = "Proportion positive (%)", colour = "SES Group", fill = "SES Group")

ggsave("figures/covida_hhsize.pdf", device = cairo_pdf, width = 7, height = 5, scale = 0.7)


# covida_hhsize %>% 
#   filter(!is.na(stratum)) %>% 
#   ggplot(aes(x = hhsize, fill = factor(stratum))) + geom_histogram() + 
#   facet_wrap(~factor(stratum), scales = "free_y")
# 
# broom::tidy(lm(positive~as.factor(stratum)-1,
#                dta_covida,
#                weights = weight_ocup), conf.int = TRUE)
# 






rm(dta, dta_covida, rates_oct, rates_jan)




# covida_rates<- covida_rates %>% mutate(grp=factor(grp, levels=c(1,2),  labels=c("November 30th","March 3rd"),ordered = TRUE),
#                    stratum=factor(stratum,levels=c(1,2,3,4), labels=c("1&2","3","4","5&6"),ordered = TRUE))
# 
# 
# covida_rates<- covida_rates %>% mutate(acumm_covid_covida=acumm_covid_covida*100,
#                    q025=q025_acumm_covid_covida*100,
#                    q975=q975_acumm_covid_covida*100) %>%
#   mutate(q975=ifelse(q975>100,100,q975),
#          q025=ifelse(q025<0,0,q025),
#   )
# 
# 
# 
# ggplot(data=covida_rates, aes(x=stratum, y=acumm_covid_covida, group=grp, col=grp))+
#   geom_point(size=1, position=position_dodge(width = .2))+
#   geom_errorbar(aes(ymin=q025, ymax=q975), width=.1, position=position_dodge(width = .2)) +
#   xlab("Socioeconomic Strata") +
#   theme_bw() +
#   scale_y_continuous("Accumulated SARS-COV-2 Cases \n as Percentage of Population",breaks =seq(0,100,20),limits=c(0,102)) +
#   theme(legend.title= element_blank() ,
#         legend.position="bottom",
#         legend.text=element_text(size=14),
#         axis.title = element_text(size=14),
#         panel.grid.major.x = element_blank(),
#         legend.background = element_rect(fill='transparent'),
#         axis.text.x =element_text( angle=0,hjust=0.5,size=14),
#         axis.text.y =element_text( size=12),
#         rect = element_rect(colour = "transparent", fill = "white")
#   ) + scale_color_manual(values=c("#3B4992B2","#EE0000B2"))




# ggsave(paste0("views/Fig2_b_",name,".pdf"),height=5,width=7)






# Probability of being group X cond on infected ---------------------------


# crossing(
#   date = c("december", "march"),
#   stratum = 1:4
# ) %>% 
#   mutate(
#     model = 
#   )
# 
# 
# data_strata <- dta_covida %>% mutate(
#   stratum1 = stratum == 1,
#   stratum2 = stratum == 2,
#   stratum3 = stratum == 3,
#   stratum4 = stratum == 4
# ) %>% 
#   count_prop(positive) %>% 
#   filter_track(positive == 1)
# 
# broom::tidy(lm(stratum1 ~ 1, data_strata %>%   filter(mes>4 & mes<12), weights = weight_ocup), conf.int = TRUE)
# broom::tidy(lm(stratum4 ~ 1, data_strata, weights = weight_ocup), conf.int = TRUE)
# 
# broom::tidy(lm(positive~as.factor(stratum)-1,dta_covida %>%   filter(mes>4 & mes<12) ,weights = weight_ocup), conf.int = TRUE)












