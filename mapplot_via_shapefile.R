library(rgdal)
library(ggplot2)
library(spdep)
library(nimble)
library(ngspatial)
library(sf)
source("Functions.R")

#================
shp_states <- rgdal::readOGR(dsn="shp5m",layer = "cb_2018_us_state_5m")
summary(shp_states@data)
shp_48statesAndDC <- shp_states[-c(40,48,49,50,52,56,28),]
summary(shp_48statesAndDC@data)

shp2 <- read_sf("cb_2018_us_state_5m.shp")
#keep 48 states and DC only
shp2_48andDC <- shp2[-c(40,48,49,50,52,56,28),]
#================Draw the map without attributes
map_states <- ggplot()+geom_polygon(data = shp_48statesAndDC,
                                    aes(x=long,
                                        y=lat,
                                        group=group),
                                    colour="black") + theme_void()
#================read covariates
df2 = read.csv("by48state0527.csv")
#================Draw the map with discretized pop density
shp2_48andDC$Pop = df2$Pop
shp2_48andDC$PopDens = df2$Popdensity %>% cut_number(n=7)
#p1 <- ggplot(shp2_48andDC)+geom_sf(aes(fill=Pop))+ theme_void() + theme(legend.position = "none")
p_PopDens <- ggplot(shp2_48andDC)+geom_sf(color="gray90",size=0.1,
                                   aes(fill=PopDens)) + theme_void() + theme(legend.position = "bottom")
p_PopDens <- p_PopDens + scale_fill_brewer(palette="Blues",
                             labels=c("6-24","25-55","57-76","88-112","148-206","212-420",
                                      "485-11011"))
p_PopDens=p_PopDens + labs(fill="Population per\nsquare mile") + guides(fill=guide_legend(nrow=2))

#================Draw the map with air pollution
shp2_48andDC$AirPol = df2$AirPol %>% cut_number(n=7)
df2$AirPol[df2$AirPol %>% order()]
p_AirPol <- ggplot(shp2_48andDC)+geom_sf(color="gray90",size=0.1,
                                   aes(fill=AirPol)) + theme_void() + theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Reds",
                    labels = levels(shp2_48andDC$AirPol)) +
  labs(fill="PM2.5 measured in \nmg per cubic meter") + guides(fill=guide_legend(nrow=2))
#================with Obesity percentage
shp2_48andDC$Obesity = df2$Obesity %>% cut_number(n=7)
df2$Obesity[df2$Obesity %>% order()]
p_Obesity <- ggplot(shp2_48andDC)+geom_sf(color="gray90",size=0.1,
                                   aes(fill=Obesity)) + theme_void() + theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Greens",
                    labels = levels(shp2_48andDC$Obesity)) +
  labs(fill="Percentage of adults\nwith BMI >= 30") + guides(fill=guide_legend(nrow=2))

#================Drug_death
shp2_48andDC$Drug_death = df2$Drug_death %>% cut_number(n=7,boundary=0.5)
df2$Drug_death[df2$Drug_death %>% order()]
levels(shp2_48andDC$Drug_death)
p_Drug_death <- ggplot(shp2_48andDC)+geom_sf(color="gray90",size=0.1,
                                   aes(fill=Drug_death)) + theme_void() + theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Purples",
                    labels = levels(shp2_48andDC$Drug_death)) +
  labs(fill="Number of deaths due to \ndrug injury per 100,000") + guides(fill=guide_legend(nrow=2))

gridExtra::grid.arrange(p_PopDens,p_Physician,p_Uninsured,p_beds,ncol=2)
gridExtra::grid.arrange(p_Smoking,p_Excessive_Drinking,p_Drug_death,p_Obesity,
                        p_Inactivity,p_MDI,ncol=2)

#================
shp2_48andDC$Smoking = df2$Smoking %>% cut_number(n=7,boundary=0.5)
#df2$Smoking[df2$Smoking %>% order()]
#levels(shp2_48andDC$Smoking)
p_Smoking <- ggplot(shp2_48andDC)+geom_sf(color="gray90",size=0.1,
                                   aes(fill=Smoking)) + theme_void() + theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "PuRd",
                    labels = levels(shp2_48andDC$Smoking)) +
  labs(fill="Percentage of adults \nwho are smokers") + guides(fill=guide_legend(nrow=2))

shp2_48andDC$Excessive_Drinking = df2$Excessive_Drinking %>% cut_number(n=7,boundary=0.5)
#================
#df2$Excessive_Drinking[df2$Excessive_Drinking %>% order()]
#levels(shp2_48andDC$Excessive_Drinking)
p_Excessive_Drinking <- ggplot(shp2_48andDC)+geom_sf(color="gray90",size=0.1,
                                   aes(fill=Excessive_Drinking)) + theme_void() + theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "OrRd",
                    labels = levels(shp2_48andDC$Excessive_Drinking)) +
  labs(fill="Percentage of adults \nwho are excessive drinking") + guides(fill=guide_legend(nrow=2))
#================
shp2_48andDC$Physician = df2$Physician %>% cut_number(n=7,boundary=0.5)
#df2$Physician[df2$Physician %>% order()]
#levels(shp2_48andDC$Physician)
p_Physician <- ggplot(shp2_48andDC)+geom_sf(color="gray90",size=0.1,
                                   aes(fill=Physician)) + theme_void() + theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "BuGn",
                    labels = levels(shp2_48andDC$Physician)) +
  labs(fill="Number of active primary care \nphysicians per 100,000") + guides(fill=guide_legend(nrow=2))
p_Physician
colnames(df2)
#[1] "State"              "AirPol"             "Obesity"            "Drug_death"        
#[5] "Excessive_Drinking" "Smoking"            "Physician"          "Uninsured"         
#[9] "Inactivity"         "beds"               "OPM"                "MDI"               
#[13] "Popdensity"         "Pop"                "Abbre"              "Confirmed"         
#[17] "Neg"                "pending"         
#================
shp2_48andDC$Uninsured = df2$Uninsured %>% cut_number(n=7,boundary=0.5)
#df2$Uninsured[df2$Uninsured %>% order()]
#levels(shp2_48andDC$Uninsured)
p_Uninsured <- ggplot(shp2_48andDC)+geom_sf(color="gray90",size=0.1,
                                            aes(fill=Uninsured)) + theme_void() + theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "BuPu",
                    labels = levels(shp2_48andDC$Uninsured)) +
  labs(fill="Percentage of the population \nwithout health insurance") + guides(fill=guide_legend(nrow=2))
p_Uninsured

#================
shp2_48andDC$Inactivity = df2$Inactivity %>% cut_number(n=7,boundary=0.5)
#df2$Inactivity[df2$Inactivity %>% order()]
#levels(shp2_48andDC$Inactivity)
p_Inactivity <- ggplot(shp2_48andDC)+geom_sf(color="gray90",size=0.1,
                                            aes(fill=Inactivity)) + theme_void() + theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Reds",
                    labels = levels(shp2_48andDC$Inactivity)) +
  labs(fill="Percentage of adults \nreporting physical inactivity") + guides(fill=guide_legend(nrow=2))
p_Inactivity

#================
shp2_48andDC$beds = df2$beds %>% cut_number(n=7,boundary=0.5)
#df2$beds[df2$beds %>% order()]
#levels(shp2_48andDC$beds)
p_beds <- ggplot(shp2_48andDC)+geom_sf(color="gray90",size=0.1,
                                             aes(fill=beds)) + theme_void() + theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Purples",
                    labels = levels(shp2_48andDC$beds)) +
  labs(fill="Hospital beds per 1,000") + guides(fill=guide_legend(nrow=2))
p_beds
#=================
shp2_48andDC$MDI = df2$MDI %>% cut_number(n=7,boundary=0.5)
#df2$MDI[df2$MDI %>% order()]
#levels(shp2_48andDC$MDI)
p_MDI <- ggplot(shp2_48andDC)+geom_sf(color="gray90",size=0.1,
                                       aes(fill=MDI)) + theme_void() + theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Blues",
                    labels = levels(shp2_48andDC$MDI)) +
  labs(fill="Multidimensional Deprivation Index") + guides(fill=guide_legend(nrow=2))
p_MDI

rook.nb = poly2nb(shp_48statesAndDC,queen = F)
A = nb2mat(rook.nb,style = "B")
n_adj=rowSums(A)
D=diag(n_adj)
adj=as.carAdjacency(A)$adj
#State = tibble(State=as.character(shp_48statesAndDC@data$NAME))
State = tibble(State = shp_48statesAndDC@data$NAME)
df = read.csv("bystate0527.csv")
df2 = left_join(State,df,by="State")
write_csv(df2,"by48state0527.csv")

#===========================
shp_48statesAndDC@data$Pop = df2$Pop


