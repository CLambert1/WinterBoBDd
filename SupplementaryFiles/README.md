Supplementary Files 
==========
for the paper **"Delayed response to environmental conditions and infra-seasonal dynamics of the short-beaked common dolphin distribution"**, by C. Lambert, M. Authier, S. Laran and J. Spitz. 


In this repository, you will find supplementary figures displaying the daily variations in environmental conditions over the course of the season, from December 1st to April 13th in 2019-2020 and 2020-2021. These figures are presented for all dynamic variables: Temp, Salinity, CurrentSpeed, Chl, NPP, Phyto, EKE, DissIC, SPCO2, SSH, MLD, ZEu, gradSSH, gradTemp, gradMLD, gradSal, gradChl, gradNPP, Dist MLD front, Dist Chl front, Dist NPP front, Dist Temp front. We also present the spatial distribution of static variables, namely the persistence of MLD, Chl, NPP and Temp fronts.

Day-to-day predictions of common dolphins are also shown, alongside the model-derived uncertainty associated to predictions. Day-to-day maps of extrapolations are also presented. 

The uncertainty and extrapolations maps derived from the seasonal model are also displayed. 

Environmental conditions - maps
---------

The time series of environmental conditions summarized over the study area for each day from December 1st to April 13th show the strong variations occurring over the season in the oceanography of the Bay of Biscay. Both the mean and the standard deviation of the conditions are displayed, in black for the winter 2019-2020, in gold for the winter 2020-2021:
![Variable Time Series](https://github.com/CLambert1/WinterBoBDd/blob/main/SupplementaryFiles/Variable_timeSeries.png)

Animated maps displaying the full range of maps over the period are visible in the folder, under the name "Variable_XXX". These animations show daily maps for the two winters with the time series running alongside. 


Predictions of common dolphin distribution - Daily model
---------
For the purpose of the study, we built a daily model of dolphin distribution, providing daily predictions. This model was our reference, used for the three main objectives of the study. Here, we present the day-to-day predicted maps (ind./m²) alongside the time series of averaged estimated densities (in black solid line), the associated averaged uncertainty (in black dotted line) and the estimated extrapolation level of the prediction map (in red solid line). The density and uncertainty time series are the averages of the daily maps.
![GIF daily predictions](https://github.com/CLambert1/WinterBoBDd/blob/main/SupplementaryFiles/DailyModel_Prediction_Uncertainty_Winter.gif)

We also computed the degree of extrapolation for each day (see main text), as a mean to quantify the distance of the prediction from the data used to fit the model (with regards to the sampled conditions). The dark grey areas are interpolations areas (that is, the occurring conditions have been sampled in the model), the lighter grey areas are extrapolation areas (that is, the occurring conditions have *not* been sampled):  
![GIf daily extrapolation](https://github.com/CLambert1/WinterBoBDd/blob/main/SupplementaryFiles/DailyModel_Extrapolation_Winter.gif)

From these daily extrapolation maps, we can see predictions within some areas are systematically interpolations (in blue), while some are always extrapolated values: 
![Mean map extrapolation daily model](https://github.com/CLambert1/WinterBoBDd/blob/main/SupplementaryFiles/DailyModel_Extrapolation_mean_maps.png)


Predictions of common dolphin distribution - Seasonal model
---------
For the purpose of the study, we also fitted a challenger model assuming no spatio-temporal variations in dolphin distribution throughout the winter (see main text for details). This model was fitted using winter-scale means and standard deviations of the same variables as depicted above. We can see the dolphin distribution (ind./km²) is predicted to be different during the two winters. The figure also display the model-associated uncertainty (middle panels) and the extrapolation level (right panels; dark grey = interpolation; light grey = extrapolation): 
![Extrapolation from the seasonal model](https://github.com/CLambert1/WinterBoBDd/blob/main/SupplementaryFiles/Seasonal_model_maps.png)



