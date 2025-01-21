##################################
## Species distribution modelling##
##################################

## 3rd Year computing practical
## tim.barraclough@biology.ox.ac.uk

## Code based on - https://rspatial.org/raster/sdm/

## Summary: This annotated R script shows you how to download a) map of the world, b) locality data for
## a chosen species from the Global Biodiversity Information Facility,GBIF https://www.gbif.org,
## and c) climate data from WorldClim https://www.worldclim.org/data/worldclim21.html
## and then to run analyses modelling the distribution of your chosen species in relation to climate. The models estimate
## which climatic variables predict current distribution of your species based on current climate, and  then to predict
## future distribution based on predict changes in climate variables into the future. This is a widely used approach for
## understanding and predicting species distributions, which can be applied at a range of scales from local to global.
## NB - there are lots of methods for species distribution modelling, we just use a simple linear model approach.


# Create project structure if it doesn't exist
project_dirs <- c("data/raw", "data/processed", "output")
sapply(project_dirs, function(dir) {
    dir_path <- here(dir)
    if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
})

## Every time you run it, need to load the packages
library(here)
library(dismo)
library(rworldmap)
library(sf)
library(geodata)

#################################
## 1) Load world map and plot it##
#################################

wrld_simpl <- getMap(resolution = "coarse")
plot(wrld_simpl)

## yay, a map!


##############################################################
## 2) Download locality data for a species from GBIF and trim##
##############################################################

## access records from GBIF - includes lots of columns, most are not useful
## I am looking at Coffea arabica, the coffee species used for arabica coffee
species.gbif_file <- here("data", "raw", "species.gbif.rds")
if (!file.exists(species.gbif_file)) {
    species.gbif <- gbif("coffea", "arabica", geo = TRUE)
    saveRDS(species.gbif, species.gbif_file)
} else {
    species.gbif <- readRDS(species.gbif_file)
}

## pulls out the lat and lon columns => there is more info if you want it, but we don't need it here
# and delete rows with missing data = NA
species.coords <- cbind(species.gbif$lon, species.gbif$lat)
species.coords <- na.omit(species.coords)
species.coords <- data.frame(species.coords)
colnames(species.coords) <- c("lon", "lat")

## plot on world country map, setting the xlim and ylim to span the range of longitudes and latitudes
plot(wrld_simpl, xlim = range(species.coords$lon), ylim = range(species.coords$lat), axes = TRUE, col = "light yellow")
## add points for this species
points(species.coords, col = "red", cex = 0.75)

## this includes the whole world - let's trim down to the data for Africa
## A bit of digging online suggests that latitudes +40 to -40 and longitudes -20 to +55
## First, here is a function to pull out just the data within a lat and lon range

trim.coords <- function(x, latmin, latmax, lonmin, lonmax) {
    x[x$lon >= lonmin & x$lon <= lonmax & x$lat >= latmin & x$lat <= latmax, ]
}

## Then use the function to make a new table of coordinates, just within the ranges specified
species.coords.trim <- trim.coords(species.coords, latmin = -40, latmax = 40, lonmin = -20, lonmax = 55)

## plot world map again now using the trimmed data
plot(wrld_simpl, xlim = range(species.coords.trim$lon), ylim = range(species.coords.trim$lat), axes = TRUE, col = "light yellow")
## add points for this species
points(species.coords.trim, col = "red", cex = 0.75)

## if it looks OK, let's use the trimmed version as our new coordinates
species.coords <- species.coords.trim

## some code for removing points that are in the sea, bit gnarly, just paste as a block and cross fingers
## this website has data specifying oceans


# Download and process ocean data
ocean_data_dir <- here("data", "raw", "ocean")
if (!dir.exists(ocean_data_dir)) dir.create(ocean_data_dir)

URL <- "https://naturalearth.s3.amazonaws.com/110m_physical/ne_110m_ocean.zip"
zip_file <- file.path(ocean_data_dir, basename(URL))

# Download only if file doesn't exist already
if (!file.exists(zip_file)) {
    download.file(URL, zip_file)
}

# Unzip to ocean data directory
fils <- unzip(zip_file, exdir = ocean_data_dir)

# Read shapefile
oceans <- read_sf(grep("shp$", fils, value = TRUE))

## convert our locality points into equivalent format
species.coords <- st_as_sf(species.coords, coords = c("lon", "lat"))
st_crs(species.coords) <- st_crs(oceans)
sf_use_s2(FALSE)  # Disable spherical geometry

## find where out points intersect with the ocean
tmp <- sapply(st_intersects(species.coords, oceans), function(z) if (length(z) == 0) NA_integer_ else z[1])

## remove points that intersect with the ocean and convert back to table of coordinates
if (sum(!is.na(tmp)) > 0) {
    species.coords <- data.frame(st_coordinates(species.coords[is.na(tmp), ]))
} else {
    species.coords <- data.frame(st_coordinates(species.coords))
}
colnames(species.coords) <- c("lon", "lat")

## now plot again to check
plot(wrld_simpl, xlim = range(species.coords$lon), ylim = range(species.coords$lat), axes = TRUE, col = "light yellow")
points(species.coords, col = "red", cex = 0.75)


#########################################################################
## 3) Extract climatic values for the localities occupied by the species##
#########################################################################

## download bioclimatic data from the worldclim database and convert to Raster format
bio.data <- worldclim_global(var = "bio", res = 10, path = here("data", "raw"))
names(bio.data) <- paste0("bio", 1:19)

## there are 15 bioclimatic variables stored in the bio.data object
## which are each defined here https://www.worldclim.org/data/bioclim.html
## e.g. bio1 is Annual Mean Temperature, bio2 is mean diurnal range
## They stored in the object as a stack of rasters
## A raster is  an image made up of pixels, in this case representing 10 minutes of a degree grid cells
## and the values of each pixel are the values for a bioclimatic variable
## you can plot them by specifying the number for each raster layer. I.e. to plot bio1 and bio2
plot(bio.data, 1, main="Annual Mean Temperature")
plot(bio.data, 2, main="Mean Diurnal Range")

## you can extract bioclimatic data for the focal localities in Africa where your species is found
bio.values <- extract(bio.data, species.coords)[, -1]
rownames(bio.values) <- rownames(species.coords)

## and plot them to see the range of environmental values the species lives in
plot(bio.values[, 1], bio.values[, 2])

## or look at how different variables correlate - here focusing on the first five
pairs(bio.values[, 1:5])

## append to lat long, remove rows with missing data, and save to file
species.data <- cbind(species.coords, bio.values)
write.csv(species.data, file = here("data", "processed", "species_data.csv"), row.names = FALSE)

## extract mean values, min and max
rbind(
    mean = colMeans(species.data),
    min = apply(species.data, 2, min),
    max = apply(species.data, 2, max)
)


############################################################################
## 4) Model the species current distribution based on bioclimatic variables##
############################################################################

## We will run a linear model - like a regression, but where Y is presence/absence
## and X variables are values of the bioclimatic variable.
## The problem is that we don't have absence data, i.e. that someone visited and said
## the species is definitely not present. So instead we create pseudo absence data
## using background points with no presence recorded from the same general region.
## I'll refer to these as background points

## Question 1: What are the limitations of this approach?

## a) Make random background points for your region
## This bit makes a raster from the world map, which is a shaded region
ext <- extent(wrld_simpl)
xy <- abs(apply(as.matrix(bbox(ext)), 1, diff))
n <- 5
r <- raster(ext, ncol = xy[1] * n, nrow = xy[2] * n)
mask <- rasterize(wrld_simpl, r)

## This plots it and shows your points
## I'm setting a plotting window using the same lat and lon as above
plot(mask, axes = TRUE, col = "grey", xlim = c(-20, 55), ylim = c(-40, 40), legend = FALSE)
points(species.coords, col = "red", pch = 20, cex = 1)

## It looks a bit too big still, lots of places with no presences, so I try a smaller view
plot(mask, axes = TRUE, col = "grey", xlim = c(-15, 50), ylim = c(-40, 20), legend = FALSE)
points(species.coords, col = "red", pch = 20, cex = 1)

## That looks better.

## Now I can set an area of extent for considering background using this new range
e <- extent(-15, 50, -40, 20)

## this generates 500 random points within that area, and plots them
bg <- randomPoints(mask, 500, ext = e, extf = 1)
colnames(bg) <- c("lon", "lat")

plot(mask, axes = TRUE, col = "grey", xlim = c(-15, 50), ylim = c(-40, 20), legend = FALSE)
plot(e, add = TRUE, col = "black")
points(species.coords, col = "red", pch = 20, cex = 1)
points(bg, col = "black", pch = 20)

## the black points are the background points with no recorded presence in GBIF
## the red points are the observed localities

## Question 2: Why is it better not to have too large a region for your background?
## Question 3: Do you think C. arabica is really absent from all the
## localities selected as background points? Will it matter?

## Now we are settled on our area of extent, we can crop the bio.data to just keep values for this region
## It isn't essential but speeds up some later steps
bio.data <- crop(bio.data, e)


## b) Next step, combine the presence data and the background data in one table
train <- rbind(species.coords, bg)
## and make a vector record 1=presence, 0=background
pb_train <- c(rep(1, nrow(species.coords)), rep(0, nrow(bg)))
## extract bioclimatic data for these points
envtrain <- extract(bio.data, train)
envtrain <- data.frame(cbind(pa = pb_train, envtrain))
## and for each set separately
testpres <- data.frame(extract(bio.data, species.coords))
testbackg <- data.frame(extract(bio.data, bg))


## c) Now we are ready to do a logistic regression to predict presence/absence
## 	 We use a general linear model assuming binomial errors, which is appropriate
## 	 for modelling a binary Y variable of 0s and 1s.
## 	 More on this in Year 2 lecture 90258 Analysis of associations Part III: Non-Linear Regression Bonsall, Michael, and in Year 4!

## 	 To start with I am including the first 5 bioclim variables as predictors, but you could include 1,2 or more.
gm1 <- glm(pa ~ bio1 + bio2 + bio3 + bio4 + bio5,
    family = binomial(link = "logit"), data = envtrain
)

## look at a summary of the results
summary(gm1)

## Question 4: Which variables contribute significantly to explaining presence/absence?

## d) Predict species distribution from the model and plot it
## Based on the model, we can predict the probability of occurrence of the species
## across the whole of the area being considered
pg <- predict(bio.data, gm1, ext = e, type = "response")
pg <- crop(pg, e)

## pg is a raster layer, like for our bioclim variables, but now representing the
## probability of occurrence from our linear model, for or area of extent e

## plot this
plot(pg, main = "GLM probability of occurrence")
## add country boundaries
plot(wrld_simpl, add = TRUE, border = "dark grey")
## add our observed locality data
points(species.coords, col = "black", pch = 4, cex = 0.5)

## Question 5: How well do you think the model has predicted the distribution?

## If you want to show a single map of distribution instead, can convert the
## probabilities by selecting a threshold that gives best match between predicted
## and observed localities
## first we evaluate how well the model predicts presence/absence at each point
ge <- evaluate(testpres, testbackg, gm1)
ge
## the output gives several metrics such as
## Area Under the Curve: AUC, and correlation coefficient between observed and predicted.
## Higher values of both metrics = better match between model predictions and observed
## presence/absences.

## then we use this evaluation to pick a threshold probability for defining presence/absence
## using the model that gives the most accurate match to observed presence/absence
tr <- threshold(ge, "prevalence")
plot(pg > tr, main = "presence/absence")
plot(wrld_simpl, add = TRUE, border = "dark grey")
points(species.coords, col = "black", pch = 4, cex = 0.5)

## e) You can construct alternative models and compare different models

## This model uses next 5 climatic variables instead
gm2 <- glm(pa ~ bio6 + bio7 + bio8 + bio9 + bio10, family = binomial(link = "logit"), data = envtrain)
summary(gm2)

## You can compare two models, even with different predictor variables, using the Akaike Information Criterion
AIC(gm1, gm2)

## AIC helps compare model fit vs complexity (lower = more parsimonious)
## While differences >2 units suggest meaningful differences between models,
## model selection should be guided by theory and research questions,
## not just by AIC or p-values. Consider multiple evaluation metrics and
## use cross-validation when possible, and give enough thought to what your aims are: do you want to 
## predict well, or understand the underlying biology?

## You can try fitting different models with different predictor variables, but
## remember that the more variables you test, the more likely you are to overfit the observed data, which can 
## lead to poor predictions on new data.

## You can also compare metrics of how well the model predicts the data
evaluate(testpres, testbackg, gm1)
evaluate(testpres, testbackg, gm2)

## Question 6: Which bioclimatic variables do you expect to be most important for the species?
## Does a model including these variables predict distribution better than alternatives?


#################################
## PREDICT FUTURE SPECIES RANGES##
#################################

## a) Download one set of future climate data for period 2061-2080
## 	 There are multiple future climate models, and versions of each, which are stored in CMIP5
## 	 This code just extracts one model and set of parameters for you to try out the methods

future.bio.data <- cmip6_world(model = "CanESM5", var = "bio", ssp = "245", res = 10, time = "2061-2080", path = getwd())
names(future.bio.data) <- names(bio.data)

## There are several different climate models you can consider, I just picked one arbitrary one

## Crop future climate data to region of interest for efficiency
future.bio.data <- crop(future.bio.data, e)

## Fit model with present data
gm1 <- glm(pa ~ bio1 + bio2 + bio3 + bio4 + bio5,
    family = binomial(link = "logit"), data = envtrain
)

## Calculate predictions for present and future
pg <- predict(bio.data, gm1, ext = e, type = "response")
pg <- crop(pg, e)
pg.future <- predict(future.bio.data, gm1, ext = e, type = "response")
pg.future <- crop(pg.future, e)

## Plot results side by side
par(mfrow = c(1, 2))

# Present distribution
plot(pg, main = "A) GLM present")
plot(wrld_simpl, add = TRUE, border = "dark grey")
points(species.coords, col = "black", pch = 4, cex = 0.5)

# Future distribution
plot(pg.future, main = "B) GLM, 2060-2081")
plot(wrld_simpl, add = TRUE, border = "dark grey")
points(species.coords, col = "black", pch = 4, cex = 0.5)

## Question 7: How is the distribution of the species expected to change in the future?

## Extract counts for range changes
## First, get predictions for present and future at all localities
predict.localities.now <- extract(pg >= tr, species.coords)[, -1]
predict.localities.future <- extract(pg.future >= tr, species.coords)[, -1]

## We can extract some numbers about how the range will change, for example:
present_range <- sum(predict.localities.now)
future_range <- sum(predict.localities.future, na.rm = TRUE)
range_expansion <- sum((predict.localities.now == 0) & (predict.localities.future == 1), na.rm = TRUE)
range_contraction <- sum((predict.localities.now == 1) & (predict.localities.future == 0), na.rm = TRUE)

## Print results in a clear format
print("Range change metrics:")
print(paste("Current suitable localities:", present_range))
print(paste("Future suitable localities:", future_range))
print(paste("Number of new suitable localities (expansion):", range_expansion))
print(paste("Number of lost suitable localities (contraction):", range_contraction))


## Question 8: Use these numbers to calculate the percentage contraction or expansion in the range.

## A final plot that can be useful to understand causes of changes is to plot how climate will change

## work out the change in bioclim variables from now to the future
change.bio.data <- future.bio.data - bio.data

## plot present, future and change in climate for bioclim 1 variable, mean annual temperature
par(mfrow = c(1, 3))
plot(bio.data[[1]], ext = e, main = "Present Day")
plot(future.bio.data[[1]], ext = e, main = "2061-2080")
plot(change.bio.data[[1]], ext = e, main = "Projected Change") 
