library(sf)
library(ggplot2)
library(viridis)
library(lwgeom)
library(raster)
library(INLA)
library(reshape2)

shapefile_path <- "D:/biostat article/st.louis project/LimaBorder/LimaBorder.shp"

# CRS
peru_sf <- st_read(shapefile_path) %>% 
  st_transform(crs = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs")

d <- read.csv('D:/biostat article/st.louis project/Peru_ModelingData_2017.csv')

d <- d[, c("Date", "GridID", "Lon", "Lat", "mean_PM25")]
names(d) <- c("date", "id", "long", "lat", "PM25value")

d$date <- as.Date(d$date)
d$date_num <- as.numeric(d$date)

p <- st_as_sf(d, coords = c("long", "lat"), crs = 4326) %>%
  st_transform(crs = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs")

d <- cbind(d, st_coordinates(p))

ind <- st_intersects(peru_sf, p, sparse = FALSE)

d <- d[rowSums(ind) > 0, ]

# Plot the map with points and zoom in
ggplot(peru_sf) + 
  geom_sf() + 
  coord_sf(datum = st_crs(peru_sf)) +
  geom_point(data = d, aes(x = X, y = Y), color = "black") + 
  theme_bw() 
# Faceted plot with zoom
ggplot(peru_sf) + 
  geom_sf() + 
  coord_sf(datum = NA) +
  geom_point(
    data = d, aes(x = X, y = Y, color = PM25value),
    size = 2
  ) +
  labs(x = "", y = "") +
  scale_color_viridis() +
  facet_wrap(~date) +
  theme_bw() 

# +
#  xlim(c(-1.075e+07, -1.04e+07)) +
#  ylim(c(-4.35e+06, -4.12e+06))

# INLA part
coo <- cbind(d$X, d$Y)
bnd <- inla.nonconvex.hull(coo, convex = -0.05, concave = 0.5)
mesh <- inla.mesh.2d(
  loc = coo, boundary = bnd,
  max.edge = c(100000, 200000), cutoff = 1000
)
mesh$n

plot(mesh)
points(coo, col = "red")

spde <- inla.spde2.pcmatern(
  mesh = mesh, alpha = 2, constr = TRUE,
  prior.range = c(10000, 0.01), # P(range < 10000) = 0.01
  prior.sigma = c(3, 0.01) # P(sigma > 3) = 0.01
)

timesn <- length(unique(d$date))
indexs <- inla.spde.make.index("s",
                               n.spde = spde$n.spde,
                               n.group = timesn
)
lengths(indexs)


group <- d$date_num - min(d$date_num) + 1
A <- inla.spde.make.A(mesh = mesh, loc = coo, group = group)

# prediction points
bb <- st_bbox(peru_sf)
x <- seq(bb$xmin, bb$xmax, length.out = 50)
y <- seq(bb$ymin, bb$ymax, length.out = 50)
dp <- as.matrix(expand.grid(x, y))
plot(dp, asp = 1)

p <- st_as_sf(data.frame(x = dp[, 1], y = dp[, 2]),
              coords = c("x", "y")
)
st_crs(p) <- st_crs("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs")
ind <- st_intersects(peru_sf, p)
dp <- dp[unlist(ind), ]
plot(dp, asp = 1)

dp <- do.call(rbind, lapply(1:timesn, function(i) cbind(dp, i)))
head(dp)

coop <- dp[, 1:2]
groupp <- dp[, 3]
Ap <- inla.spde.make.A(mesh = mesh, loc = coop, group = groupp)

print("Dimensions of A matrix:")
print(dim(A))
print("Dimensions of Ap matrix:")
print(dim(Ap))

print("Number of rows in d (estimation data):")
print(nrow(d))
print("Number of rows in dp (prediction data):")
print(nrow(dp))

print("Structure of indexs:")
str(indexs)

print("Number of rows in indexs$s:")
print(length(indexs$s))





effects_est <- list(data.frame(b0 = rep(1, nrow(d))), s = indexs$s[1:ncol(A)])
effects_pred <- list(data.frame(b0 = rep(1, nrow(dp))), s = indexs$s[1:ncol(Ap)])

stk.e <- inla.stack(
  tag = "est",
  data = list(y = d$PM25value),
  A = list(1, A),
  effects = effects_est
)


print("Dimensions of stk.e:")
str(stk.e)

stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(1, Ap),
  effects = effects_pred
)

stk.full <- inla.stack(stk.e, stk.p)

print("Indexs:")
print(indexs)

rprior <- list(theta = list(prior = "pccor1", param = c(0, 0.9)))

formula <- y ~ 0 + b0 + f(s, model = spde, control.group = list(model = "ar1", hyper = rprior))

# Fit the model
res <- inla(formula,
            data = inla.stack.data(stk.full),
            control.predictor = list(
              compute = TRUE,
              A = inla.stack.A(stk.full)
            )
)
summary(res)




index <- inla.stack.index(stack = stk.full, tag = "pred")$data

dp <- data.frame(dp)
names(dp) <- c("x", "y", "time")

dp$pred_mean <- res$summary.fitted.values[index, "mean"]
dp$pred_ll <- res$summary.fitted.values[index, "0.025quant"]
dp$pred_ul <- res$summary.fitted.values[index, "0.975quant"]


specific_time_points <- c(1,50, 164)

# Filter the data based on specific time points
dp <- dp[dp$time %in% specific_time_points, ]

dpm <- melt(dp,
            id.vars = c("x", "y", "time"),
            measure.vars = c("pred_mean", "pred_ll", "pred_ul")
)

ggplot(peru_sf) + 
  geom_sf() + 
  coord_sf(datum = NA) +
  geom_tile(data = dpm, aes(x = x, y = y, fill = value)) +
  labs(x = "", y = "") +
  facet_wrap(variable ~ time) +
  scale_fill_viridis("PM2.5") +
  theme_bw() 