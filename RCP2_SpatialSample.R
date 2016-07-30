
require(rgdal)
require( sp)
require( randtoolbox)
require( class)
require( fields)
require( maptools)

#state boundaries
AusBorders <- readShapeLines("./aust_states/australia_states.shp")
NSWborder <- AusBorders[AusBorders$STATE %in% c("NSW"),]

#point data (sampling locations)
data <- read.csv("NSW_FF_locations_breif.csv")
data$uniqueID <- 1:nrow( data)
#unique ID needed as the text id has duplicates (199+ of them):
sum( duplicated( data$SiteNo))

plot( NSWborder)
points( data[,3:2], pch='.')

#there must be a better way to do this!
#dropped ACT border and Island
NSWcoords <- coordinates( NSWborder@lines[[1]]@Lines[[2]])
NSWboundingBox <- matrix( c( range( NSWcoords[,1]), range( NSWcoords[,2])), ncol=2)

plot( NSWcoords, type='l')
points( data[,3:2], pch='.')

#crop data outside NSW
ids <- point.in.polygon( data[,3], data[,2], NSWcoords[,1], NSWcoords[,2])
data <- data[ids==1,]

my.seed <- 5
invisible( halton( sample( 1:1000, size=1), dim=2, init=TRUE))	#initialise the sequence -- perverse way of doing it.
#the quasi-random (spatially balanced) set of points
qn <- halton( 2*nrow( data), dim=2, init=FALSE)	#yes a million is a lot, but they will need to be subsetted yet.
qn[,1] <- NSWboundingBox[1,1] + (NSWboundingBox[2,1]-NSWboundingBox[1,1]) * qn[,1]
qn[,2] <- NSWboundingBox[1,2] + (NSWboundingBox[2,2]-NSWboundingBox[1,2]) * qn[,2]

ids <- point.in.polygon( qn[,1], qn[,2], NSWcoords[,1], NSWcoords[,2])

plot( NSWcoords, type='l')
points( qn[ids==1,], col=3, pch=20)
points( qn[ids==0,], col=2, pch=20)

#crop sample points outside NSW
qn <- qn[ids==1,]	#just those points in NSW
qn <- qn[1:nrow( data),]	#a spatially balanced set

#finding nearest point to each qn location (and removing taken ones from subset)
spatsamp <- data.frame( siteID=rep( NA, nrow( data)), dist=rep( NA, nrow( data)))	#for the row numbers of data and the distance to sampled point
#loop takes a while, but only has to be done once
for( ii in 1:nrow( spatsamp)){
	cat( ii," ")
	trainingData <- data[!as.character( data[,4]) %in% spatsamp$siteID,]
	spatsamp$siteID[ii] <- trainingData[knn1( train=trainingData[,3:2], test=qn[ii,], cl=1:nrow( trainingData)),4]
	spatsamp$dist[ii] <- rdist.earth( data[data[,4]==spatsamp$siteID[ii],3:2], matrix( qn[ii,], nrow=1), miles=FALSE)	#kms
}
cat( "\n")
spatsamp$ord <- 1:nrow( spatsamp)

#so spatial sampling becomes less and less efficient
plot( 1:nrow( spatsamp), spatsamp[,2], pch='.', main="Distance from Spat. Sample to Observed")

#match up (well order) the original samples
tmpdata <- merge( data, spatsamp, by.x="uniqueID", by.y="siteID", sort=FALSE, all=TRUE)
tmpdata <- tmpdata[order( tmpdata$ord),]

plot( NSWcoords, type='l')
#still a little screwy out west
points( tmpdata[1:2000,4:3], pch=20)
#add sequentially
points( tmpdata[1:10000,4:3], pch=20)
points( tmpdata[1:20000,4:3], pch=20)
points( tmpdata[1:30000,4:3], pch=20)
points( tmpdata[1:40000,4:3], pch=20)
points( tmpdata[1:50000,4:3], pch=20)
points( tmpdata[1:60000,4:3], pch=20)
points( tmpdata[1:70000,4:3], pch=20)
points( tmpdata[1:80000,4:3], pch=20)

write.csv(tmpdata, file="spatialSample1.csv")



