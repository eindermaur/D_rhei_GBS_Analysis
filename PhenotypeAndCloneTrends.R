# This script will summarize trends from clone correction of D. rhei samples
# Output:
  # A table summarizing number of isolates and clones found across fields
  # A map of NY showing sample locations

#### Import data and packages ####

library(RColorBrewer)
library(sp)
library(raster)
library(lme4)

didy <- read.csv("Isolate_metadata.csv", na.strings = "NA", stringsAsFactors = F)
clones <- read.table("clone_assignments_both_years.txt", header = T, stringsAsFactors = F)
indvs.cc <- read.table("rhei_diversity_cc_bothyears.012.indv", stringsAsFactors = F)
indvs.cc <- unlist(indvs.cc$V1)

#Clean up all sample names --- get rid of extra fields output by TASSEL
shorten_name <- function(x){
  return(unlist(strsplit(x,":"))[1])
}
clones$Sample <- sapply(clones$Sample, shorten_name)
indvs.cc <- sapply(indvs.cc, shorten_name)

#### Isolate number table and plotting metadata ####

#Make table assigning counties to geographical regions of NY
regions <- rbind(c("Cattaraugus", "WNY"),
                 c("Erie", "WNY"),
                 c("Montgomery", "CD"),
                 c("Ontario", "CNY"))
regions <- as.data.frame(regions)
colnames(regions) <- c("County", "Region")

#Get 'populations' designation: combination of year and field
populations <- paste(didy$Field, " (", didy$year, ")", sep = "")

#Now make table assigning fields to counties, regions, and plotting colors/symbols
field_info <- data.frame("Field" = unique(populations),
                         "County" = NA,
                         "Region" = NA)
for(i in 1:nrow(field_info)){
  field <- as.character(field_info$Field[i])
  county <- NA
  if(field == "NY unknown"){
    region <- "NY unknown"
  }else if(field %in% c("WA", "CO")){
    region <- "Non NY"
  }else{
    county <- unlist(strsplit(field," #", fixed = T))[1]
    region <- as.character(regions$Region)[regions$County==county]
  }
  field_info$Region[i] <- region
  field_info$County[i] <- county
}

#Put in alphabetical order with non-NY locations at end
region_order <- c("CD", "CNY", "WNY")
field_unordered <- field_info
start <- 1
for(region in region_order){
  region.size <- sum(field_unordered$Region==region)
  field_info[start:(start+region.size-1),] <- field_unordered[field_unordered$Region==region,]
  start <- start + region.size
}

#Add color and pch
field_info$color <- NA
#color_types <- brewer.pal(7,"Spectral")[-3] #changed from Greg's code
color_types <- brewer.pal(7,"Accent")[-6] #changed from Greg's code
field_info$color <- color_types[match(field_info$Region, unique(field_info$Region))]

#order sites within fields and add pch info
field_info$pch <- NA
for(region in unique(field_info$Region)){
  n.sites <- sum(field_info$Region==region)
  field_info.region <- field_info[field_info$Region==region,]
  field_info[field_info$Region==region,] <- field_info.region[order(field_info.region$Field),]
  field_info$pch[field_info$Region==region] <- 1:n.sites
}

#Save table assigning isolates to fields, counties, regions, and plotting info
isolate_metadata <- cbind("SamplesSZ" = didy$SampleSZ,
                          field_info[match(populations, field_info$Field),])
write.csv(isolate_metadata, "isolate_plotting_metadata.csv", quote = F, row.names = F)
write.csv(field_info, "field_plotting_metadata.csv", quote = F, row.names = F)

#Get total isolates per field-year and total unique genotypes per field-year
populations[populations %in% c("WA", "CO")] <- "Non NY"
population_counts <- table(populations)
field_sums <- data.frame("Field" = gsub("(.*) \\(.*", "\\1", rownames(population_counts), perl = T),
                         "Year" = as.integer(gsub(".*\\((.*)\\)", "\\1", rownames(population_counts), perl = T)),
                         "NumberIsolates" = as.integer(population_counts))
field_counties <- sapply(as.character(field_sums$Field), function(x) unlist(strsplit(x, " "))[1])
field_sums$Region <- regions$Region[match(field_counties, regions$County)]

#### Sample Map ####

#Note this includes all isolates including those that didn't pass genotyping

usa <- raster::getData('GADM', country='USA', level=2)
ny <- usa[usa$NAME_1=="New York",]
ny.clean <- ny[ny@data$TYPE_2=='County',]
counties.unique <- unique(didy$County)
counties.unique <- counties.unique[!is.na(counties.unique)]
#pdf("plots/fieldmap.pdf", width=3.25, height = 2.141176) ##unsure of origin in Greg's files?
old.par <- par(no.readonly = T)
par(ann = FALSE,
    bg = "white",
    bty = "n",
    mai = c(0,0,0,0),
    mgp = c(0,0,0),
    oma = c(0,0,0,0),
    omd = c(1,1,1,1),
    omi = c(0,0,0,0),
    usr = c(-79.5, -71, 40, 45.6),
    pin = c(3.25,2.141176),
    plt = c(0,1,0,1),
    pty = "m",
    xaxs = 'i',
    xaxt = 'n',
    xpd = FALSE,
    yaxs = 'i',
    yaxt = 'n')
plot(ny.clean, xlim = c(-79.5, -71), ylim=c(40,45.6))
set.seed(823)
for(county in counties.unique){
  county.polygon <- ny[ny$NAME_2==county,]
  county.color <- unique(field_info$color[field_info$County==county])
  county.color <- county.color[!is.na(county.color)]
  for(field in unique(didy$Field[didy$County == county & !is.na(didy$County)])){
    if(field != "Other NY"){
      field.point <- spsample(county.polygon,1,type = 'random')
      for(year in unique(didy$year[didy$Field == field])){
        n.isolates <- sum(didy$year == year & didy$Field == field, na.rm = T)
        if(year < 2022){
          point.pch = 17
        }else{
          point.pch = 16
        }
        points(field.point, cex=((n.isolates-.95)^.20)+0.4,
               col = 'black',
               pch = point.pch)
        points(field.point, cex=(n.isolates-.95)^.20,
               col = county.color,
               pch = point.pch)
      }
    }
  }
}

legend(-72.8, 45.35,
       pt.cex = 1.5,
       pch = c(17,16),
       legend = c("2021", "2022"),
       cex=1.5, bty='n', y.intersp = 1.5)
legend(-72.8, 44.2,
       legend = unique(field_info$Region)[1:3],
       fill = unique(field_info$color)[1:3],
       cex=1.5, bty='n', y.intersp = 1.5)
par(old.par)
dev.off()

#### Random Calculations ####

#How many clonal lineages are there?
ug_sizes <- table(didy$UniqueGenotype)
clineages <- as.integer(names(ug_sizes[ug_sizes>1]))
length(clineages)
cl_sizes <- ug_sizes[names(ug_sizes) %in% clineages]
sum(cl_sizes)
max(cl_sizes)
median(cl_sizes)
