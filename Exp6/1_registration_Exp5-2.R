
# e15 Embryo Registrations in Parallel
# Authors: R.A. Roston & A.M. Maga
# Date: 2023-01-25

# Exp5 Heartx0.8

# Input working directory, specify atlas, specify threads
## The working directory should have two subfolders: 
## "Input" with the .nrrd & .mrk.json files and "Output" for the transform files

experiment = "Exp5"
dir.imgs <- "/home2/rachel/P01/KOMP_Analysis/Experiments/Exp5_HeartShrink/CT/"
dir.lms <- "/home2/rachel/P01/KOMP_Analysis/Experiments/baselines/LMs/"

setwd("/home2/rachel/P01/KOMP_Analysis/Experiments/Exp5_HeartShrink/")
dir.transforms <- paste0(getwd(), "/Transforms/")
dir.totaltransforms <- paste0(getwd(), "/TotalTransforms/")
dir.imgs.transformed <- paste0(getwd(), "/CT_transformed/")

dir.create(dir.transforms)
dir.create(dir.totaltransforms)
dir.create(dir.imgs.transformed)

ref <- "Embryo_Atlas" # atlas aka fixed image
Sys.setenv(ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS = 48)
ncluster <- 4 # set up how many jobs run in parallel
nthreads <- 12 # limit number of cores for each task


library(doParallel)
library(foreach)
library(SlicerMorphR)
library(patchMatchR)
library(ANTsR)


# FUNCTIONS

doRegistration = function(specimen, 
                          ref, 
                          dir.imgs, 
                          dir.lms, 
                          dir.transforms,
                          dir.totaltransforms,
                          dir.imgs.transformed) {
  

  library(SlicerMorphR)
  library(patchMatchR)
  library(ANTsR)
  
  # INPUT SPECIMEN FILES
  tmp.lm <- read.markups.json(paste0(dir.lms, specimen, "_baseline.mrk.json"))
  tmp.img <- antsImageRead(paste0(dir.imgs, specimen, "_Exp5.nrrd"))
  tmp.img <- iMath(tmp.img, "Normalize")
  
  # INPUT ATLAS / REFERENCE FILES and MOVING IMG FILES
  ref.lm <- read.markups.json(paste0(dir.lms, ref, ".mrk.json"))
  ref.img <- antsImageRead(paste0(dir.imgs, ref, ".nii.gz"))
  ref.img <- iMath(ref.img, operation = "Normalize") #normalizes intensities to 0 1 range
  
  # REGISTRATION 
  lm.tx <- fitTransformToPairedPoints(fixedPoints = ref.lm, 
                                      movingPoints = tmp.lm,
                                      transformType = "Similarity", 
                                      domainImage = ref.img)
  
  #lm.check = 10
  #if (lm.tx$error < lm.check) {
    
    affine <- antsRegistration(fixed=ref.img, moving = tmp.img, typeofTransform = "Rigid", initialTransform = lm.tx$transform)
    syn1 <- antsRegistration(fixed=ref.img, moving = tmp.img, typeofTransform = "SyN", initialTransform = affine$fwdtransforms)
    
    composite.transforms <- antsApplyTransforms(fixed = ref.img, 
                                                moving = tmp.img,
                                                transformlist = c(syn1$fwdtransforms[1], syn1$fwdtransforms[2]),
                                                compose = paste0(dir.totaltransforms, specimen, "_", experiment, "_Total-"))
    
    composite.img <- antsApplyTransforms(fixed = ref.img,
                                        moving = tmp.img,
                                        transformlist = composite.transforms)
    
    file.copy(from = syn1$fwdtransforms[1],to = paste0(dir.transforms, specimen, "_", experiment, "_Warp", ".nii.gz"))
    file.copy(from = syn1$fwdtransforms[2],to = paste0(dir.transforms, specimen, "_", experiment, "_Affine", ".mat"))
    file.copy(from = syn1$invtransforms[2],to = paste0(dir.transforms, specimen, "_", experiment, "_InverseWarp", ".nii.gz"))
    antsImageWrite(image = composite.img, filename = paste0(dir.imgs.transformed, specimen, "_", experiment, "_transformed.nrrd"))
    
#  } else print( paste0("check the landmarks for this sample: ", specimen))
  
}

registrationresults <- function(specimens){
  
  library(ANTsRCore)
  library(ANTsR)
  
  dir.create("CT_transformed_sag")
  for(i in 1:length(specimens)){
    jpeg(filename = paste0("CT_transformed_sag", "/", specimens[i], "_", experiment, ".jpg"))
    tmp.img.transformed <- antsImageRead(paste0(dir.imgs.transformed, specimens[i], "_", experiment, "_transformed.nrrd"))
    txt <- list(x = 6, 
                y = 2,
                label = paste0(specimens[i], "_", experiment),
                cex = 0.7,
                col = 'white')
    plot(tmp.img.transformed,
         axis = 1,
         slices = 170,
         text = txt)
    dev.off()
  }
  dev.off()
}

# CREAT VECTORS OF SPECIMEN NAMES AND IMAGE FILES TO ANALYZE
# Landmarks needs to be saved in JSON format from Slicer. 

f <- dir(path = dir.lms, patt="json") #get file names for all json files
f <- f[-grep(ref, f)]
done <- dir(path=dir.transforms, patt=".mat")
done <- gsub("Exp5_Affine.mat", "baseline.mrk.json", fixed=T, done ) #subs _lms.mrk.json for -label.nii.gz
f <- f[which(! f %in% done)] # ! means not f in done, so this subsets elements n f that are not in done
specimens <- gsub("baseline.mrk.json", "Exp5", fixed=TRUE, f) #make specimen list by subbing all ".mark.json" with "" and return values (not indices)
imgs <- paste0(specimens, ".nrrd") #make an image list by pasting ".nrrd" to specimens

remove(done) #removes the variable done from the environment
remove(f)

specimens = gsub("_Exp5", "", specimens)





# CHECK THAT LIST OF .nrrds MATCHES ACTUAL .nrrd FILE NAMES 
# because this is within the larger function, it might not actually print outside of the function. I need a way to deal with this
if (!all(file.exists(paste0(dir.lms, "/", specimens, "_baseline.mrk.json")))) {
  print("something's mismatched, check files")
} #checks if list of .nrrds from _lms.mrk.jason list matches the actual .nrrd file names




# DO PARALLEL REGISTRATIONS
# This will the registration 8 times (1:8) in batches of 4 (PSOCK cluster setting.)
# Set up how many jobs are going to run in parallel & working directory
cl <- makePSOCKcluster(ncluster)
Sys.setenv(ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS = nthreads)
registerDoParallel(cl)

foreach (i = 1:length(specimens[2:length(specimens)])) %dopar% {doRegistration(specimen = specimens[i],
                                                          ref = ref,
                                                          dir.imgs = dir.imgs,
                                                          dir.lms = dir.lms,
                                                          dir.transforms = dir.transforms,
                                                          dir.totaltransforms = dir.totaltransforms, 
                                                          dir.imgs.transformed = dir.imgs.transformed)}



stopCluster(cl)

# make sure you execute stopCluster(cl) line. Otherwise resources assigned to the job will not be released
# and you might get confused about how many jobs you are running. 


doRegistration(specimen = specimens,
               ref = ref,
               dir.imgs = dir.imgs,
               dir.lms = dir.lms,
               dir.transforms = dir.transforms,
               dir.totaltransforms = dir.totaltransforms, 
               dir.imgs.transformed = dir.imgs.transformed)
