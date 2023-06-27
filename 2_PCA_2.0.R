
PCA.analysis <- function(runID, # string
                   transform.type, # options: "Total", "Syn1"
                   mask, # options: "head", "limbs", "torso", "wholebody", "thorax", "abdomen"
                   mouselines, # vector containing mouselines to include in analysis
                   wd = "/home2/rachel/P01/KOMP_Analysis/AAA/",
                   dir.masks = paste0(wd, "Masks/"),
                   dir.output = paste0(wd, "PCA/")){
  
  # LIBRARIES
  library(ANTsR)
  library(ANTsRCore)                   
                     
  # FUNCTIONS
  print("Writing functions")
  doPCA <- function(dir.PCAresults, transformlist, dir.transforms, pca_mask){
    
    library(ANTsR)
    library(ANTsRCore)
    
    # Read transform files file transformlist and downsample
    print("Reading transforms into tlist")
    tlist <- list()
    for (i in 1:length(transformlist)) {
      tlist[[i]] <- resampleImage(antsImageRead(paste0(dir.transforms, transformlist[i])), c(0.054,0.054,0.054))
    }
    
    # Run PCA
    print("PCA")
    pca <- multichannelPCA(x = tlist, mask = pca_mask, pcaOption = "randPCA")
    
    # SAVE PCA RESULTS
    # save() does not work for ANTs images, so pcaWarps has to be saved separately
    "Saving PCA results"
    dir.create(dir.PCAresults)
    dir.create(paste0(dir.PCAresults, "/pcaWarps"))
    
    save(objects = transformlist, file = paste0(dir.PCAresults, "/transformlist"))
    print(paste0(dir.PCAresults, "/transformlist"))
    print(transformlist)
    save(objects = pca, file = paste0(dir.PCAresults, "/PCA"))
    for(i in 1:length(pca$pcaWarps)){
      antsImageWrite(image = pca$pcaWarps[[i]], filename = paste0(dir.PCAresults, "/pcaWarps/pcaWarp_", i, ".nrrd"))

    }
  }
  
  
  # Create transformlist
  print("creating transformlist")
  if(transform.type == "Total"){
    dir.transforms = paste0(wd, "TotalTransforms/")
    transformlist = vector()
    for(i in 1:length(mouselines)){
      tmp = dir(dir.transforms, pattern = mouselines[i])
      transformlist = c(transformlist, tmp)
    }
  } else if (transform.type == "Syn1") {
    dir.transforms = paste0(wd, "Transforms/")
    transformlist = vector()
    for(i in 1:length(mouselines)){
      tmp = dir(dir.transforms, pattern = mouselines[i])
      transformlist = c(transformlist, tmp)
    }
    transformlist = transformlist[grep(transformlist, pattern = "_Warp")]
  } else {print("Error: transform.type does not exist")}
  
  # Create output folders
  if (dir.exists(dir.output) == FALSE){
    print("Creating dir.PCAresults folder")
    dir.create(path = dir.output)
    dir.create(path = paste0(dir.output, runID))
    dir.PCAresults = paste0(dir.output, runID, "/", transform.type, "-", mask)
  } else { 
    print("Getting dir.PCAresults folder")
    dir.create(path = paste0(dir.output, runID))
    dir.PCAresults = paste0(dir.output, runID, "/", transform.type, "-", mask)
  }
  
  # Read Mask & do PCA
  print("Reading mask file")
  maskfile = paste0(dir.masks, "Embryo_Atlas_mask_", mask, ".nrrd")
  
  if( file.exists(maskfile)){
    
    pca_mask <- resampleImage(antsImageRead(maskfile), c(0.054,0.054,0.054), interpType = 1)
    
    print("Doing PCA")
    doPCA(dir.PCAresults = dir.PCAresults, 
          transformlist = transformlist, 
          dir.transforms = dir.transforms,
          pca_mask = pca_mask )
    print(dir.PCAresults)
    
    print("PCA analysis completed")
    
  } else {print("Error: Mask does not exist")}
  
}

# to do one PCA
start.time = Sys.time()
PCA.analysis(runID = "Psph-WT", 
             transform.type = "Syn1",
             mask = "wholebody",
             mouselines = c("baseline", "Psph"))
end.time = Sys.time()
print(end.time - start.time)

# to do multiple PCAs: same mouselines, diff masks and transforms
## when using doparallel / foreach, messages do not print out

ncluster = 10
nthreads = 4

runID = "Psph-WT"
mouselines = c("baseline", "Psph")
transform.type = c("Total", "Syn1")
masks = c("wholebody", "head", "torso", "thorax", "abdomen")

library(doParallel)
library(foreach)

cl <- makePSOCKcluster(ncluster)
Sys.setenv(ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS = nthreads)
registerDoParallel(cl)

for (i in 1:length(transform.type)) {
  foreach (j = 1:length(masks)) %dopar%  {
    print(paste0("Starting ", transform.type[i], "-", masks[j], " : ", Sys.time()))
    PCA.analysis(runID = runID,
                 transform.type = transform.type[i], 
                 mask = masks[j], 
                 mouselines = mouselines)
    print(paste0("Completed ", transform.type[i], "-", masks[j], " : ", Sys.time()))
  }
}

stopCluster(cl)