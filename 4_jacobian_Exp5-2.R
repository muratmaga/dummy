# Jacobian Analysis
# Author: R.A. Roston, Ph.D., A.M. Maga, Ph.D.
# Date: 2023-04-17

# function to save jacobian determinant images
# function to do jacobian analysis, producing heatmap plots in a pdf doc
# current version does not downsample transforms, mask, or atlas (this option will be added later)



saveJacobian <- function(mask, 
                         mouseline,
                         atlas.file = "/home2/rachel/P01/KOMP_Analysis/CT/Embryo_Atlas.nii.gz",
                         dir.transforms = "/home2/rachel/P01/KOMP_Analysis/AAA/Transforms", 
                         dir.mask = "/home2/rachel/P01/KOMP_Analysis/AAA/Masks/",
                         dir.output = "/home2/rachel/P01/KOMP_Analysis/AAA/JacobianAnalysis/"
                         ){
  
    library(ANTsR)
    library(ANTsRCore)
  
    setwd(dir.transforms)
    files = dir (patt = mouseline)
    files = files[grep('_Warp.nii.gz', files, fixed=TRUE)]
    
    warps = jacs= list()
    for (i in 1:length(files)) warps[[i]] = antsImageRead(files[i]) 

    atlas = antsImageRead(atlas.file)
    atlas_mask = antsImageRead(paste0(dir.mask, "Embryo_Atlas_mask_", mask, ".nrrd"))
    # atlas_mask = thresholdImage(atlas_mask, 1, 103) #merges all the labels together. You should build your own masks, or use the ones from PCA. 
    #Statistics will be only evaluated within the mask

  for (i in 1:length(files)){
    tmp = createJacobianDeterminantImage( atlas, warps[[i]], doLog = TRUE, geom = TRUE )
    jacs[[i]] = smoothImage( tmp, 2*antsGetSpacing(atlas)[1], FALSE ) #smooth the edges.  
  }
    
    
# save jacobians
  jacs.files = gsub(files, pattern ='Warp.nii.gz', replacement = paste0("Jacobian_", mask, ".nrrd"))
  setwd(dir.output)
  
  dir.create(paste0(dir.output, mouseline, "_", mask, "/"))
    
  for (i in 1:length(jacs)){
    antsImageWrite(image = jacs[[i]], 
                   filename = paste0(dir.output, mouseline, "_", mask, "/", jacs.files[i]))
  }  
  
}


jacobian.analysis <- function(mask,
                              ko.line,
                              run,
                              P = 0.05,
                              ctrl = "baseline",
                              atlas.file = "/home2/rachel/P01/KOMP_Analysis/CT/Embryo_Atlas.nii.gz",
                              dir.transforms = "/home2/rachel/P01/KOMP_Analysis/Experiments/Exp5_HeartShrink/Transforms/", 
                              dir.mask = "/home2/rachel/P01/KOMP_Analysis/Masks/",
                              dir.output = "/home2/rachel/P01/KOMP_Analysis/Experiments/Exp5_HeartShrink/Jacobian/"){

  library(ANTsR)
  library(ANTsRCore)
  library(viridis)
  
  
  
  calc.jacs = function(files, mask1 = mask, atlas.file1 = atlas.file, dir.transforms1 = dir.transforms, dir.mask1 = dir.mask){
    warps = jacs= list()
    for (i in 1:length(files)) warps[[i]] = antsImageRead(paste0(dir.transforms1, files[i]))
    
    atlas = antsImageRead(atlas.file1)
    
    # atlas_mask = thresholdImage(atlas_mask, 1, 103) #merges all the labels together. You should build your own masks, or use the ones from PCA. 
    #Statistics will be only evaluated within the mask
    
    for (i in 1:length(files)){
      tmp = createJacobianDeterminantImage( atlas, warps[[i]], doLog = TRUE, geom = TRUE )
      jacs[[i]] = smoothImage( tmp, 2*antsGetSpacing(atlas)[1], FALSE )
    }
    
    return(jacs)
  }
  
  setwd(dir.transforms)
  muts = dir (patt = ko.line)
  muts = muts[grep('_Warp.nii.gz', muts, fixed=TRUE)]
  
  if(mask == "wholebody"){
    ctrls = dir(paste0(dir.output, "baseline_wholebody/"))
    print("File names obtained")
    
    print("Loading saved Jacobian determinants for ctrls")
    jacs.ctrls = list()
    for(i in 1:length(ctrls)) jacs.ctrls[[i]] = antsImageRead(paste0(dir.output, "baseline_wholebody/", ctrls[i]))
    
    print("Calculating Jacobian determinants for muts")
    jacs.muts = calc.jacs(muts)
    
    jacs = c(jacs.ctrls, jacs.muts)
    
  } else {
    ctrls = dir (patt = ctrl)
    ctrls = ctrls[grep('_Warp.nii.gz', ctrls, fixed=TRUE)]
    
    files = c(ctrls, muts)
    print("File names obtained")
    
    print("Calculating Jacobian determinants for ctrls and muts")
    jacs = calc.jacs(files)
    
    }
  
  groups = as.factor(c(rep("C", length(ctrls)), rep("M", length(muts))))
 
  
  #set the regression
  print("Starting statistical analysis")
  atlas_mask = antsImageRead(paste0(dir.mask, "Embryo_Atlas_mask_", mask, ".nrrd"))
  j_mat = imageListToMatrix( jacs, atlas_mask )
  j_mdl = lm( j_mat ~ groups)    # raw volume data
  j_bmdl = bigLMStats( j_mdl , 1.e-5 )
  print( min( p.adjust( j_bmdl$beta.pval, 'none' ) , na.rm=T ) )
  
  #do multiple corrections and prepare the plots
  #P defines the significance level we want to evaluate things after FDR correction
  volstats = j_bmdl$beta.pval
  correct.ps= p.adjust(volstats,'fdr')
  pImg = makeImage( atlas_mask, correct.ps)
  pMask = thresholdImage(pImg, 10^-16, P, 1, 0)
  betaImg = makeImage( atlas_mask, j_bmdl$beta ) * pMask
  
  # plots will display only the voxels (within the mask) that remain significantly different after the FDR 
  # the color map is the magnitude of the effect 
  # positive values expansion, negative values contraction with respect to the atlas.
  
  print("Creating images")
  atlas = antsImageRead(atlas.file)
  plot.folder = dir.output
  
  jpeg(file = paste0(plot.folder, ctrl, "vs", ko.line, "_", run, "_pos.jpg"), quality = 100)
  plot( atlas, betaImg, axis=1, 
        nslices=30, ncolumns=5, 
        window.overlay = c(0.001,1), 
        color.overlay = "viridis",
        alpha=.8, colorbar = T)
  dev.off()
  
  jpeg(file = paste0(plot.folder, ctrl, "vs", ko.line, "_", run, "_neg.jpg"), quality = 100)
  plot( atlas, betaImg, axis=1, 
        nslices=30, ncolumns=5, 
        window.overlay=c(-1, -.001), 
        color.overlay = "viridis",
        alpha=.8, colorbar = T)
  dev.off()
  
  print("Jacobian analysis completed")
}


jacobian.analysis(mask = "wholebody",
                  ko.line = )
