
setwd("/home2/rachel/P01/KOMP_Analysis/Experiments/Exp5_HeartShrink/")
dir.baselineCT = "/home2/rachel/P01/KOMP_Analysis/Experiments/baselines/CT/"
dir.baselineTransforms = "/home2/rachel/P01/KOMP_Analysis/Experiments/baselines/Transforms/"
dir.baselineCTtransformed = "/home2/rachel/P01/KOMP_Analysis/Experiments/baselines/CT_transformed/"
atlas = antsImageRead("/home2/rachel/P01/KOMP_Analysis/Experiments/baselines/CT/Embryo_Atlas.nii.gz")
experiment = "Exp5" # Heartx0.8

library(ANTsR)
library(ANTsRCore)
library(SlicerMorphR)

# Create experimental deformation from original LM set and modified LM set using bSpline
originalLMs = read.markups.json("Exp5_deformation/mergedMarkups_Heartx1.mrk.json")
modLMs = read.markups.json("Exp5_deformation/mergedMarkups_Heartx0.8.mrk.json")


ex.tx = fitTransformToPairedPoints(movingPoints = originalLMs, 
                                   fixedPoints = modLMs, 
                                   transformType = "bspline", 
                                   domainImage = atlas, 
                                   meshSize = 4, 
                                   numberOfFittingLevels = 6)

deformedatlas = applyAntsrTransformToImage(transform = ex.tx$transform, image = atlas, reference = atlas, interpolation = "bspline")
antsImageWrite(deformedatlas, filename = paste0("Embryo_Atlas_", experiment, "-x0.8.nrrd"))
# Now check that atlas is deformed in expected way in Slicer (it is)


# Apply Exp5 deformation to subjects in atlas space
subjects = dir(dir.baselineCTtransformed, pattern = "baseline")
length(subjects)

landmarksdone = dir(path = "/home2/rachel/P01/KOMP_Analysis/Experiments/baselines/LMs/")
subjects = gsub("_transformed.nrrd", "", subjects)
landmarksdone = gsub(".mrk.json","", landmarksdone)
setdiff(landmarksdone, subjects) # do registration on missing specimens


subjects = dir(dir.baselineCTtransformed, pattern = "baseline")
length(subjects)

imgs.ex.tx = list()
for(i in 1:length(subjects)) {
  tmp.img = antsImageRead(paste0(dir.baselineCTtransformed, subjects[i]))
  imgs.ex.tx[[i]] = applyAntsrTransformToImage(transform = ex.tx$transform,
                                     image = tmp.img, 
                                     reference = deformedatlas,
                                     interpolation = "linear")
}

# Transform subjects back to subject space
subjects = gsub("baseline_transformed.nrrd", "", subjects)

invtransforms = vector()
for(i in 1:length(subjects)){
  invtransforms[i] = paste0(dir.baselineTransforms, subjects[i], "baseline_InverseWarp.nii.gz")
}

affine = vector()
for(i in 1:length(subjects)){
  affine[i] = paste0(dir.baselineTransforms, subjects[i], "baseline_Affine.mat")
}


for(i in 1:length(subjects)){
  tmp = antsApplyTransforms(fixed = antsImageRead(paste0(dir.baselineCT, subjects[i], "baseline.nrrd")),
                            moving = imgs.ex.tx[[i]], 
                            transformlist = c(affine[i], invtransforms[i]),
                            "linear")
  antsImageWrite(tmp, filename = paste0("CT/", subjects[i], experiment , ".nrrd"))
}


### Correct intensity range for MEMOS
subjects = dir("CT")
dir.create("CT_255")
dir.create(paste0("CT_255/Exp5"))
for(i in 1:length(subjects)){
  i1=antsImageRead(paste0("CT/", subjects[i]))
  i1=i1*255
  antsImageWrite(antsImageClone(i1, "unsigned char"), paste0("CT_255/Exp5/", subjects[i]))
}

dir.create(paste0("CT_255/baselines"))
subjects = gsub("Exp5.nrrd", "", subjects)
for(i in 1:length(subjects)){
  if(file.exists(paste0(dir.baselineCT, subjects[i], "baseline.nrrd"))){
    file.copy(from = paste0(dir.baselineCT, subjects[i], "baseline.nrrd"),
              to = paste0("CT_255/baseline/", subjects[i], "baseline.nrrd"))
   } else {print("something's wrong")}
}
# Now run MEMOS on CT_255/ data

MEMOS.dir = "/home2/rachel/P01/KOMP_Analysis/Experiments/baselines/MEMOS/"
MEMOS.baselines = dir(MEMOS.dir)

for(i in 1:length(MEMOS.baselines)){
file.rename(from = paste0(MEMOS.dir, MEMOS.baselines[i]),
            to = paste0(MEMOS.dir, gsub("-e15.seg.nrrd", "-e15.5_baseline.seg.nrrd", MEMOS.baselines[i])))
}

MEMOS.dir = "/home2/rachel/P01/KOMP_Analysis/Experiments/Exp5_HeartShrink/MEMOS/"
MEMOS.experiment = dir(MEMOS.dir)
for(i in 1:length(MEMOS.experiment)){
  file.rename(from = paste0(MEMOS.dir, MEMOS.experiment[i]),
              to = paste0(MEMOS.dir, gsub("-e15.seg.nrrd", "-e15.5_Exp5.seg.nrrd", MEMOS.experiment[i])))
}   

MEMOS.base.dir = "/home2/rachel/P01/KOMP_Analysis/Experiments/baselines/MEMOS/"
MEMOS.baselines = dir(MEMOS.base.dir)
MEMOS.exp.dir = "/home2/rachel/P01/KOMP_Analysis/Experiments/Exp5_HeartShrink/MEMOS/"
MEMOS.experiment = dir(MEMOS.exp.dir)

length(setdiff(gsub("_baseline.seg.nrrd", "", MEMOS.baselines), gsub("_Exp5.seg.nrrd", "", MEMOS.experiment)))
       