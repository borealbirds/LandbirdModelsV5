

#17. Package and save----
out <- data.table::rbindlist(out.list, fill=TRUE)

write.csv(out, file.path(root, "output", "validation", "ModelValidation_BCR.csv"), row.names = FALSE)

#TO DO: FINAL PREDICTION RASTER STACKS
#REPROJECT TO LAEA HERE

#TO DO: AVERAGE THE VALIDATION METRICS####

#TO DO: ADD ATTRIBUTES FOR SP ETC

#TO DO: CHECK EVERYTHING IS THERE#######

#LIST OF PRODUCTS PER SPECIES:
#2 extents: boreal, regional
#Products in both: mean, cv, extrapolation, sampling

#LIST OF SAMPLING PRODUCT:
#All the extrapolation
#Sampling density
#Distance to sampling