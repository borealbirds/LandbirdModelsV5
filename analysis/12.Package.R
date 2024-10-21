#TO DO: PACKAGE VALIDATION: MEAN & SD OF VALIDATIONS######

#17. Package and save----
out <- data.table::rbindlist(out.list, fill=TRUE)

write.csv(out, file.path(root, "output", "validation", "ModelValidation_BCR.csv"), row.names = FALSE)


#TO DO: SPECIES GUILDS STACKS#######

#TO DO: FINAL PREDICTION RASTER STACKS
#REPROJECT TO LAEA HERE

#LIST OF PRODUCTS PER SPECIES:
#2 extents: boreal, regional
#Products in both: mean, sd, extrapolation, sampling

#LIST OF SAMPLING PRODUCT:
#All the extrapolation
#Sampling density
#Distance to sampling