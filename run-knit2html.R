.PATH = sprintf('.tmp/%s', strsplit(x=basename(OUT), split="\\.")[[1]][1])
dir.create(.PATH, showWarnings = F, recursive= T)
file.copy(from=IN, to=.PATH, overwrite=T)
rmarkdown::render(input=paste(.PATH,basename(IN),sep="/"), output_dir=dirname(OUT), 
                  output_file=basename(OUT), clean=T)