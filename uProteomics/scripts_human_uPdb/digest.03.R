library(cleaver)
library(seqinr)
library(fs)
library(data.table)

# this is relatively faster than digest.02.R

#inp <- "db_unique"
#out <- "db_unique_digested"
inp <- "ucat.1"
out <- "ucat.1_digested"

minl <- 8       # minimum peptide length (val. included)
maxl <- 30      # maximum peptide length (val. included)

db <- list.files(inp, "*.faa$")

for (fasta in db) {
    indb <- paste(inp, fasta, sep='/')
    oudb <- paste(out, fasta, sep='/')
    if ( file.exists(oudb) ) {
        # already done, skip
        next
    }
    cat('Processing', indb, '\n')
    file_create(oudb)         # ensure out file exists so we can run in parallel
    # read the input database
    seqs <- read.fasta(indb, seqtype="AA", as.string=TRUE, set.attributes=TRUE)
    cat('    ', length(seqs), '..')
    # open the output database
    #fa.file <- file(oudb)
    # reset the database contents
    #fa.lines <- c()
    fa.list <- list()
    # for each sequence in the database
    npep=0
    for (seq in seqs) {
        #cat('Name:', attr(seq, 'name'), '\n')
        #cat('Annot:', attr(seq, 'Annot'), '\n')
        #cat('        ', seq[1], '\n')
        cleavage <- cleave(seq[1], enzym="trypsin")
        peptides <- cleavage[[1]]
        i=0
        for (pep in peptides) {
            len <- nchar(pep)
            i <- i + 1  #; if (i > 3) break
            #cat('    ', i, pep, '(', len, ') ')
            # save all peptides within range, we'll clean later
            if (len < minl || len > maxl) {
                #cat('KO\n')
                next
            } else {
                npep <- npep + 1
                #cat('OK\n')
                #fa.head <- paste('>', attr(seq, 'name'), '_', i, gsub('>', ' : ', attr(seq, 'Annot')), sep='')
                #fa.seq <- pep
                ##fa.lines <- c(fa.lines, fa.head, fa.seq)
                fa.entry <- paste('>', attr(seq, 'name'), '_', i, gsub('>', ' : ', attr(seq, 'Annot')), '\n', pep, sep='')
                fa.list[[npep]] <- fa.entry
            }
        }
    }
    cat(' ', npep, '\n', sep='')
    # save peptide database
    #writeLines(fa.lines, fa.file)
    #close(fa.file)
    cat('    Saving to ', oudb, '\n', sep='')
    fwrite(fa.list, file=oudb, sep='\n', quote=FALSE)
    #stop()
}
