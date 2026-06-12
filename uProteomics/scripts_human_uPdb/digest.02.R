library(cleaver)
library(seqinr)
library(fs)

# this is relatively faster than digest.01.R, but only saves at the end

inp <- "db_unique"
out <- "db_unique_digested"

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
    cat('   ', length(seqs), '..')
    # open the output database
    fa.file <- file(oudb)
    # reset the database contents
    fa.lines <- c()
    #fa.list <- list()
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
            # save all peptides, we'll clean later
            if (len < 8 || len > 30) {
                #cat('KO\n')
                next
            } else {
                #cat('OK\n')
                npep <- npep + 1
                fa.head <- paste('>', attr(seq, 'name'), '_', i, gsub('>', ' : ', attr(seq, 'Annot')), sep='')
                fa.seq <- pep
                # add header an sequence to output
                fa.lines <- c(fa.lines, fa.head, fa.seq)
                #fa.entry <- paste('>', attr(seq, 'name'), '_', i, gsub('>', ' : ', attr(seq, 'Annot')), '\n', pep, sep='')
                #fa.list[[npep]] <- fa.entry
            }
        }
    }
    cat('', npep, '\n')
    # save peptide database
    writeLines(fa.lines, fa.file)
    #fwrite(fa.list, sep='\n')
    close(fa.file)
    #stop()
}
