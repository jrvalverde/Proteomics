library(cleaver)
library(seqinr)
library(fs)


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
    # read the database
    seqs <- read.fasta(indb, seqtype="AA", as.string=TRUE, set.attributes=TRUE)
    cat('   ', length(seqs), '..')
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
            i=i+1  #; if (i > 3) break
            #cat('    ', i, pep, '(', len, ') ')
            # save all peptides, we'll clean later
            if (len < 8 || len > 30) {
                #cat('KO\n')
                next
            } else {
                #cat('OK\n')
                npep <- npep + 1
                cat('>', attr(seq, 'name'), '_', i, gsub('>', ' : ', attr(seq, 'Annot')), '\n',
                    pep, '\n',
                    sep='',
                    file=oudb,
                    append=TRUE)
            }
        }
    }
    cat('', npep, '\n')
}
