// vim: ts=4:expandtab:sw=4:cindent
//////////////////////////////////////////////////////////////////
// 
// Pipeline stage to run CODEX on exome data
//
//////////////////////////////////////////////////////////////////

CODEX_HG38_RELOAD_CODE = """

        library(BSgenome.Hsapiens.UCSC.hg38)

        getgc = function (chr, ref) {
          if (chr == "X" | chr == "x" | chr == "chrX" | chr == "chrx") {
            chrtemp <- 23
          } else if (chr == "Y" | chr == "y" | chr == "chrY" | chr == "chry") {
            chrtemp <- 24
          } else {
            chrtemp <- as.numeric(mapSeqlevels(as.character(chr), "NCBI")[1])
          }
          if (length(chrtemp) == 0) message("Chromosome cannot be found in NCBI Homo sapiens database!")
          chrm <- unmasked(Hsapiens[[chrtemp]])
          seqs <- Views(chrm, ref)
          af <- alphabetFrequency(seqs, baseOnly = TRUE, as.prob = TRUE)
          gc <- round((af[, "G"] + af[, "C"]) * 100, 2)
          gc
        }

""".stripIndent()

codex_call_cnvs_combined = {
    
    var batch_name : false,
        k_offset : 0,
        max_k : 9,
        codex_deletion_threshold: 1.7,
        codex_duplication_threshold: 2.3,
        codex_copynumber_mode : 'integer', // integer or fraction
        build : 19

    def outputFile = batch_name ? batch_name + '.cnvs.tsv' : input.bam + '.codex.cnvs.tsv'

    from(analysable_target) produce(outputFile) {
        
        def bamDir = file(input.bam).parentFile.absoluteFile.absolutePath
        
        def hg38_reload = ""
        
        if(build == 38) {
            hg38_reload = CODEX_HG38_RELOAD_CODE
        }
            
    
        R({"""
            source("$TOOLS/r-utils/cnv_utils.R")
            library(CODEX)

            $hg38_reload

            source("$TOOLS/codex/codex_segment_targeted.R")

            bam.files = c('${inputs.bam.join("','")}')

            bam.samples = as.matrix(data.frame(sample=sapply(bam.files, function(file.name) {
                # Scan the bam header and parse out the sample name from the first read group
                read.group.info = strsplit(scanBamHeader(file.name)[[1]]$text[["@RG"]],":")
                names(read.group.info) = sapply(read.group.info, function(field) field[[1]]) 
                return(read.group.info$SM[[2]])
            })))

            bedFile <- "$input.bed"
            targ.chr <- unique(as.matrix(read.table(bedFile, sep = "\\t")[,1]))
            gene.all=as.matrix(read.table(bedFile,head=F,sep='\\t')[,4])

            chr=targ.chr[[1]]

            bambedObj <- getbambed(bamdir = bam.files, bedFile = bedFile,
                                   sampname = bam.samples, 
                                   projectname = "ximmer", chr)
            
            bamdir <- bambedObj$bamdir; sampname <- bambedObj$sampname
            
            ref <- bambedObj$ref; 
            projectname <- bambedObj$projectname; 
            
            bambedObj <- getbambed(bamdir = bam.files, bedFile = bedFile,
                                   sampname = bam.samples, 
                                   projectname = "ximmer", chr)
            
            bamdir=bambedObj$bamdir; sampname=bambedObj$sampname; ref=bambedObj$ref; projectname=bambedObj$projectname;chr=bambedObj$chr
            # get raw depth of coverage
            coverageObj=getcoverage(bambedObj,mapqthres=20)
            Y=coverageObj$Y; readlength=coverageObj$readlength
            # get gc content
            gc=getgc(chr,ref)
            # get mappability
            mapp=getmapp(chr,ref)
            
            ref.all=bambedObj$ref
            Y.all=coverageObj$Y
            gc.all=gc
            mapp.all=mapp
            chr.all=rep(chr,length=length(mapp))
            
            for(chr in targ.chr[2:length(targ.chr)]) {
              print(chr)
              if(!is.element(chr,targ.chr)) next
              # get bam directories, read in bed file, get sample names
              bambedObj <- getbambed(bamdir = bam.files, bedFile = bedFile,
                                     sampname = bam.samples, 
                                     projectname = "ximmer", chr)
              
              bamdir=bambedObj$bamdir; sampname=bambedObj$sampname; ref=bambedObj$ref; projectname=bambedObj$projectname;chr=bambedObj$chr
              # get raw depth of coverage
              coverageObj=getcoverage(bambedObj,mapqthres=20)
              Y=coverageObj$Y; readlength=coverageObj$readlength
              # get gc content
              gc=getgc(chr,ref)

              # get mappability
              if("$build" == "38") {
                  mapp=rep(1,length(gc)) 
              } else {
                  mapp <- getmapp(chr, ref)
              }
              
              ref.all=c(ref.all,bambedObj$ref)
              Y.all=rbind(Y.all,coverageObj$Y)
              gc.all=c(gc.all,gc)
              mapp.all=c(mapp.all,mapp)
              chr.all=c(chr.all,rep(chr,length=length(mapp)))
            }

            Y=Y.all
            ref=ref.all
            gc=gc.all
            mapp=mapp.all
            gene = gene.all

            #------------------------------------

            qcObj <- qc(Y, sampname, chr, ref, mapp, gc, cov_thresh = c(20, 8000),
                        length_thresh = c(20, 4000), 
                        mapp_thresh = 0.9, 
                        gc_thresh = c(20, 80))
            

            Y_qc <- qcObj$Y_qc; 
            sampname_qc <- qcObj$sampname_qc; 
            gc_qc <- qcObj$gc_qc
            mapp_qc <- qcObj$mapp_qc; 
            ref_qc <- qcObj$ref_qc; 
            qcmat <- qcObj$qcmat

            gene_qc=gene[which(as.logical(qcmat[,4])==TRUE)]
            chr_qc=chr.all[which(as.logical(qcmat[,4])==TRUE)]

            
            # Rarely normalisation fails if k is too high
            # Here we loop while reducing max_k until it works
            # If this is happening consistently then better to reduce the configured initial max_k
            # via the parameter.
            max_k=$max_k; 
            normObj = F;  
            while(typeof(normObj) == 'logical') { 
                try(normObj <- normalize(Y_qc, gc_qc, K = 1:max_k)); max_k=max_k-1; 
            }
            
            Yhat <- normObj$Yhat; AIC <- normObj$AIC; BIC <- normObj$BIC

            RSS <- normObj$RSS; 
            K <- normObj$K 
            
            optK = min(9,max(1,K[which.max(BIC)] + $k_offset))


            finalcall=matrix(ncol=14)
            
            #--------------------------------------
            
            colnames(finalcall)=c('sample_name','chr','gene','cnv',
                                  'st_bp','ed_bp','length_kb',
                                  'st_exon','ed_exon','raw_cov',
                                  'norm_cov','copy_no','lratio',
                                  'mBIC')
            for(genei in unique(gene_qc)){
              cat('Segmenting gene',genei,'\\n')
              geneindex=which(gene_qc==genei)
              yi=Y_qc[geneindex,]
              yhati=Yhat[[optK]][geneindex,]
              refi=ref_qc[geneindex]
              chri=chr_qc[geneindex][1]
              finalcalli=segment_targeted(yi, yhati, sampname_qc, refi, genei, chri, lmax=length(geneindex), mode='$codex_copynumber_mode') 
              finalcall=rbind(finalcall,finalcalli)
            }
            
            finalcall=finalcall[-1,]
            cn=(as.numeric(as.matrix(finalcall[,'copy_no'])))
            if("$codex_copynumber_mode" == "fraction") {
                cn.filter=(cn<=$codex_deletion_threshold)|(cn>=$codex_duplication_threshold)
                finalcall=finalcall[cn.filter,]
            } else {
                finalcall=finalcall[cn!=2,]
            }

            write.table(finalcall[,-3], # to be concordant with chr-split results (below), remove gene column
                        file="$output.tsv",
                        row.names=F,
                        col.names=T,
                        quote=F,
                        sep="\\t"
                        )
        """}, "codex")
    }
    
    branch.caller_result = output.tsv
}

codex_call_cnvs = {

    var batch_name : false,
        k_offset : 0,
        max_k : 9,
        build: 19

    def chr = branch.name
    
    
    if(!(chr in analysable_chromosomes))
        succeed "Chromosome $chr is not analysable by CODEX"

    if(filter_to_sex == "FEMALE" && chr == "Y" || chr == "chrY") {
        succeed "Ignoring Y chromosome due to sex filtering to female samples"
    }
    else {
        println "Analysing $chr with CODEX"
    }

    def hg38_reload = ""
    
    if(build == 38) {
        hg38_reload = CODEX_HG38_RELOAD_CODE
    }

    
    def outputFile = batch_name ? batch_name + '.codex.' + chr + '.cnvs.tsv' : input.bam + '.codex.' + chr + '.cnvs.tsv'

    produce(outputFile) {
        
        def bamDir = file(input.bam).parentFile.absoluteFile.absolutePath

        R({"""

            source("$TOOLS/r-utils/cnv_utils.R")

            library(CODEX)

            $hg38_reload

            bam.files = c('${inputs.bam.join("','")}')

            bam.samples = as.matrix(data.frame(sample=sapply(bam.files, function(file.name) {
                # Scan the bam header and parse out the sample name from the first read group
                read.group.info = strsplit(scanBamHeader(file.name)[[1]]$text[["@RG"]],":")
                names(read.group.info) = sapply(read.group.info, function(field) field[[1]]) 
                return(read.group.info$SM[[2]])
            })))

            bedFile <- "$input.bed"
            chr <- "$chr"

            bambedObj <- getbambed(bamdir = bam.files, bedFile = bedFile,
                                   sampname = bam.samples, 
                                   projectname = "ximmer", chr)
            
            bamdir <- bambedObj$bamdir; sampname <- bambedObj$sampname
            
            
            ref <- bambedObj$ref; 
            projectname <- bambedObj$projectname; 
            chr <- bambedObj\$chr
            
            
            coverageObj <- getcoverage(bambedObj, mapqthres = 20)
            
            Y <- coverageObj$Y; readlength <- coverageObj$readlength
            
            gc <- getgc(chr, ref)
            
            if("$build" == "38") {
                mapp=rep(1,length(gc)) 
            } else {
                mapp <- getmapp(chr, ref)
            }
            
            qcObj <- qc(Y, sampname, chr, ref, mapp, gc, cov_thresh = c(20, 4000),
                        length_thresh = c(20, 2000), 
                        mapp_thresh = 0.9, 
                        gc_thresh = c(20, 80))
            
            Y_qc <- qcObj$Y_qc; 
            sampname_qc <- qcObj$sampname_qc; 
            gc_qc <- qcObj$gc_qc
            
            mapp_qc <- qcObj$mapp_qc; ref_qc <- qcObj$ref_qc; qcmat <- qcObj$qcmat
            
            # Rarely normalisation fails if k is too high
            # Here we loop while reducing max_k until it works
            # If this is happening consistently then better to reduce the configured initial max_k
            # via the parameter.
            max_k=$max_k; 
            normObj = F;  
            while(typeof(normObj) == 'logical') { 
                try(normObj <- normalize(Y_qc, gc_qc, K = 1:max_k)); max_k=max_k-1; 
            }
            
            Yhat <- normObj$Yhat; AIC <- normObj$AIC; BIC <- normObj$BIC

            RSS <- normObj$RSS; 
            K <- normObj$K 
            
            optK = min(9,max(1,K[which.max(BIC)] + $k_offset))
            
            finalcall <- segment(Y_qc, Yhat, optK = optK, K = K, sampname_qc,
                                 ref_qc, chr, lmax = 200, mode = "integer")

            write.table(finalcall,
                        file="$output.tsv",
                        row.names=F,
                        col.names=T,
                        quote=F,
                        sep="\\t"
                        )

        """}, "codex")
    }
    branch.caller_result = output.tsv
}

merge_codex = {
        exec """
            head -1 $input.cnvs.tsv > $output.tsv; 

            for i in $inputs.cnvs.tsv; do echo "Merge $i"; grep -v '^sample_name' $i  >> $output.tsv || true ; done

            echo "Done."
        """
        
        branch.caller_result = output.tsv
    
}

if(codex_split_chrs) {
    codex_pipeline = segment {
        chromosomes * [ codex_call_cnvs ] >> merge_codex
    }
}
else {
    codex_pipeline = segment {
         codex_call_cnvs_combined 
    }
}
