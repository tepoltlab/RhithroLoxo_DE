## run with this snakemake command: 

## snakemake --jobs 20 --use-conda --cluster-config envs/cluster.yaml --cluster "sbatch --parsable --partition={cluster.queue} --job-name=RhithroLoxo_DE.{rule}.{wildcards} --mem={cluster.mem}gb --time={cluster.time} --ntasks={cluster.threads} --nodes={cluster.nodes}"

## DEFINE VARIABLES 

#configfile: 'envs/trinity_config.yaml'

WILDCARDS = glob_wildcards('trimmed_reads/{sample}.fq') #define wildcards object that will contain sample names
FASTQC_ZIP = expand('outputs/fastqc/{wc1}_fastqc.zip', wc1=WILDCARDS.sample)
FASTQC_HTML = expand('outputs/fastqc/{wc1}_fastqc.html', wc1=WILDCARDS.sample)
BLASTDB_EXT = ["nhr", "nin", "nog", "nsd", "nsi", "nsq"]
LOXO_DB = expand('txms/loxo_db/loxo_db.{wc1}', wc1=BLASTDB_EXT)
LOXO_BLAST_RESULTS = expand('outputs/blast2loxo/{wc1}_loxohits.tsv', wc1=WILDCARDS.sample)
CONTAM_LIST = expand('outputs/blast2loxo/{wc1}_contam_IDs.txt', wc1=WILDCARDS.sample)
CONTAM_RATES = 'outputs/blast2loxo/contam_rates.txt'
CLEANED_READS = expand('reads_decontam/{wc1}_clean.fq', wc1=WILDCARDS.sample)
RHITHRO_CLEAN_CAT = 'cat_seq/cat_rhithro_clean.fq'
RHITHRO_TXMS = expand('txms/rhithro/trinity_kmer{wc1}cov{wc2}.Trinity.fasta', wc1=["21", "25", "29"], wc2=["1", "2"])
RHITHRO_TXM_GENETRANSMAPS = expand('txms/rhithro/trinity_kmer{wc1}cov{wc2}.Trinity.fasta.gene_trans_map', wc1=["21", "25", "29"], wc2=["1", "2"])
QUANT_CAT = expand('outputs/trinity_stats/rhithro_kmer{wc1}cov{wc2}/quant.sf{wc3}', wc1=["21", "25", "29"], wc2=["1", "2"], wc3=["",".genes"])
PRELIM_TXM_IDXS = expand('txms/rhithro/trinity_kmer{wc1}cov{wc2}.Trinity.fasta.salmon_quasi.idx', wc1=["21", "25", "29"], wc2=["1", "2"])
EXN50 = expand('outputs/trinity_stats/rhithro_kmer{wc1}cov{wc2}/ExN50.stats', wc1=["21", "25", "29"], wc2=["1", "2"])
N50 = expand('outputs/trinity_stats/rhithro_kmer{wc1}cov{wc2}/N50.txt', wc1=["21", "25", "29"], wc2=["1", "2"])
TXM_LONG = 'txms/rhithro/rhithro_txm_long.fasta'
PADDED_RANGE = [str(item).zfill(2) for item in range(0,50)]
DIRTY_CHUNKS = expand('txms/rhithro/dirty_chunks/rhithro_txm_long.{wc1}.fasta', wc1=PADDED_RANGE)
BLAST_CONTAM = expand('outputs/blast/blast_contam.{wc1}.tsv', wc1=PADDED_RANGE)
BLAST_CONTAM_MERGED = 'outputs/blast/blast_contam.merged.tsv'
TXM_LONG_CLEAN = 'txms/rhithro/rhithro_txm_long_clean.fasta'
BLAST_CONTAM_LIST = 'outputs/blast/contam_list.txt'
TRANSDECODER_OUT = expand('outputs/transdecoder/rhithro_txm_long_clean.fasta.transdecoder.{ext}', ext=["pep","cds","gff3","bed"])
RHITHRO_TXM_IDX = 'txms/rhithro/rhithro_txm_long_clean.fasta.salmon_quasi.idx'
DE_SAMPLES = ['AP_C_1','AP_C_2','AP_C_3','AP_C_4','AP_C_5','AP_C_6','AP_P_1','AP_P_2','AP_P_3','AP_P_6','FP_C_10','FP_C_11','FP_C_12','FP_C_13','FP_C_5','FP_C_9','FP_P_10','FP_P_4','FP_P_5','FP_P_7','FP_P_8','FP_P_9','LA_C_1','LA_C_2','LA_C_3','LA_C_4','LA_C_6','LA_C_8','LA_P_1','LA_P_2','MA_C_1','MA_C_2','MA_C_4','MD_C_10','MD_C_11','MD_C_12','MD_C_1','MD_C_4','MD_C_7','MD_P_1','ML_C_10','ML_C_2','ML_C_3','ML_C_5','ML_C_7','ML_C_9','ML_P_1','ML_P_2','NH_C_11','NH_C_12','NH_C_13','NH_C_5','NH_C_8','NH_C_9','NH_P_1','NH_P_2','NH_P_3','NH_P_4','NH_P_5','NH_P_6','NJ_C_10','NJ_C_11','NJ_C_12','NJ_C_13','NJ_C_14','NJ_C_6','NJ_P_1','NJ_P_3','NJ_P_4','NJ_P_5','NJ_P_6','NJ_P_7','SC_C_12','SC_C_14','SC_C_2','SC_C_6','SC_C_7','SC_C_9','SC_P_1','SC_P_2','SC_P_3']
QUANT_OUT = expand('outputs/quant/{wc1}/quant.sf', wc1=DE_SAMPLES)
QUANT_MAT = expand('outputs/quant/salmon.isoform.{wc1}', wc1 = ['counts.matrix','TPM.not_cross_norm'])
DOWNLOAD_DB_OUT = expand('db/{wc1}', wc1=['uniref90.fasta.gz','uniprot_sprot.pep','uniprot_trembl.fasta.gz','nr.gz','Pfam-A.hmm.h3f','uniprot_sprot.fasta.gz','refseq_complete.faa.gz'])
DIAMOND_DB = expand('db/{wc1}.dmnd', wc1 = ['nr','uniref90','trembl','sprot'])
DIAMONDP_ANNOT_OUT = expand('outputs/trinotate/diamondp.{wc1}.tsv', wc1=['nr','uniref90','trembl','sprot'])
DIAMONDX_ANNOT_OUT = expand('outputs/trinotate/diamondx.{wc1}.tsv', wc1=['nr','uniref90','trembl','sprot'])
DIAMONDX_XML = 'outputs/trinotate/diamondx.nr.xml'
HMMSCAN_OUT = 'outputs/trinotate/hmmscan.out'
GENE_TRANS_MAP = 'txms/rhithro/rhithro_txm_long_clean.gene_trans_map'
TRINOTATE_INIT = 'outputs/track/trinotate_init.done'
TRINOTATE_LOAD = 'outputs/track/trinotate_load.done'
TRINOTATE_ANNOT = 'outputs/trinotate/trinotate_annotations.tsv'
TRINOTATE_STATS = 'outputs/trinotate/trinotate_stats.txt'


#DIAMOND_ANNOT_OUT_SPLIT = expand('outputs/trinotate/{wc2}/diamond{wc3}.{wc1}.{wc2}.tsv',  wc1=PADDED_RANGE2, wc2=['nr','uniref90','trembl','sprot'], wc3=['x','p'])
#CLEAN_CHUNKS_NT = expand('txms/rhithro/clean_chunks/rhithro_txm_long_clean.{wc1}.fasta', wc1=PADDED_RANGE2)
#CLEAN_CHUNKS_AA = expand('txms/rhithro/clean_chunks/rhithro_txm_long_clean.{wc1}.fasta.transdecoder.pep', wc1=PADDED_RANGE2)


## GET SNAKEMAKEY

rule all:
    input: FASTQC_ZIP, FASTQC_HTML, LOXO_DB, LOXO_BLAST_RESULTS, CONTAM_LIST, CONTAM_RATES, CLEANED_READS, RHITHRO_CLEAN_CAT, RHITHRO_TXMS, RHITHRO_TXM_GENETRANSMAPS, PRELIM_TXM_IDXS, QUANT_CAT, EXN50, N50, TXM_LONG, DIRTY_CHUNKS, BLAST_CONTAM, BLAST_CONTAM_MERGED, BLAST_CONTAM_LIST, TXM_LONG_CLEAN, TRANSDECODER_OUT, RHITHRO_TXM_IDX, QUANT_OUT, QUANT_MAT, DOWNLOAD_DB_OUT, DIAMOND_DB, DIAMONDP_ANNOT_OUT, DIAMONDX_ANNOT_OUT, DIAMONDX_XML, HMMSCAN_OUT, GENE_TRANS_MAP, TRINOTATE_INIT, TRINOTATE_LOAD, TRINOTATE_ANNOT, TRINOTATE_STATS,
    
#DIAMOND_ANNOT_OUT_SPLIT
#CLEAN_CHUNKS_NT
#CLEAN_CHUNKS_AA

#get fastqc results on all samples

rule fastqc:
    input: 'trimmed_reads/{file}.fq'
    output:
        fastqc_zip = 'outputs/fastqc/{file}_fastqc.zip',
        fastqc_html = 'outputs/fastqc/{file}_fastqc.html'
    conda:
        "envs/fastqc.yaml" 
    log: 'logs/fastqc/{file}_fastqc.log'
    shell:
        """
        mkdir -p outputs/fastqc 
        mkdir -p logs/fastqc 
        fastqc -o outputs/fastqc {input} 1> {log} 2>&1 #run fastqc and output to folder
        """

#looking at HTML reports shows that all samples are good except for MA_C_3. samples.txt has this row removed

#blast trimmed reads to Carolyn's Loxo txm

#make blast database from the Loxo txm

rule makeloxodb:
    input: 
        loxo_txm = 'txms/Loxo_assembly_clean.fasta'
    output:
        loxo_db = LOXO_DB
    conda: 
        "envs/magicblast.yaml"
    shell:
        """
        mkdir -p txms/loxo_db
        makeblastdb -in {input.loxo_txm} \
            -out txms/loxo_db/loxo_db \
            -dbtype nucl \
            -parse_seqids \
            -title "loxo_db"
        """

#blast against loxo db. using magic blast because it can handle fastq

rule blast2loxo:
    input:
        db = LOXO_DB,
        reads = 'trimmed_reads/{file}.fq'
    output:
        results = 'outputs/blast2loxo/{file}_loxohits.tsv'
    conda: 
        "envs/magicblast.yaml"
    shell:
        """
        mkdir -p outputs/blast2loxo
        magicblast -query {input.reads} \
            -db txms/loxo_db/loxo_db \
            -infmt fastq \
            -reftype transcriptome \
            -outfmt tabular \
            -no_unaligned \
            -num_threads 35 \
            -out {output.results}
        """

#get mapping IDs
rule list_contam:   
    input:
        blast_results = 'outputs/blast2loxo/{file}_loxohits.tsv'
    output:
        seq_ids = 'outputs/blast2loxo/{file}_contam_IDs.txt'
    shell:
        """
        #print seq headers that align over entire length of read w/ 98% ID
        #field 8 is query alignment end, field 7 is query alignment start, field 16 is read length, field 3 is PID
        #so $8-$7+1 is the alignment length. if equals read length ($16) and PID>=98, print header and sort/uniq
        tail -n+4 {input.blast_results} | awk '$8-$7+1>=$16 && $3>=98  {{print $1}}' | sort | uniq > {output.seq_ids}
        """

rule contam_rates:
    input: CONTAM_LIST
    output: CONTAM_RATES
    shell:
        """
        for file in `ls outputs/blast2loxo/*contam*`
        do
            wc -l $file >> outputs/blast2loxo/contam_count.txt
        done
        
        sort -k2 -o outputs/blast2loxo/contam_count.txt outputs/blast2loxo/contam_count.txt
        
        for file in `ls trimmed_reads/*.fq`
        do
            awk '{{OFS="\t"}}{{s++}}END{{print FILENAME,s/4}}' $file >> outputs/blast2loxo/read_count.txt
        done
        
        sort -k1 -o outputs/blast2loxo/read_count.txt outputs/blast2loxo/read_count.txt
        
        paste outputs/blast2loxo/read_count.txt outputs/blast2loxo/contam_count.txt > outputs/blast2loxo/combined.txt
        
        awk '{{OFS="\t"}}{{print $1,$2,$3,$3/$2}}' outputs/blast2loxo/combined.txt > {output}
        
        rm outputs/blast2loxo/combined.txt
        rm outputs/blast2loxo/read_count.txt
        rm outputs/blast2loxo/contam_count.txt
        """

rule remove_contam:
    input: 
        contam_list = 'outputs/blast2loxo/{file}_contam_IDs.txt',
        trimmed_seq = 'trimmed_reads/{file}.fq'
    output:
        cleaned = 'reads_decontam/{file}_clean.fq'
    conda:
        "envs/bbtools.yaml"
    shell:
        """
        mkdir -p reads_decontam
        filterbyname.sh in={input.trimmed_seq} out={output.cleaned} names={input.contam_list} 
        """

rule cat_clean:
    input:
        reads = CLEANED_READS
    output:
        rhithro_clean = RHITHRO_CLEAN_CAT
    shell:
        """
        mkdir -p cat_seq
        cat $(awk '$1== "uninfected" {{print $3}}' metadata/samples.txt) >> {output.rhithro_clean}
        """
       
#make rhithro txms using different parameter values 
       
rule rhithro_txm:
    input: 
        seq = RHITHRO_CLEAN_CAT
    conda: 
        "envs/trinity.yaml"
    log:
        'logs/trinity/rhithro_kmer{ksize}cov{kcov}.log'
    output: 
        txm = 'txms/rhithro/trinity_kmer{ksize}cov{kcov}.Trinity.fasta',
        gene_trans_map = 'txms/rhithro/trinity_kmer{ksize}cov{kcov}.Trinity.fasta.gene_trans_map'
    shell:
        """
        mkdir -p logs/trinity
        
        Trinity --seqType fq \
            --SS_lib_type R \
            --max_memory 120G \
            --CPU 35 \
            --single {input.seq} \
            --KMER_SIZE {wildcards.ksize} \
            --min_kmer_cov {wildcards.kcov} \
            --full_cleanup \
            --output txms/rhithro/trinity_kmer{wildcards.ksize}cov{wildcards.kcov} 1> {log} 2>&1
        
        rm -rf txms/rhithro/trinity_kmer{wildcards.ksize}cov{wildcards.kcov}
        """

#prep the reference directly in salmon b/c can't pass kmer options to salmon index from trinity

rule prep_ref:
    input:
        txm = 'txms/rhithro/trinity_kmer{wc1}cov{wc2}.Trinity.fasta'
    conda:
        "envs/trinity.yaml"
    output:
        idx = directory('txms/rhithro/trinity_kmer{wc1}cov{wc2}.Trinity.fasta.salmon_quasi.idx')
    log:
        'logs/trinity/prep_ref_kmer{wc1}cov{wc2}.log'
    shell:
        """
        salmon index -t {input.txm} \
            --keepDuplicates \
            -i {output.idx} \
            --type quasi \
            -k 25 \
            -p 16 1> {log} 2>&1
        """

#run salmon w/in trinity to get count data for ExN50 comparisons

rule quant_cat:
    input:
        txm = 'txms/rhithro/trinity_kmer{wc1}cov{wc2}.Trinity.fasta',
        cat = RHITHRO_CLEAN_CAT,
        gtm = 'txms/rhithro/trinity_kmer{wc1}cov{wc2}.Trinity.fasta.gene_trans_map',
        idx = 'txms/rhithro/trinity_kmer{wc1}cov{wc2}.Trinity.fasta.salmon_quasi.idx'
    conda:
        "envs/trinity.yaml"
    output:
        quant1 = 'outputs/trinity_stats/rhithro_kmer{wc1}cov{wc2}/quant.sf',
        quant2 = 'outputs/trinity_stats/rhithro_kmer{wc1}cov{wc2}/quant.sf.genes'
    log:
        'logs/trinity/quant_cat_kmer{wc1}cov{wc2}.log'
    shell:
        """
        mkdir -p outputs/trinity_stats
        
        align_and_estimate_abundance.pl --transcripts {input.txm} \
            --seqType fq \
            --single {input.cat} \
            --est_method salmon \
            --gene_trans_map {input.gtm} \
            --output_dir outputs/trinity_stats/rhithro_kmer{wildcards.wc1}cov{wildcards.wc2} \
            --SS_lib_type R \
            --thread_count 16 \
            --fragment_length 150 \
            --fragment_std 20 \
            --salmon_add_opts "--validateMappings " 1> {log} 2>&1 
        """

rule txm_stats:
    input:
        quant = 'outputs/trinity_stats/rhithro_kmer{wc1}cov{wc2}/quant.sf'
    conda:
        "envs/trinity.yaml"
    params:
        dir = directory('outputs/trinity_stats/rhithro_kmer{wc1}cov{wc2}')
    output:
        ExN50 = 'outputs/trinity_stats/rhithro_kmer{wc1}cov{wc2}/ExN50.stats',
        N50 = 'outputs/trinity_stats/rhithro_kmer{wc1}cov{wc2}/N50.txt' 
    log:
        'logs/trinity/txm_stats_kmer{wc1}cov{wc2}.log'
    shell:
        """
        cd {params.dir}
        abundance_estimates_to_matrix.pl --est_method salmon \
            --gene_trans_map none quant.sf
        contig_ExN50_statistic.pl salmon.isoform.TPM.not_cross_norm \
            ../../../txms/rhithro/trinity_kmer{wildcards.wc1}cov{wildcards.wc2}.Trinity.fasta | tee ExN50.stats
        TrinityStats.pl ../../../txms/rhithro/trinity_kmer{wildcards.wc1}cov{wildcards.wc2}.Trinity.fasta > N50.txt 
        """

#okay now we have the ExN50 and N50 statistics. See jupyter_notebooks/txm_compare.ipynb for plots and comparisons.
#the default of kmer25cov1 looks the best. highest Ex90N50.

#retain only longest isoform per 'gene'

rule txm_long:
    input:
        txm = 'txms/rhithro/trinity_kmer25cov1.Trinity.fasta',
        script1 = {'scripts/list-trinity-dups-to-remove.py'},
        script2 = {'scripts/fasta_subsetter.py'}
    output:
        txm_long = TXM_LONG
    shell:
        """
        python {input.script1} {input.txm}
        python {input.script2} {input.txm} txms/rhithro/trinity_kmer25cov1.Trinity.fasta_dups.txt REMOVE
        mv txms/rhithro/trinity_kmer25cov1.Trinity.fasta_dups_REMOVE.fasta {output.txm_long}
        """
        
#split the unfiltered txm for blasting contam

rule split_dirty:
    input:
        txm = 'txms/rhithro/rhithro_txm_long.fasta'
    output:
        chunks = DIRTY_CHUNKS
    params:
        dir = directory('txms/rhithro/dirty_chunks')
    conda:
        "envs/pyfasta.yaml"
    shell:
        """
        mkdir -p {params.dir}
        pyfasta split -n 50 {input.txm}
        mv txms/rhithro/rhithro_txm_long.*.fasta {params.dir}
        mv txms/rhithro/rhithro_txm_long.fasta.flat {params.dir}
        mv txms/rhithro/rhithro_txm_long.fasta.gdx {params.dir}
        """

#blast chunks against nt. db release Dec 31 2018

rule blast_contam:
    input:
        txm = 'txms/rhithro/dirty_chunks/rhithro_txm_long.{wc1}.fasta'
    output:
        results = 'outputs/blast/blast_contam.{wc1}.tsv'
#    conda:
#        "envs/blast.yaml" #gave up on trying to use posedion blastdb with with conda install of blast. using poseidon's module. 2.7.1
    shell:
        """
        module load bio blast/2.7.1
        
        mkdir -p outputs/blast
        blastn -query {input.txm} \
            -db nt \
            -evalue 1e-10 \
            -max_target_seqs 1 \
            -max_hsps 1 \
            -num_threads 35 \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids" \
            -out {output.results}       
        """
        
#merge blast results

rule merge_blast:
    input: 
        BLAST_CONTAM
    output:
        BLAST_CONTAM_MERGED
    params:
        dir = directory('outputs/blast')
    shell:
        """
        cd {params.dir}
        cat blast_contam.*.tsv >> blast_contam.merged.tsv
        """

#id contam: bacteria,archaea,viruses,platyhelminthes,nematoda,fungi,alveolata,viridiplantae,rhodophyta,amoebozoa,rhizaria,stramenopiles,rhizocephala,entoniscidae

rule id_contam:
    input: 
        BLAST_CONTAM_MERGED,
        'scripts/map_contam_ids.py'
    output:
        BLAST_CONTAM_LIST
    conda:
        "envs/ete3.yaml"
    shell:
        """
        python {input[1]} {input[0]} {output}
        """

#remove contam contigs

rule remove_blast_contam:
    input:
        txm = TXM_LONG,
        contam_list = BLAST_CONTAM_LIST,
        script = {'scripts/fasta_subsetter.py'}
    output:
        txm_long = TXM_LONG_CLEAN
    shell:
        """
        python {input.script} {input.txm} {input.contam_list} REMOVE
        mv outputs/blast/contam_lis_REMOVE.fasta {output.txm_long} #'contam_lis' instead of 'contam_list' because .strip() bug in fasta_subsetter. ignore for now, fix later
        """

rule transdecoder:
    input:
        txm = TXM_LONG_CLEAN
    output:
        output = TRANSDECODER_OUT
    conda:
        "envs/transdecoder.yaml"
    params:
        dir = directory('outputs/transdecoder')
    shell:
        """
        TransDecoder.LongOrfs -t {input.txm} -O {params.dir}
        TransDecoder.Predict -t {input.txm} -O {params.dir}
        mv rhithro_txm_long_clean.fasta.transdecoder.* outputs/transdecoder
        """
#quantifying all except for MA_C_3 (b/c bad seq). leaving in NJ_P_4 even though it clustered with NH in pop gen, probably mislabeled. even if it's misassigned to population, still in the absent range and still infected. shouldn't matter

rule prep_rhithro_ref:
    input:
        txm = TXM_LONG_CLEAN
    conda:
        "envs/trinity.yaml"
    output:
        idx = directory('txms/rhithro/rhithro_txm_long_clean.fasta.salmon_quasi.idx')
    log:
        'logs/trinity/prep_rhithro_ref.log'
    shell:
        """
        salmon index -t {input.txm} \
            --keepDuplicates \
            -i {output.idx} \
            --type quasi \
            -k 25 \
            -p 16 1> {log} 2>&1
        """

rule quant:
    input:
        txm = TXM_LONG_CLEAN,
        metadata = 'metadata/trinity_samples.txt',
        idx = RHITHRO_TXM_IDX,
        reads = expand('reads_decontam/{sample}_clean.fq', sample=DE_SAMPLES)
    output:
        output = QUANT_OUT
    conda:
        "envs/trinity.yaml"
    params:
        dir = 'outputs/quant'
    log:
        'logs/trinity/quant.log'
    shell:
        """
        mkdir -p {params.dir}
        cd {params.dir}
        align_and_estimate_abundance.pl --transcripts ../../{input.txm} \
            --seqType fq \
            --samples_file ../../{input.metadata} \
            --est_method salmon \
            --SS_lib_type R \
            --thread_count 16 \
            --fragment_length 150 \
            --fragment_std 20 \
            --salmon_add_opts "--validateMappings " 1> ../../{log} 2>&1 
        """

rule quant_mat:
    input:
        quant_files = QUANT_OUT,
        quant_paths = 'metadata/quant_path.txt' #change this file if you want to include/exclude certain samples, or just keep all and remove desired columns downstream in R, since they aren't cross-sample normed
    output:
        QUANT_MAT
    conda:
        "envs/trinity.yaml"
    params:
        dir = 'outputs/quant'
    log:
        'logs/trinity/quant_mat.log'
    shell:
        """
        mkdir -p {params.dir}
        cd {params.dir}
        abundance_estimates_to_matrix.pl --est_method salmon \
            --gene_trans_map 'none'\
            --cross_sample_norm none \
            --name_sample_by_basedir \
            --quant_files ../../{input.quant_paths} 1> ../../{log} 2>&1
        """

#get pfam and uniprot db and build Trinotate sqlite db. current releases as of 12/29/19
#if i want to run EnTap, i have to redownload the sprot in fasta format (trinotate changed headers)

rule download_db:
    input:
    output:
        DOWNLOAD_DB_OUT
    params:
        dir = directory('db/')
    conda:
        "envs/trinotate.yaml"
    shell:
        """
        cd {params.dir}
        Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
        gunzip Pfam-A.hmm.gz
        hmmpress Pfam-A.hmm
        wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
        wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
        wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
        wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
        wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/complete.nonredundant_protein.*.protein.faa.gz
        cat complete.nonredundant_protein.*.protein.faa.gz > refseq_complete.faa.gz
        rm complete.nonredundant_protein.*.protein.faa.gz
        """

#prepare diamond databases

rule diamond_db: 
    input:
        nr = 'db/nr.gz',
        sprot = 'db/uniprot_sprot.pep',
        uniref = 'db/uniref90.fasta.gz',
        trembl = 'db/uniprot_trembl.fasta.gz'
    output:
        DIAMOND_DB
    conda:
        "envs/diamond.yaml"
    params:
        dir = directory('db')
    shell:
        """
        cd {params.dir}
        diamond makedb --in nr.gz --db nr
        diamond makedb --in uniref90.fasta.gz --db uniref90
        diamond makedb --in uniprot_trembl.fasta.gz --db trembl
        diamond makedb --in uniprot_sprot.pep --db sprot
        """

rule diamondp_annot:
    input:
        aa = 'outputs/transdecoder/rhithro_txm_long_clean.fasta.transdecoder.pep',
        db = 'db/{db}.dmnd'
    output:
        hits = 'outputs/trinotate/diamondp.{db}.tsv'
    conda:
        "envs/diamond.yaml"
    log:
        aa = 'logs/trinotate/diamondp_{db}.log',
    shell:
        """
        diamond blastp --query {input.aa} \
            --db {input.db} \
            --out {output.hits} \
            --threads 36 \
            --outfmt 6 \
            --evalue 0.001 \
            --sensitive \
            --block-size 8 \
            --index-chunks 1 \
            --max-target-seqs 1 1> {log.aa} 2>&1
        """

rule diamondx_annot:
    input:
        nt = TXM_LONG_CLEAN,
        db = 'db/{db}.dmnd'
    output:
        hits = 'outputs/trinotate/diamondx.{db}.tsv',
    conda:
        "envs/diamond.yaml"
    log:
        nt = 'logs/trinotate/diamondx_{db}.log'
    shell:
        """
        diamond blastx --query {input.nt} \
            --db {input.db} \
            --out {output.hits} \
            --threads 36 \
            --outfmt 6 \
            --evalue 0.001 \
            --sensitive \
            --block-size 8 \
            --index-chunks 1 \
            --max-target-seqs 1 1> {log.nt} 2>&1
        """

rule diamondx_xml: #for blast2go
    input:
        nt = TXM_LONG_CLEAN,
        db = 'db/nr.dmnd'
    output:
        hits = 'outputs/trinotate/diamondx.nr.xml',
    conda:
        "envs/diamond.yaml"
    log:
        nt = 'logs/trinotate/diamondx_nr_xml.log'
    shell:
        """
        diamond blastx --query {input.nt} \
            --db {input.db} \
            --out {output.hits} \
            --threads 36 \
            --outfmt 5 \
            --evalue 0.001 \
            --sensitive \
            --block-size 8 \
            --index-chunks 1 \
            --max-target-seqs 1 1> {log.nt} 2>&1
        """

rule hmmscan:
    input:
        aa = 'outputs/transdecoder/rhithro_txm_long_clean.fasta.transdecoder.pep',
        pfam = DOWNLOAD_DB_OUT #for tracking. #Trinotate.sqlite removed from list b/c of later rules. assume download works
    output:
        hits = HMMSCAN_OUT
    conda:
        "envs/trinotate.yaml"
    log:
        'logs/trinotate/hmmscan.log'
    shell:
        """
        hmmscan --cpu 36 \
            --domtblout {output.hits} \
            db/Pfam-A.hmm {input.aa} 1> {log} 2>&1
        """

rule make_gene_trans_map:
    input:
        TXM_LONG_CLEAN
    output:
        GENE_TRANS_MAP
    conda:
        "envs/trinity.yaml"
    shell:
        """
        get_Trinity_gene_to_trans_map.pl {input} > {output}
        """

rule trinotate_init:
    input:
        txm_nt = TXM_LONG_CLEAN,
        txm_aa = 'outputs/transdecoder/rhithro_txm_long_clean.fasta.transdecoder.pep',
        gtm = GENE_TRANS_MAP,
        sqlite = ancient('db/Trinotate.sqlite') #ignore timestamp. output is same as input so do wan't to run again
    output:
        touch(TRINOTATE_INIT)
    conda:
        "envs/trinotate.yaml"
    log:
        "logs/trinotate/trinotate_init.txt"
    shell:
        """
        Trinotate {input.sqlite} init --gene_trans_map {input.gtm} \
            --transcript_fasta {input.txm_nt} \
            --transdecoder_pep {input.txm_aa} 1> {log} 2>&1
        """

rule trinotate_load:
    input:
        init_done = TRINOTATE_INIT,
        sqlite = ancient('db/Trinotate.sqlite'),
        blastp_sprot = 'outputs/trinotate/diamondp.sprot.tsv',
        blastx_sprot = 'outputs/trinotate/diamondx.sprot.tsv',
        blastp_trembl = 'outputs/trinotate/diamondp.trembl.tsv',
        blastx_trembl = 'outputs/trinotate/diamondx.trembl.tsv',
        blastp_uniref = 'outputs/trinotate/diamondp.uniref90.tsv',
        blastx_uniref = 'outputs/trinotate/diamondx.uniref90.tsv',
        blastp_nr = 'outputs/trinotate/diamondp.nr.tsv',
        blastx_nr = 'outputs/trinotate/diamondx.nr.tsv',
        hmmer = 'outputs/trinotate/hmmscan.out'
    output:
        touch(TRINOTATE_LOAD)
    conda:
        "envs/trinotate.yaml"
    log:
        "logs/trinotate/trinotate_load.txt"
    shell:
        """
        Trinotate {input.sqlite} \
            LOAD_swissprot_blastp \
            {input.blastp_sprot} 1> {log} 2>&1
        Trinotate {input.sqlite} \
            LOAD_swissprot_blastx \
            {input.blastx_sprot} 1>> {log} 2>&1
        Trinotate {input.sqlite} \
            LOAD_custom_blast \
            --outfmt6 {input.blastp_trembl} \
            --prog blastp \
            --dbtype blastp_trembl 1>> {log} 2>&1
        Trinotate {input.sqlite} \
            LOAD_custom_blast \
            --outfmt6 {input.blastx_trembl} \
            --prog blastx \
            --dbtype blastx_trembl 1>> {log} 2>&1
        Trinotate {input.sqlite} \
            LOAD_custom_blast \
            --outfmt6 {input.blastp_uniref} \
            --prog blastp \
            --dbtype blastp_uniref90 1>> {log} 2>&1
        Trinotate {input.sqlite} \
            LOAD_custom_blast \
            --outfmt6 {input.blastx_uniref} \
            --prog blastx \
            --dbtype blastx_uniref90 1>> {log} 2>&1
        Trinotate {input.sqlite} \
            LOAD_custom_blast \
            --outfmt6 {input.blastp_nr} \
            --prog blastp \
            --dbtype blastp_nr 1>> {log} 2>&1
        Trinotate {input.sqlite} \
            LOAD_custom_blast \
            --outfmt6 {input.blastx_nr} \
            --prog blastx \
            --dbtype blastx_nr 1>> {log} 2>&1
        Trinotate {input.sqlite} \
            LOAD_pfam \
            {input.hmmer} 1>> {log} 2>&1
        """

rule trinotate_report:
    input:
        load_done = TRINOTATE_LOAD,
        sqlite = ancient('db/Trinotate.sqlite') 
    output:
        trinotate_out = TRINOTATE_ANNOT
    conda:
        "envs/trinotate.yaml"
    log:
        "logs/trinotate/trinotate_annot.txt"
    shell:
        """
        Trinotate {input.sqlite} \
            report \
            -E 1e-5 \
            > {output.trinotate_out} 2> {log}
        """

rule trinotate_stats:
    input:
        TRINOTATE_ANNOT
    output:
        TRINOTATE_STATS
    conda:
        "envs/trinotate.yaml"
    log:
        "logs/trinotate/trinotate_stats.txt"
    shell:
        """
        count_table_fields.pl {input} > {output} 2> {log}
        """


#rule clean: #clean up
#    shell:
#        """
#        rm -rf outputs/*
#        rm -rf logs/*
#        rm -rf cat_seq/
#        """