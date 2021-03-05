configfile: "config/config.yml"

rule all:
    input:
        expand("results/VirSorter/{sample}/VIRSorter_global-phage-signal.csv", sample=config["samples"]),
        expand("results/VirSorter2/{sample}/final-viral-score.tsv", sample=config["samples"]),
        expand("results/VirFinder/{sample}.txt", sample=config["samples"]),
        expand("results/DeepVirFinder/{sample}_dvf.txt", sample=config["samples"]),
        expand("results/seeker/{sample}_seeker.txt", sample=config["samples"]),
        expand("results/Vibrant/{sample}_renamed/VIBRANT_{sample}_renamed/VIBRANT_log_{sample}_renamed.log", sample=config["samples"]),
        expand("input/{sample}_DEPTH.txt", sample=config["samples"]),
        expand("results/marvel/{sample}_bin/results/phage_bins/list_bins.txt", sample=config["samples"]),
        expand("results/viral_sequences/prophages/{sample}_prophages_selection1.fna", sample=config["samples"]),
        expand("results/viral_sequences/phages/{sample}_phages_selection1.fna", sample=config["samples"]),


##### workflow starts here

rule virsorter:
    input:
        f="input/{base}_renamed.fasta",
        scr="scripts/VirSorter/wrapper_phage_contigs_sorter_iPlant.pl",
    params:
        outdir= "results/VirSorter/{base}",
        datadir= "/xdisk/bhurwitz/mig2020/rsgrps/bhurwitz/alise/tools/virsorter-data",
    output:
        "results/VirSorter/{base}/VIRSorter_global-phage-signal.csv",
    shell:
        """
        set +eu
        source ~/.bashrc
        source activate viral_env
        {input.scr} -f {input.f} --db 1 --wdir {params.outdir} --ncpu 20 --data-dir {params.datadir}
        """

rule virsorter2:
    input:
        f="input/{base}_renamed.fasta",    
    params:
        outdir= "results/VirSorter2/{base}",
    output:
        "results/VirSorter2/{base}/final-viral-score.tsv",
    shell:
        """
        set +eu
        source ~/.bashrc
        source activate virSorter2
        virsorter run -w {params.outdir} -i {input.f} -j 20 all 
        """

rule deepvirFinder:
    input:
        f="input/{base}_renamed.fasta",
    params:
        temp="results/DeepVirFinder/{base}_renamed.fasta_gt500bp_dvfpred.txt",
        scr="scripts/DeepVirFinder/dvf.py",
        outdir="results/DeepVirFinder",
        outfile="results/DeepVirFinder/{base}_dvf.txt",
    output:
        "results/DeepVirFinder/{base}_dvf.txt",
    shell:
        """
        set +eu
        source ~/.bashrc 
        source activate dvf
        python {params.scr} -i {input.f} -o {params.outdir} -l 500
        mv {params.temp} {params.outfile}
        """

rule seeker:
    input:
        f="input/{base}_renamed.fasta",
    params:
        scr="scripts/seeker/filter_size.py",
        filt="input/{base}_filtrenamed.fasta",
        outdir= "results/seeker",
    output:
        "results/seeker/{base}_seeker.txt",
    shell:
        """
        set +eu
        source ~/.bashrc
        source activate seeker
        python {params.scr} {input.f} {params.filt} 200
        predict-metagenome {params.filt} >> {output}
        """


rule virfinder:
    input:
        f="input/{base}_renamed.fasta",
        scr="scripts/VirFinder/eval_default.r",
    output:
        "results/VirFinder/{base}.txt",
    shell:
        """
        set +eu
        source ~/.bashrc
        module load R
        Rscript {input.scr} {input.f} {output}
        """

rule vibrant:
    input:
        f="input/{base}_renamed.fasta",
    params:
        outdir="results/Vibrant/{base}_renamed",
    output:
        "results/Vibrant/{base}_renamed/VIBRANT_{base}_renamed/VIBRANT_log_{base}_renamed.log",
    shell:
        """
        set +eu
        source ~/.bashrc
        source activate viral_env
        cd {params.outdir}
        VIBRANT_run.py -i ../../../{input.f}
        cd ../../..
        """

rule metabat:
    input:
        f="input/{base}_renamed.fasta",
        bam="input/{base}.bam"
    params:
        bin_dir="results/marvel/{base}_bin",
    output:
        "input/{base}_DEPTH.txt",
    shell:
        """
        source ~/.bashrc
        set +eu
        source activate viral_env
        mkdir -p {params.bin_dir}
        jgi_summarize_bam_contig_depths --outputDepth {output} {input.bam}
        metabat2 -i {input.f} -a {output} -m 1500 -s 10000 -o {params.bin_dir}/bin
        """


rule marvel:
    input:
        f="input/{base}_DEPTH.txt",
    params:
        marvel_script="scripts/MARVEL",
        bin_dir="results/marvel/{base}_bin",
    output:
        "results/marvel/{base}_bin/results/phage_bins/list_bins.txt"
    shell:
        """
        source ~/.bashrc
        set +eu
        source activate viral_env
        cd {params.marvel_script}
        python3 marvel_bins.py -i ../../{params.bin_dir} -t 28
        cd ../.. 
        find {params.bin_dir}/results/phage_bins -name "*.fasta" >> {output}
        """

rule prophage:
    input:
        f_marvel="results/marvel/{base}_bin/results/phage_bins/list_bins.txt",
        f_virsorter="results/VirSorter/{base}/VIRSorter_global-phage-signal.csv",
        f_virfinder="results/VirFinder/{base}.txt",
        f_vibrant="results/Vibrant/{base}_renamed/VIBRANT_{base}_renamed/VIBRANT_log_{base}_renamed.log",
    params:
        prophage_script="scripts/parsing_results/get_prophages.sh",
        output_dir="results/viral_sequences/prophages",
        sample_name="{base}",
    output:
        "results/viral_sequences/prophages/{base}_prophages_selection1.fna"
    shell:
        """
        mkdir -p {params.output_dir}
        bash {params.prophage_script} {params.sample_name} {params.output_dir} 
        """

rule phage:
    input:
        f_marvel="results/marvel/{base}_bin/results/phage_bins/list_bins.txt",
        f_virsorter="results/VirSorter/{base}/VIRSorter_global-phage-signal.csv",
        f_virfinder="results/VirFinder/{base}.txt",
        f_vibrant="results/Vibrant/{base}_renamed/VIBRANT_{base}_renamed/VIBRANT_log_{base}_renamed.log",
    params:
        phage_script="scripts/parsing_results/get_viral_sequences.sh",
        output_dir="results/viral_sequences/phages",
        sample_name="{base}",
    output:
        "results/viral_sequences/phages/{base}_phages_selection1.fna"
    shell:
        """
        mkdir -p {params.output_dir}
        bash {params.phage_script} {params.sample_name} {params.output_dir}
        """




