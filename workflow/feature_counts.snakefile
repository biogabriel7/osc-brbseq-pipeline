from snakemake.utils import min_version
min_version("6.0")

# Define SAMPLES as an empty dict if not already defined (for syntax checking)
if 'SAMPLES' not in globals():
    SAMPLES = {}
    config = {}

# No configfile import - config is loaded in main Snakefile

# Rule to run feature counts for all samples
rule count_all:
    input:
        # Individual count files
        expand("Analysis/Counts/FeatureCounts/{sample}/{sample}.counts.txt", 
               sample=SAMPLES.keys()),
        # Merged count matrix
        "Analysis/Counts/FeatureCounts/merged_gene_counts.txt",
        # MultiQC report
        "Analysis/Counts/MultiQC/multiqc_report.html",
        # Completion marker
        "Analysis/Counts/.counts_complete"

# Rule to run featureCounts on each BAM file
rule feature_counts:
    input:
        bam="Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bai="Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai",
        alignment_complete="Analysis/Alignment/.main_alignment_complete"
    output:
        counts="Analysis/Counts/FeatureCounts/{sample}/{sample}.counts.txt",
        summary="Analysis/Counts/FeatureCounts/{sample}/{sample}.counts.txt.summary"
    log:
        "logs/counts/{sample}.featurecounts.log"
    params:
        gtf=config.get("gtf_file", "resources/genome/Homo_sapiens.GRCh38.90.gtf"),
        # Feature type to count (default: exon)
        feature_type=config.get("feature_type", "exon"),
        # Attribute to use for grouping features (default: gene_id)
        attribute=config.get("attribute", "gene_id"),
        # Strand specificity: 0 (unstranded), 1 (stranded), 2 (reversely stranded)
        strand=config.get("strand_specificity", "0"),
        # Additional parameters
        extra_params=config.get("featurecounts_extra", "-p -B -C --fracOverlap 0.2")
    threads: config.get("featurecounts_threads", 4)
    resources:
        mem_mb=config.get("featurecounts_memory", 8000)
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.counts})
        
        # Run featureCounts
        echo "Running featureCounts on {wildcards.sample}..." > {log} 2>&1
        
        featureCounts \
            -a {params.gtf} \
            -o {output.counts} \
            -t {params.feature_type} \
            -g {params.attribute} \
            -s {params.strand} \
            -T {threads} \
            {params.extra_params} \
            {input.bam} \
            >> {log} 2>&1
            
        # Check if output was created successfully
        if [ ! -f "{output.counts}" ]; then
            echo "ERROR: featureCounts failed to create output file" >> {log}
            exit 1
        fi
        
        # Print summary statistics
        echo "Summary statistics:" >> {log}
        cat {output.summary} >> {log}
        """

# Rule to merge individual count files into a single matrix
rule merge_counts:
    input:
        counts=expand("Analysis/Counts/FeatureCounts/{sample}/{sample}.counts.txt", 
                     sample=SAMPLES.keys())
    output:
        merged="Analysis/Counts/FeatureCounts/merged_gene_counts.txt",
        normalized="Analysis/Counts/FeatureCounts/merged_gene_counts.normalized.txt"
    log:
        "logs/counts/merge_counts.log"
    resources:
        mem_mb=config.get("merge_counts_memory", 4000)
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.merged})
        
        # Run the merge script
        echo "Merging count files..." > {log} 2>&1
        
        # Create header with sample names
        echo -n "Gene_ID" > {output.merged}
        for f in {input.counts}; do
            sample=$(basename $(dirname $f))
            echo -n -e "\\t$sample" >> {output.merged}
        done
        echo "" >> {output.merged}
        
        # Extract the counts (skip header lines)
        # Assuming the first column is gene ID and the 7th column is the count
        # This works for standard featureCounts output format
        
        # Get the list of all genes from the first file (skip header lines)
        tail -n +3 $(echo {input.counts} | cut -d' ' -f1) | cut -f1 > gene_list.tmp
        
        # For each gene, extract counts from all files
        while read gene; do
            echo -n "$gene" >> {output.merged}
            for f in {input.counts}; do
                count=$(tail -n +3 $f | grep -w "^$gene" | cut -f7)
                if [ -z "$count" ]; then
                    count="0"
                fi
                echo -n -e "\\t$count" >> {output.merged}
            done
            echo "" >> {output.merged}
        done < gene_list.tmp
        
        # Remove temporary file
        rm gene_list.tmp
        
        # Create normalized counts (CPM)
        echo "Generating normalized counts (CPM)..." >> {log}
        
        # Copy header
        head -n 1 {output.merged} > {output.normalized}
        
        # Calculate library sizes
        declare -A lib_sizes
        for f in {input.counts}; do
            sample=$(basename $(dirname $f))
            lib_size=$(tail -n +3 $f | awk '{{sum+=$7}} END {{print sum}}')
            lib_sizes["$sample"]=$lib_size
            echo "Library size for $sample: $lib_size" >> {log}
        done
        
        # Normalize counts (CPM)
        tail -n +2 {output.merged} | while read line; do
            gene=$(echo "$line" | cut -f1)
            echo -n "$gene" >> {output.normalized}
            
            # Get counts for each sample and normalize
            for sample in $(head -n 1 {output.merged} | cut -f2-); do
                count=$(echo "$line" | grep -w "^$gene" | cut -f$(head -n 1 {output.merged} | tr '\\t' '\\n' | grep -n "^$sample$" | cut -d: -f1))
                lib_size=${lib_sizes["$sample"]}
                
                # Calculate CPM: (count * 1,000,000) / library_size
                cpm=$(awk -v count="$count" -v lib="$lib_size" 'BEGIN {{printf "%.2f", (count * 1000000) / lib}}')
                echo -n -e "\\t$cpm" >> {output.normalized}
            done
            echo "" >> {output.normalized}
        done
        
        echo "Count merging completed successfully" >> {log}
        echo "Output files:" >> {log}
        echo "  Raw counts: {output.merged}" >> {log}
        echo "  Normalized counts (CPM): {output.normalized}" >> {log}
        """

# Rule to run MultiQC on feature counts results
rule multiqc_counts:
    input:
        summaries=expand("Analysis/Counts/FeatureCounts/{sample}/{sample}.counts.txt.summary", 
                        sample=SAMPLES.keys())
    output:
        html="Analysis/Counts/MultiQC/multiqc_report.html",
        data_dir=directory("Analysis/Counts/MultiQC/multiqc_data")
    log:
        "logs/counts/multiqc.log"
    params:
        outdir="Analysis/Counts/MultiQC"
    resources:
        mem_mb=config.get("multiqc_counts_memory", 4000)
    shell:
        """
        # Create output directory
        mkdir -p {params.outdir}
        
        # Run MultiQC
        echo "Running MultiQC on feature counts results..." > {log} 2>&1
        
        multiqc --force \
            --outdir {params.outdir} \
            Analysis/Counts/FeatureCounts/ \
            >> {log} 2>&1
            
        # Verify output
        if [ ! -f "{output.html}" ]; then
            echo "ERROR: MultiQC report was not created" >> {log}
            exit 1
        fi
        
        echo "MultiQC completed successfully" >> {log}
        """

# Rule to mark completion of counting stage
rule counts_complete:
    input:
        counts=expand("Analysis/Counts/FeatureCounts/{sample}/{sample}.counts.txt", 
                     sample=SAMPLES.keys()),
        merged="Analysis/Counts/FeatureCounts/merged_gene_counts.txt",
        normalized="Analysis/Counts/FeatureCounts/merged_gene_counts.normalized.txt",
        multiqc="Analysis/Counts/MultiQC/multiqc_report.html"
    output:
        touch("Analysis/Counts/.counts_complete")
    log:
        "logs/workflow/counts_complete.log"
    shell:
        """
        echo "Feature counts completed successfully at $(date)" > {log}
        echo "Generated count files for $(echo {input.counts} | wc -w) samples" >> {log}
        echo "Merged count matrix: {input.merged}" >> {log}
        echo "Normalized count matrix: {input.normalized}" >> {log}
        """ 