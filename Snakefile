# Snakefile modifié pour accepter soit un FASTQ, soit un BAM
configfile: "config/config.yaml"

rule all:
    input:
        expand("results/variants/{sample}.snvs.vcf", sample=config["samples"]),
        expand("results/variants/{sample}.svs.vcf", sample=config["samples"]),
        "results/report.csv"

rule prepare_bam:
    input:
        fastq="data/raw/{sample}.fastq",
        reference=config["reference"]
    output:
        bam="results/bam/{sample}.sorted.bam"
    log:
        "logs/{sample}.align.log"
    shell:
        """
        minimap2 -ax map-ont {input.reference} {input.fastq} | \
        samtools sort -o {output.bam} && \
        samtools index {output.bam}
        """

rule handle_input:
    input:
        fastq="data/raw/{sample}.fastq",
        bam="data/raw/{sample}.bam"  # BAM pré-aligné (optionnel)
    output:
        bam="results/bam/{sample}.sorted.bam"
    log:
        "logs/{sample}.input.log"
    run:
        import os
        # Si le BAM est fourni, l'utiliser tel quel
        if os.path.exists(input.bam):
            shell(f"ln -s {input.bam} {output.bam}; ln -s {input.bam}.bai {output.bam}.bai")
        else:
            # Sinon, aligner le FASTQ pour produire le BAM
            shell(f"""
                minimap2 -ax map-ont {config['reference']} {input.fastq} | \
                samtools sort -o {output.bam} && \
                samtools index {output.bam}
            """)

rule call_snvs:
    input:
        bam="results/bam/{sample}.sorted.bam",
        reference=config["reference"]
    output:
        vcf="results/variants/{sample}.snvs.vcf"
    log:
        "logs/{sample}.snvs.log"
    shell:
        """
        deepvariant --model_type=ONT \
                    --ref {input.reference} \
                    --reads {input.bam} \
                    --output_vcf {output.vcf}
        """

rule call_svs:
    input:
        bam="results/bam/{sample}.sorted.bam",
        reference=config["reference"]
    output:
        vcf="results/variants/{sample}.svs.vcf"
    log:
        "logs/{sample}.svs.log"
    shell:
        """
        sniffles --input {input.bam} --reference {input.reference} --vcf {output.vcf}
        """

rule summarize:
    input:
        snvs=expand("results/variants/{sample}.snvs.vcf", sample=config["samples"]),
        svs=expand("results/variants/{sample}.svs.vcf", sample=config["samples"])
    output:
        "results/report.csv"
    run:
        import pandas as pd
        summary = []
        for snvs, svs in zip(input.snvs, input.svs):
            summary.append({"sample": snvs.split("/")[-1].split(".")[0],
                            "snvs_vcf": snvs,
                            "svs_vcf": svs})
        df = pd.DataFrame(summary)
        df.to_csv(output[0], index=False)
