//do not change these
params.input_dir="none"
params.help = false
params.step= false

//main script
if(params.help){
    println "SMELT: a Slurm wrapper for MELT"
    println "Usage: nextflow smelt.nf --step <STEP> --working_dir <working_dir> --input_dir <input_dir> --te <TE_zip_path> --ref {reference_fasta} -c smelt.conf"
    println "--step parameters"
    println "prep - run the MELT prep procedure on the bam files in the input_dir (step1)"
    println "indiv - run the indiv analysis on the extracted signals (step2)"
    println "group  - run the group step of the MELT split in input_dir (step3)"
    println "genotype  - run the MELT split analysis on the extracted signals in input_dir (step4)"
    println "makevcf  - run the MELT split analysis on the extracted signals in input_dir (step5)"
    println "single - run melt single on all vcf files in a folder"
    println ""
    println "note! do not change the working_dir and input dir inbetween the analysis steps!"
    println ""
    println "other parameters:"
    println "--te - path to the transposable element zip file"

}else if (params.step == "prep"){

    input_dir=file(params.input_dir)
    if(!input_dir.exists()) exit 1, "Error: The input directory was not found"

    tmp="mkdir -p ${params.working_dir}/prep".execute().text
    bam_files=Channel.fromPath("${params.input_dir}/*.cram").map{
        line ->
        ["${file(line).baseName}".replaceFirst(/.cram/,""),file(line),file("${file(line)}.crai"),line]
    }
    Channel.fromPath("${params.input_dir}/*.cram").subscribe{
        println it
        "ln -s ${it} ${params.working_dir}/prep".execute().text
        "ln -s ${it}.bai ${params.working_dir}/prep".execute().text
    }

    process prep{

        publishDir "${params.working_dir}/prep", mode: 'copy', overwrite: true
        scratch false
        errorStrategy 'ignore'

   
        cpus 1
       

        input:
            set ID, file(bam), file(bai), bam_path from bam_files

        output:
            set "${bam.baseName}.cram.disc", "${bam.baseName}.cram.disc.bai", "${bam.baseName}.cram.fq" into prep_output 


        script:
        """
            java -Xmx8G -jar ${params.melt} Preprocess -bamfile ${bam} -h ${params.ref}
        """
    }

}else if (params.step == "indiv"){

    input_dir=file(params.input_dir)
    if(!input_dir.exists()) exit 1, "Error: The input directory was not found"

    bam_files=Channel.fromPath("${params.working_dir}/prep/*.cram").map{
        line ->
        ["${file(line).baseName}",file("${file(line)}.disc"),file("${file(line)}.disc.bai"),file("${file(line)}.fq"),file(line),file("${file(line)}.bai")]

    }

    process indiv{

        publishDir "${params.working_dir}/indiv", mode: 'copy', overwrite: true
        scratch false
        errorStrategy 'ignore'

        cpus 1
        memory '6 GB'

        input:
            set ID, file(disc_bam), file(disc_index), file(fq), file(bam), file(bai) from bam_files

        output:
            file "*.final.sorted.bam" into final_sorted
            file "*.final.sorted.bam.bai" into final_sorted_bai
            file "*.hum_breaks.sorted.bam" into hum_breaks
            file "*.hum_breaks.sorted.bam.bai" into hum_breaks_bai
            file "*.pulled.sorted.bam" into pulled_sorted
            file "*.pulled.sorted.bam.bai" into pulled_sorted_bai
            file "*.tmp.bed" into temp_bed

        script:
        """
            mkdir tmp_indiv_out
            java -Xmx8G -jar ${params.melt} IndivAnalysis -bamfile ${bam} -b ${params.exclude} -h ${params.ref} -t ${params.te} -c ${params.cov} -r ${params.readlen} -w tmp_indiv_out
            cp -r tmp_indiv_out/* .
        """
    }

}else if (params.step == "group"){

    input_dir=file("${params.working_dir}/indiv")
    if(!input_dir.exists()) exit 1, "Error: The input directory was not found"

    process group{

        publishDir "${params.working_dir}/group", mode: 'copy', overwrite: true
        scratch false
        errorStrategy 'ignore'

        cpus 1
       
        
        input:
            file input_dir
        
        output:
            file "*" into outs

        script:


        """
        mkdir Group
        java -Xmx8G -jar ${params.melt} GroupAnalysis -h ${params.ref} -discoverydir ${input_dir} -t ${params.te} -w Group -n ${params.genes} -r ${params.readlen}
        mv Group/* .
        rm -rf tmp

        """
    }

}else if (params.step == "genotype"){

    input_dir=file(params.input_dir)
    if(!input_dir.exists()) exit 1, "Error: The input directory was not found"

    bam_files=Channel.fromPath("${params.input_dir}/*.cram").map{
        line ->
        ["${file(line).baseName}".replaceFirst(/.cram/,""),file(line),file("${file(line)}.crai"),line]
    }

    process genotype{

        publishDir "${params.working_dir}/genotype", mode: 'copy', overwrite: true
        scratch false
        errorStrategy 'ignore'


        cpus 1
        memory '6 GB'

        input:
            set ID, file(bam),file(bai), bam_path from bam_files

        output:
            file "*.tsv" into genotype_output 


        script:
        """
            mkdir genotype_dir
            java -Xmx8G -jar ${params.melt} Genotype -h ${params.ref} -bamfile ${bam} -t ${params.te} -w genotype_dir -p ${params.working_dir}/group/ -e ${params.insert_size}
            mv genotype_dir/*tsv .
        """

    }
}else if (params.step == "makevcf"){

    input_dir=file("${params.working_dir}/prep")
    if(!input_dir.exists()) exit 1, "Error: The input directory was not found"

    genotype_dir=file("${params.working_dir}/genotype/")
    if(!genotype_dir.exists()) exit 1, "Error: The input genotype directory was not found (did you run the genotype step?)"


    group_dir=file("${params.working_dir}/group")
    if(!group_dir.exists()) exit 1, "Error: The input genotype directory was not found (did you run the group step?)"

    process makevcf{

        publishDir "${params.working_dir}/", mode: 'copy', overwrite: true
        scratch false
        errorStrategy 'ignore'

        cpus 1
       

        input:
            file genotype_dir
            file group_dir

        output:
            file "*.vcf" into final_vcf 


        script:

        """
        printf "" > tsv_list.txt

        for tsvfile in ${params.working_dir}/genotype/*.tsv
            do
            printf "%s\n" \$tsvfile  >> tsv_list.txt
        done

        mkdir vcf_output
        java -Xmx8G -jar ${params.melt} MakeVCF -ac -h ${params.ref} -genotypingdir ${params.working_dir}/genotype -t ${params.te} -w ${genotype_dir} -p ${group_dir} -o vcf_output
        mv vcf_output/*vcf .
        """

    }
}else if (params.step == "single"){

    input_dir=file(params.input_dir)
    if(!input_dir.exists()) exit 1, "Error: The input directory was not found"

    bam=Channel.fromPath("${params.input_dir}/*.bam").map{
        line ->
        ["${file(line).baseName}".replaceFirst(/.bam/,""),file(line),file("${file(line)}.bai"),line]
    }



    process single{
       publishDir "${params.working_dir}/single", mode: 'copy', overwrite: true
     
       cpus 1
      

       input:
            set ID, file(bam), file(bai), bam_path from bam

       output:
             file "${bam.baseName}.vcf" into single_vcf_folder
       
       """
       java -Xmx8G -jar ${params.melt}  Single -h ${params.ref} -b ${params.exclude} -bamfile ${bam} -n ${params.genes} -t ${params.te} -c ${params.cov} -r ${params.readlen} -e ${params.insert_size} -w ${bam.baseName}
       vcf-concat ${bam.baseName}/*.vcf | vcf-sort -c > ${bam.baseName}.vcf
       """


    }

}else{
    println  "type: nextflow smelt.nf --help"
}
