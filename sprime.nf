#!/usr/bin/env nextflow
/* Pipeline to identify and characterize introgression tracts from archaic  *
 *  hominins (e.g. Neandertals and Denisovans)                              *
 * Core steps:                                                              *
 *  Extraction of target populations + outgroup from main VCF for Sprime -> *
 *  Concatenation of per-autosome VCFs for Sprime ->                        *
 *  Sprime for each target population ->                                    *
 *  (Calculate match rates && Project S' tracts onto individuals )          */

nextflow.enable.dsl=2

//Default paths, globs, and regexes:
//Jointly genotyped VCFs:
params.vcf_glob = "${projectDir}/perchrom_VCFs/*.vcf.gz"
//Regex for extracting chromosome from VCF filename:
params.vcf_regex = ~/_chr(\p{Alnum}+)$/
//Recombination rate maps:
params.recmap_glob = "${projectDir}/"
//Regex for extracting chromosome from recombination rate map filename:
params.recmap_regex = ~/chr([0-9XY_par]+)/
//Archaic VCFs:
params.arcvcf_glob = "${projectDir}/archaic_VCFs/*.vcf.gz"

//Include/filter expression to apply to input VCFs:
params.input_filter_str = ""
//Sample ID file for exclusion:
params.samples_to_exclude = "<(echo)"

//Reference-related parameters for the pipeline:
params.ref_prefix = "/gpfs/gibbs/pi/tucci/pfr8/refs"
params.ref = "${params.ref_prefix}/1kGP/hs37d5/hs37d5.fa"
params.autosomes = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"
params.nuclearchroms = "${params.autosomes},X,Y"

//Default parameter values:
//Sample ID column name in metadata file:
params.id_colname = "Sample"
//Outgroup name to use:
params.sprime_outgroup_colname = "Region"
params.sprime_outgroup = "Africa"
//Target group column name:
params.sprime_target_colname = "AnalysisGroup"

//Archaic hominin groups represented in archaic VCFs:
params.archaic_groups = "Neandertal,Denisovan"
//Archaic hominin samples to consider individually in archaic VCFs:
params.archaic_indivs = "AltaiNeandertal,Vindija33.19,Chagyrskaya-Phalanx"
//Short names for these archaic hominin samples:
params.archaic_shortnames = "Altai,Vindija,Chagyrskaya"

//Internal gap length in phased projection allowed to be spanned:
//This gap length is a count of Sprime sites at which the haplotype
// in question does not match the Sprime-identified putative archaic
// allele.
//Strongly recommend setting to 0, as this matches haplotype
// predictions from other methods better.
params.tract_max_gap = 0

//Defaults for cpus, memory, and time for each process:
//VCF subsetting for Sprime
params.vcfsubset_cpus = 1
params.vcfsubset_mem = 4
params.vcfsubset_mem_increment = 8
params.vcfsubset_timeout = '24h'
params.vcfsubset_retry_timeout = '24h'
params.vcfsubset_queue = 'day'
params.vcfsubset_retry_queue = 'day'
//VCF concatenating for Sprime
params.concatvcf_cpus = 1
params.concatvcf_mem = 4
params.concatvcf_mem_increment = 16
params.concatvcf_timeout = '24h'
params.concatvcf_retry_timeout = '48h'
params.concatvcf_queue = 'day'
params.concatvcf_retry_queue = 'week'
//Sprime
params.sprime_cpus = 1
params.sprime_mem = 8
params.sprime_mem_increment = 16
params.sprime_timeout = '24h'
params.sprime_retry_timeout = '24h'
params.sprime_queue = 'day'
params.sprime_retry_queue = 'day'
//Sprime match rates
params.matchrate_cpus = 1
params.matchrate_mem = 4
params.matchrate_mem_increment = 16
params.matchrate_timeout = '24h'
params.matchrate_retry_timeout = '24h'
params.matchrate_queue = 'day'
params.matchrate_retry_queue = 'day'
//Sprime projection
params.projection_cpus = 1
params.projection_mem = 1
params.projection_mem_increment = 16
params.projection_timeout = '24h'
params.projection_retry_timeout = '48h'
params.projection_queue = 'day'
params.projection_retry_queue = 'week'
//Combine outputs across targets and chromosomes
params.catouts_cpus = 1
params.catouts_mem = 1
params.catouts_mem_increment = 16
params.catouts_timeout = '2h'
params.catouts_retry_timeout = '24h'
params.catouts_queue = 'day'
params.catouts_retry_queue = 'day'

process vcfsubset {
   tag "${pop} chr${chrom}"

   cpus params.vcfsubset_cpus
   memory { params.vcfsubset_mem.plus(task.attempt.minus(1).multiply(params.vcfsubset_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.vcfsubset_retry_timeout : params.vcfsubset_timeout }
   queue { task.attempt >= 2 ? params.vcfsubset_retry_queue : params.vcfsubset_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*_samples.tsv', saveAs: { 'chr'+chrom+'_'+it }

   input:
   tuple val(pop), val(chrom), path(input_vcf), path(input_tbi)
   path metadata
   path excluded_samples

   output:
   tuple path("bcftools_view_selectpops_Sprime_${pop}_chr${chrom}.stderr"), path("bcftools_view_nomissinggenos_Sprime_${pop}_chr${chrom}.stderr"), path("bcftools_view_nomissinggenos_Sprime_${pop}_chr${chrom}.stdout"), emit: logs
   tuple val(pop), path("${pop}_chr${chrom}_nomissinggenos.vcf.gz"), path("${pop}_chr${chrom}_nomissinggenos.vcf.gz.tbi"), emit: vcf
   tuple val(pop), path("${pop}_outgroup_samples.tsv"), emit: outgroup
   tuple val(pop), path("${pop}_samples.tsv"), emit: samples

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.sprime_outgroup_colname}" -v "select=!{params.sprime_outgroup}" !{metadata} | \
      !{projectDir}/HumanPopGenScripts/excludeSamples.awk !{excluded_samples} - > !{pop}_outgroup_samples.tsv
   !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.sprime_target_colname}" -v "select=!{pop}" !{metadata} | \
      !{projectDir}/HumanPopGenScripts/excludeSamples.awk !{excluded_samples} - > !{pop}_samples.tsv
   cat !{pop}_samples.tsv !{pop}_outgroup_samples.tsv > !{pop}_plusOutgroup_samples.tsv
   #We have to separate the sample selection from the genotype missingness filter
   # because of a bug in bcftools view when -S and -g are used in tandem:
   #Way too many sites get filtered if you use them together, possibly more than
   # if -g was used alone.
   bcftools view -S !{pop}_plusOutgroup_samples.tsv -Ou !{input_vcf} 2> bcftools_view_selectpops_Sprime_!{pop}_chr!{chrom}.stderr | \
   bcftools view -i '!{params.input_filter_str}' -g ^miss -Oz -o !{pop}_chr!{chrom}_nomissinggenos.vcf.gz 2> bcftools_view_nomissinggenos_Sprime_!{pop}_chr!{chrom}.stderr > bcftools_view_nomissinggenos_Sprime_!{pop}_chr!{chrom}.stdout
   tabix -f !{pop}_chr!{chrom}_nomissinggenos.vcf.gz
   '''
}

process concatvcf {
   tag "${pop}"

   cpus params.concatvcf_cpus
   memory { params.concatvcf_mem.plus(task.attempt.minus(1).multiply(params.concatvcf_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.concatvcf_retry_timeout : params.concatvcf_timeout }
   queue { task.attempt >= 2 ? params.concatvcf_retry_queue : params.concatvcf_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'

   input:
   tuple val(pop), path("*"), path("*")
   val(autosome_num_list)
   val(num_autosomes)

   output:
   tuple path("bcftools_concat_Sprime_${pop}.stderr"), path("bcftools_concat_Sprime_${pop}.stdout"), emit: logs
   tuple val(pop), path("${pop}_nomissinggenos.vcf.gz"), path("${pop}_nomissinggenos.vcf.gz.tbi"), emit: vcf

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_htslib}
   if [[ "!{num_autosomes}" -gt "1" ]]; then
      bcftools concat --threads !{task.cpus} -Oz -o !{pop}_nomissinggenos.vcf.gz !{pop}_chr{!{autosome_num_list}}_nomissinggenos.vcf.gz 2> bcftools_concat_Sprime_!{pop}.stderr > bcftools_concat_Sprime_!{pop}.stdout
   else
      ln -s !{pop}_chr!{autosome_num_list}_nomissinggenos.vcf.gz !{pop}_nomissinggenos.vcf.gz
      echo "Only one chromosome provided, skipping bcftools concat for !{pop}" > bcftools_concat_Sprime_!{pop}.stderr
      echo "" > bcftools_concat_Sprime_!{pop}.stdout
   fi
   tabix -f !{pop}_nomissinggenos.vcf.gz
   '''
}

process sprime {
   tag "${pop} chr${chrom}"

   cpus params.sprime_cpus
   memory { params.sprime_mem.plus(1).plus(task.attempt.minus(1).multiply(params.sprime_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.sprime_retry_timeout : params.sprime_timeout }
   queue { task.attempt >= 2 ? params.sprime_retry_queue : params.sprime_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.std{err,out}'
   publishDir path: "${params.output_dir}/logs", mode: 'copy', pattern: '*.log'
   publishDir path: "${params.output_dir}/Sprime", mode: 'copy', pattern: '*.score'

   input:
   tuple val(pop), path("${pop}_nomissinggenos.vcf.gz"), path("${pop}_nomissinggenos.vcf.gz.tbi"), path("${pop}_outgroup_samples.tsv"), val(chrom), path(recmap)
   val(autosome_map)

   output:
   tuple path("Sprime_${pop}_chr${chrom}.stderr"), path("Sprime_${pop}_chr${chrom}.stdout"), path("${pop}_chr${chrom}_Sprime.log"), emit: logs
   tuple val(pop), path("${pop}_chr${chrom}_Sprime.score"), emit: score

   shell:
   sprime_retry_mem = params.sprime_mem.plus(task.attempt.minus(1).multiply(params.sprime_mem_increment))
   chrid = autosome_map[chrom]
   '''
   module load !{params.mod_Sprime}
   java -Xmx!{sprime_retry_mem}g -jar ${SPRIME} gt=!{pop}_nomissinggenos.vcf.gz outgroup=!{pop}_outgroup_samples.tsv map=!{recmap} chrom=!{chrid} out=!{pop}_chr!{chrom}_Sprime 2> Sprime_!{pop}_chr!{chrom}.stderr > Sprime_!{pop}_chr!{chrom}.stdout
   '''
}

process matchrate {
   tag "${pop}"

   cpus params.matchrate_cpus
   memory { params.matchrate_mem.plus(task.attempt.minus(1).multiply(params.matchrate_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.matchrate_retry_timeout : params.matchrate_timeout }
   queue { task.attempt >= 2 ? params.matchrate_retry_queue : params.matchrate_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/Sprime", mode: 'copy', pattern: '*.tsv.gz'

   input:
   tuple val(pop), path("*")
   path("*")
   path archaic_metadata
   val(autosome_num_list)

   output:
   path("${pop}_autosomes_Sprime_matches.tsv.gz"), emit: matches
   path("${pop}_autosomes_Sprime_match_rates.tsv.gz"), emit: matchrates

   shell:
   '''
   module load !{params.mod_bcftools}
   module load !{params.mod_bedtools}
   #Identify Sprime alleles that match or mismatch the archaics:
   #The extra comma is not a typo, it's a way to handle only one autosome with brace expansion.
   for chr in {!{autosome_num_list},};
      do
      #Establish the set of Sprime sites for this population and chromosome:
      fgrep -v "CHROM" !{pop}_chr${chr}_Sprime.score | \
         cut -f1,2 > !{pop}_chr${chr}_Sprime_sites.tsv
      !{projectDir}/HumanPopGenScripts/Sprime/archaicMatchSprime.awk -v "spop=!{pop}" -v "groups=!{params.archaic_groups}" -v "arcs=!{params.archaic_indivs}" -v "short=!{params.archaic_shortnames}" !{archaic_metadata} \
         <(bcftools view -R !{pop}_chr${chr}_Sprime_sites.tsv -T !{pop}_chr${chr}_Sprime_sites.tsv *_chr${chr}.vcf.gz) \
         !{pop}_chr${chr}_Sprime.score
   done | \
      awk 'BEGIN{FS="\t";OFS=FS;}/^CHROM/&&NR==1{print;}!/^CHROM/{print;}' | \
      gzip -9 > !{pop}_autosomes_Sprime_matches.tsv.gz
   #Calculate the match rates from this file:
   gzip -dc !{pop}_autosomes_Sprime_matches.tsv.gz | \
      !{projectDir}/HumanPopGenScripts/Sprime/SprimeArchaicMatchRate.awk -v "species=!{params.archaic_groups}" -v "arcs=!{params.archaic_shortnames}" | \
      gzip -9 > !{pop}_autosomes_Sprime_match_rates.tsv.gz
   '''
}

process projection {
   tag "${pop}"

   cpus params.projection_cpus
   memory { params.projection_mem.plus(task.attempt.minus(1).multiply(params.projection_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.projection_retry_timeout : params.projection_timeout }
   queue { task.attempt >= 2 ? params.projection_retry_queue : params.projection_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/Sprime", mode: 'copy', pattern: '*.bed'

   input:
   tuple val(pop), path("*"), path("*"), path("*"), path("*")
   val(autosome_list)
   val(has_chr_prefix)

   output:
   path("${pop}_Sprime_phased_tracts_perSample_maxgap${params.tract_max_gap}.bed"), emit: projections

   shell:
   '''
   module load !{params.mod_bcftools}
   #Project Sprime alleles onto each individual to identify tracts from genotypes:
   #The extra comma is not a typo, it's a way to handle only one autosome with brace expansion.
   for chr in {!{autosome_list},};
      do
      if [[ "!{has_chr_prefix}" == "true" ]]; then
         chrnum=${chr#chr};
      else
         chrnum=${chr};
      fi;
      bcftools query -f '%CHROM:%POS[\t%GT]\n' -r ${chr} -H -S ^!{pop}_outgroup_samples.tsv !{pop}_chr${chrnum}_nomissinggenos.vcf.gz | \
         !{projectDir}/HumanPopGenScripts/Sprime/SprimePerSampleTracts.awk -v "spop=!{pop}" -v "allout=1" -v "phased=1" !{pop}_chr${chrnum}_Sprime.score - | \
         !{projectDir}/HumanPopGenScripts/Sprime/SprimeTractBED.awk -v "phased=1" -v "max_gap=!{params.tract_max_gap}" | \
         sort -k1,1V -k2,2n -k3,3n
   done > !{pop}_Sprime_phased_tracts_perSample_maxgap!{params.tract_max_gap}.bed
   '''
}

process cat_outs {
   cpus params.catouts_cpus
   memory { params.catouts_mem.plus(task.attempt.minus(1).multiply(params.catouts_mem_increment))+' GB' }
   time { task.attempt >= 2 ? params.catouts_retry_timeout : params.catouts_timeout }
   queue { task.attempt >= 2 ? params.catouts_retry_queue : params.catouts_queue }
   errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
   maxRetries 1

   publishDir path: "${params.output_dir}/Sprime", mode: 'copy', pattern: '*.t{sv,csv}.gz'

   input:
   tuple path(matchrates), path(projections), path(targetpops)
   path metadata
   path excluded_samples

   output:
   tuple path("${params.run_name}_perPop_Sprime_autosomal_match_rates.tsv.gz"), path("${params.run_name}_Sprime_perChrom_perIndiv_perPop_tract_lengths.tsv.gz"), path("${params.run_name}_Sprime_allPopPairs_tract_freqs.tsv.gz"), path("${params.run_name}_Sprime_gene_name_lists.tcsv.gz"), emit: outputs

   shell:
   '''
   #Symlink the many input files to avoid SLURM SBATCH script size limits:
   while IFS=$'\t' read -a a;
      do
      ln -s ${a[2]} ${a[1]};
   done < <(cat matchrate_paths.tsv projection_paths.tsv)
   #Concatenate the match rate files across populations:
   header=1
   while IFS=$'\t' read p;
      do
      gzip -dc ${p}_autosomes_Sprime_match_rates.tsv.gz | \
         awk -v "header=${header}" 'BEGIN{FS="\t";OFS=FS;}NR==1&&header{print;}NR>1{print;}'
      header=0;
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_perPop_Sprime_autosomal_match_rates.tsv.gz
   unset header
   #Calculate total per-individual tract length and
   # concatenate files across populations:
   header=1
   while IFS=$'\t' read p;
      do
      !{projectDir}/HumanPopGenScripts/selectSubsamples.awk -v "idcol=!{params.id_colname}" -v "selectcol=!{params.sprime_target_colname}" -v "select=${p}" !{metadata} | \
         !{projectDir}/HumanPopGenScripts/excludeSamples.awk !{excluded_samples} - | \
         !{projectDir}/HumanPopGenScripts/Sprime/SprimeTractBEDtoLengths.awk -v "header=${header}" -v "pop=${p}" -v "phased=1" -v "addpop=0" - ${p}_Sprime_phased_tracts_perSample_maxgap!{params.tract_max_gap}.bed
      header="";
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_Sprime_perChrom_perIndiv_perPop_tract_lengths.tsv.gz
   unset header
   #Concatenate tract frequency files and add a header:
   header=1
   while IFS=$'\t' read p;
      do
      if [[ "${header}" == "1" ]]; then
         printf "Chromosome\tStart\tEnd\tTractID\tQueryPop\tMedianAF\tMinAF\tMaxAF\tNumSites\n"
         header=""
      fi
      while IFS=$'\t' read q;
         do
         gzip -dc ${p}_Sprime_${q}_tract_freqs.tsv.gz | \
            sort -k1,1V -k2,2n -k3,3n;
      done < !{targetpops};
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_Sprime_allPopPairs_tract_freqs.tsv.gz
   #Concatenate the gene name lists across populations:
   header=1
   while IFS=$'\t' read p;
      do
      if [[ "${header}" == "1" ]]; then
         printf "Chromosome\tStart\tEnd\tTractID\tOverlappingGenes\n"
         header=""
      fi
      gzip -dc ${p}_Sprime_gene_name_lists.tcsv.gz
      header="";
   done < !{targetpops} | \
      gzip -9 > !{params.run_name}_Sprime_gene_name_lists.tcsv.gz
   unset header
   '''
}

workflow {
   //Detection mechanism for "chr" prefixes of autosomes:
   has_chr_prefix = params.autosomes.count('chr') > 0
   autosome_list = params.autosomes
   autosome_num_list = has_chr_prefix ? params.autosomes.replaceAll('chr', '') : autosome_list
   Map autosome_map = [params.autosomes.replaceAll('chr', '').tokenize(','), params.autosomes.tokenize(',')]
      .transpose()
      .collectEntries()

   num_autosomes = params.autosomes.tokenize(',').size()

   //Set up the channels of per-chromosome jointly-genotyped VCFs and their indices:
   Channel
      .fromPath(params.vcf_glob, checkIfExists: true)
      .ifEmpty { error "Unable to find VCFs matching glob: ${params.vcf_glob}" }
      .map { a -> [ (a.getSimpleName() =~ params.vcf_regex)[0][1], a] }
      .ifEmpty { error "Regex failed to extract chromosome ID from VCF: ${params.vcf_regex}" }
      .tap { perchrom_vcfs }
      .subscribe { println "Added ${it[1]} to perchrom_vcfs channel" }

   Channel
      .fromPath(params.vcf_glob+'.tbi', checkIfExists: true)
      .ifEmpty { error "Unable to find VCF indices matching glob: ${params.vcf_glob}.tbi" }
      .map { a -> [ (a.getSimpleName() =~ params.vcf_regex)[0][1], a] }
      .ifEmpty { error "Regex failed to extract chromosome ID from VCF index: ${params.vcf_regex}" }
      .tap { perchrom_tbis }
      .subscribe { println "Added ${it[1]} to perchrom_tbis channel" }

   //Preprocess the per-chromosome VCF channel to include the indices:
   perchrom_vcfs_subset = perchrom_vcfs
      .join(perchrom_tbis, by: 0, failOnDuplicate: true, failOnMismatch: true)
      .filter({ autosome_num_list.tokenize(',').contains(it[0]) })

   //Set up the channels of archaic VCFs and their indices:
   Channel
      .fromPath(params.arcvcf_glob, checkIfExists: true)
      .ifEmpty { error "Unable to find archaic VCFs matching glob: ${params.arcvcf_glob}" }
      .map { a -> [ (a.getSimpleName() =~ params.vcf_regex)[0][1], a] }
      .ifEmpty { error "Regex failed to extract chromosome ID from archaic VCF: ${params.vcf_regex}" }
      .tap { archaic_vcfs }
      .subscribe { println "Added ${it[1]} to archaic_vcfs channel" }

   Channel
      .fromPath(params.arcvcf_glob+'.tbi', checkIfExists: true)
      .ifEmpty { error "Unable to find archaic VCF indices matching glob: ${params.arcvcf_glob}.tbi" }
      .map { a -> [ (a.getSimpleName() =~ params.vcf_regex)[0][1], a] }
      .ifEmpty { error "Regex failed to extract chromosome ID from archaic VCF index: ${params.vcf_regex}" }
      .tap { archaic_tbis }
      .subscribe { println "Added ${it[1]} to archaic_tbis channel" }

   //Do the same for the archaic VCF channel:
   archaic_vcfs_autosomes = archaic_vcfs
      .join(archaic_tbis, by: 0, failOnDuplicate: true, failOnMismatch: true)
      .filter({ autosome_num_list.tokenize(',').contains(it[0]) })
      .flatMap({ a -> [a[1], a[2]] })
      .collect()

   //Set up the channel of recombination rate map files for Sprime:
   Channel
      .fromPath(params.recmap_glob, checkIfExists: true)
      .ifEmpty { error "Unable to find recombination rate map files matching glob: ${params.recmap_glob}" }
      .map { a -> [ (a =~ params.recmap_regex)[0][1], a] }
      .ifEmpty { error "Regex failed to extract chromosome ID from recombination rate map: ${params.recmap_glob} =~ ${params.recmap_regex}" }
      .tap { recmaps }
      .subscribe { println "Added chr${it[0]} recombination rate map (${it[1]}) to recmaps channel" }

   //Retain only the autosomes for the recombination rate maps:
   recmaps_autosomes = recmaps
      .filter({ autosome_num_list.tokenize(',').contains(it[0]) })

   //Set up the channel of target populations for Sprime:
   Channel
      .fromPath(params.target_pops_file, checkIfExists: true)
      .ifEmpty { error "Unable to find target populations file for Sprime: ${params.target_pops_file}" }
      .splitCsv(sep: "\t")
      .map { it[0] }
      .tap { sprime_pops }
      .tap { sprime_target_pops }
      .subscribe { println "Added ${it} to sprime_pops channel" }

   //Set up the file channels for the metadata files and the annotation GFF:
   metadata = file(params.metadata_file, checkIfExists: true)
   archaic_metadata = file(params.arc_metadata_file, checkIfExists: true)
   //And a file channel for the sample IDs to exclude:
   excluded_samples = file(params.samples_to_exclude, checkIfExists: true)

   //Subset the input VCFs for use with Sprime (i.e. retain target and outgroup pops and no missing genotypes):
   vcfsubset_inputs = sprime_pops.combine(perchrom_vcfs_subset)
   vcfsubset(vcfsubset_inputs, metadata, excluded_samples)
//vcfsubset.out.logs into sprime_vcf_subset_logs
//vcfsubset.out.vcf into Sprime_perchrom_vcfs,Sprime_perchrom_vcfs_project,Sprime_perchrom_vcfs_tractfreqs
//vcfsubset.out.outgroup into Sprime_outgroup,Sprime_outgroup_project
//vcfsubset.out.samples into Sprime_pop_samples

   //Prepare VCFs for concatvcf and projection steps:
   vcfs_for_concatvcf = vcfsubset.out.vcf
      .groupTuple(by: 0, size: num_autosomes)
      .tap { vcfs_for_projection }

   //Prepare single outgroup files for Sprime and projection steps:
   outgroups_for_sprime = vcfsubset.out.outgroup
      .groupTuple(by: 0, size: num_autosomes)
      .map({ a -> [ a[0], a[1].take(1) ] })
      .tap { outgroups_for_projection }

   //Concatenate the per-chromosome subsetted VCFs to serve as Sprime input:
   concatvcf(vcfs_for_concatvcf, autosome_num_list, num_autosomes)
//concatvcf.out.logs into concat_input_vcf_logs
//concatvcf.out.vcf into Sprime_vcfs

   //Run Sprime per chromosome but with concatenated VCF:
   sprime_inputs = concatvcf.out.vcf
      .join(outgroups_for_sprime, by: 0, failOnDuplicate: true, failOnMismatch: true)
      .combine(recmaps_autosomes)
   sprime(sprime_inputs, autosome_map)
//sprime.out.log into Sprime_logs
//sprime.out.score into Sprime_scores_matchrate,Sprime_scores_project,Sprime_scores_tractfreqs

   //Calculate the match rates of each tract to the high-depth archaics:
   matchrate_inputs = sprime.out.score
      .groupTuple(by: 0, size: num_autosomes)
      .tap { sprime_scores_bypop }
   matchrate(matchrate_inputs, archaic_vcfs_autosomes, archaic_metadata, autosome_num_list)
//matchrate.out.matches into Sprime_matches
//matchrate.out.matchrates into Sprime_matchrates

   //Project Sprime population-level tracts onto haplotypes:
   projection_inputs = sprime_scores_bypop
      .join(vcfs_for_projection, by: 0, failOnDuplicate: true, failOnMismatch: true)
      .join(outgroups_for_projection, by: 0, failOnDuplicate: true, failOnMismatch: true)
   projection(projection_inputs, autosome_list, has_chr_prefix)
//projection.out.projections into Sprime_project_BEDs

   //Concatenate outputs for simplicity of analysis:
   cat_outs_inputs = matchrate.out.matchrates
      .collectFile() { a -> [ "matchrate_paths.tsv", a.getSimpleName()+'\t'+a.getName()+'\t'+a+'\n' ] }
      .combine(projection.out.projections.collectFile() { a -> [ "projection_paths.tsv", a.getSimpleName()+'\t'+a.getName()+'\t'+a+'\n' ] })
      .combine(sprime_target_pops.collectFile(name: 'Sprime_target_populations.txt', newLine: true, sort: true))
   cat_outs(cat_outs_inputs, metadata, excluded_samples)
//catouts.out.outputs into catouts_outputs
}
