#!/usr/bin/env nextflow

process insurveyor_calling {
	cpus 8
	maxForks 4
	module 'INSurVeyor'
	memory params.ins_mem
	publishDir "insurveyor_vcf"
	input:
		tuple val(core), path(f)
	output:
		tuple val(core), path(ins_out), path("${ins_out}.tbi")
	script:
	ins_out="${core}.insurveyor.vcf.gz"
	"""
	set -euxo pipefail
	insurveyor.py --threads 8 ${f[0]} ./ ${params.genome}
	mv out.pass.vcf.gz ${ins_out}
	tabix -p vcf ${ins_out}
	"""	
}


process svaba_calling {
	errorStrategy 'finish'
	cpus 24
	maxForks 4
	module 'svaba'
	memory params.sva_mem
	publishDir "svaba_vcf"
	input:
		tuple val(core), path(f)
	output:
		tuple val(core), path(sva_conc_sort)
	script:
	sva_sv="${core}.svaba.sv.vcf.gz"
	sva_ind="${core}.svaba.indel.vcf.gz"
	sva_conc_sort="${core}.svaba.vcf.gz"
	"""
	set -euxo pipefail  
	/usr/bin/time -o "resources.t" -f "%e %M" svaba run -p 24 -G ${params.genome} -t ${f[0]} -z 
	mv no_id.svaba.sv.vcf.gz ${sva_sv}
	mv no_id.svaba.indel.vcf.gz ${sva_ind}
	bcftools index -t ${sva_sv}
	bcftools index -t ${sva_ind}   
	bcftools concat ${sva_sv} ${sva_ind} -a -D -Oz -o - | bcftools sort -Oz -o ${sva_conc_sort}
	"""	
}


process fix_files{
	publishDir "svaba_fixed_vcf"
	input:
		tuple val (core), path(fs)
	output:
		tuple val(core), path(sva_fix), path("${sva_fix}.tbi")  
	script:
	sva_fix="${core}.svaba.fixed.vcf.gz"
	"""
	set -euxo pipefail
	zcat "${core}.svaba.vcf.gz" | sd "ID=GQ,Number=1,Type=String" "ID=GQ,Number=1,Type=Integer" | sd 'ID=PL,Number=[^,]*' 'ID=PL,Number=G' | bgzip > ${sva_fix}
	tabix -p vcf ${sva_fix}	
	"""
}


process survclusterer_ensemble {
	cpus 4
	maxForks 4
	memory params.clu_mem
	publishDir "clusterer_vcf"
	input:
		tuple val(core), path(files)
	output:
		tuple path(clu_out), path("${clu_out}.tbi")
	script:
	clu_out = "${core}.clusterer.vcf.sv.gz"
	"""
	#!/bin/bash
	set -euxo pipefail
	
	if [ -f "file.txt" ]; then
		rm "file.txt"
	fi
	
	for file in ${files}; do
		if [[ \$file == *.vcf.gz ]]; then
			echo \$file >> file.txt
		fi
	done
	
	clusterer -d 1000 -t 4 file.txt ${params.genome} -o ${core}.clusterer.vcf
	bgzip -k ${core}.clusterer.vcf.sv
	tabix ${clu_out}
	"""	
}


process survivor_ensemble {
	cpus 2
	maxForks 4
	memory params.sur_mem
	publishDir "survivor_vcf"
	input:
		tuple val(core), path(files)
	output:
		tuple path(sur_out), path("${sur_out}.tbi")
	script:
	sur_out = "${core}.surv_rem.vcf.gz"
	"""
	#!/bin/bash
	set -euxo pipefail

	if [ -f "file.txt" ]; then
		rm "file.txt"
	fi
	
	for file in ${files}; do
		if [[ \$file == *.vcf.gz ]]; then
			zcat \$file > "\$file.vcf"
			echo "\$file.vcf" >> file.txt
		fi
	done
	
	SURVIVOR merge file.txt 1000 1 1 -1 -1 -1 ${core}.survivor.vcf
	bgzip ${core}.survivor.vcf
	zgrep -v HLA ${core}.survivor.vcf.gz | bgzip |  bcftools sort  -Oz -o ${core}.survivor_sorted.vcf.gz
	bcftools view -e 'INFO/END < POS' ${core}.survivor_sorted.vcf.gz -o ${sur_out}
	bcftools index -t ${sur_out}
	"""	
}


process truvari_ensemble {
	cpus 2
	maxForks 4
	module 'truvari'
	memory params.tru_mem
	publishDir "truvari_vcf"
	input:
		path(files)
	output:
		tuple path(tru_out), path("${tru_out}.tbi")
	script:
	tru_out = "tru_coll_sorted.vcf.gz"
	"""
	#!/bin/bash
	set -euxo pipefail
	vcf_files="" 
	for file in ${files}; do
		if [[ \$file == *.vcf.gz ]]; then
			vcf_files="\$vcf_files \$file"
		fi
	done
	bcftools merge -m none \$vcf_files -Oz -o bcftools.merge.vcf.gz --force-samples
	bcftools index -t bcftools.merge.vcf.gz
	truvari collapse -i bcftools.merge.vcf.gz -o truvari.vcf.gz -c tru_collapse.vcf.gz
	vcf-sort tru_collapse.vcf.gz | bgzip -c > tru_coll_sorted.vcf.gz
	tabix -p vcf ${tru_out} 
	"""	
 }


workflow {
	
	input = Channel.fromFilePairs(params.input+params.samples) { f -> f.getSimpleName()}
	prev_vcf =  Channel.fromPath(params.directories.collect{it+params.ensamples}).map { f -> [f.simpleName, f]}
	
		main:
		//input.view()
		//prev_vcf.view()
		
		ins_res = insurveyor_calling(input)
		sva_res = svaba_calling(input)
		|  fix_files
		all_res = sva_res.mix(prev_vcf).mix(ins_res)
		|  groupTuple
		//all_res.view()
		
		survclusterer_ensemble(all_res)
		survivor_ensemble(all_res)
		| collect
		| truvari_ensemble
		
}
