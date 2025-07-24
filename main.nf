process getCellMeta {
	tag "$study_name"
	label 'getCellMeta'

	conda "/home/rschwartz//anaconda3/envs/scanpyenv"
	publishDir "${params.outdir}/cell_meta"

	input:
	val(study_name)

	output:
	tuple val(study_name), path("**cell_level_meta.tsv.gz"), emit: cell_level_meta

	script:
	"""
	python $projectDir/bin/get_gemma_cellmeta.py --study_name $study_name
	"""
}

process getFullMatrix {
	tag "$study_name"
	label 'getFullMatrix'

	publishDir "${params.outdir}/full_matrices/${study_name}"

	conda "/home/rschwartz//anaconda3/envs/scanpyenv"

	input:
	val(study_name)

	output:
	tuple val(study_name), path("**full_matrix.tsv.gz"), emit: full_matrix

	script:
	"""
	gemma-cli-staging getSingleCellDataMatrix -e ${study_name} \\
						--format CELL_BROWSER \\
						--use-bioassay-ids \\
						--output-file ${study_name}_full_matrix.tsv.gz 
	"""
}

process runCbScanpy {
	tag "$study_name"
	label 'runCbScanpy'

	conda "/home/rschwartz//anaconda3/envs/scanpyenv"

	//publishDir "${params.cb_build_dir}/"
	publishDir "${params.outdir}"

	input:
	tuple val(study_name), path(cell_level_meta), path(full_matrix)

	output:
	tuple val(study_name), path("${new_study_name}/"), emit: cb_dir_channel

	script:
	new_study_name = study_name.replaceAll("\\.", "_")
	"""
	cbScanpy -e ${full_matrix} \\
		  -m ${cell_level_meta} \\
		  -o ${new_study_name} \\
		  -n ${new_study_name} \\
		  -c ${params.cbscanpy_config}
	"""
	//	  -c ${params.cbscanpy_config}
}

process runCbBuild {
	tag "$study_name"
	label 'cbBuild'
	conda "/home/rschwartz//anaconda3/envs/scanpyenv"

	publishDir "${params.cb_outdir}"

	input:
	tuple val(study_name), path(cb_dir)

	output:
	path "${new_study_name}/"

	script:
	new_study_name = study_name.replaceAll("\\.", "_")
	""" 
	cbBuild -i ${cb_dir}/cellbrowser.conf \\
		-o ${new_study_name} 
	"""

}

workflow {
	// Define the study name
	Channel.from(params.study_names)
	.set { study_names }	

	// Call the process to get cell metadata
	getCellMeta(study_names)

	getFullMatrix(study_names)

	getCellMeta.out.cell_level_meta
	.set { cell_level_meta }

	getFullMatrix.out.full_matrix
	.set { full_matrices }

	cell_level_meta.combine(full_matrices, by: 0)
	.set { cbInput }

	// Run the cbScanpy process
	runCbScanpy(cbInput)

	runCbScanpy.out.cb_dir_channel
	.set { cb_dir_channel }

	runCbBuild(cb_dir_channel)
}