process getCellMeta {
	tag "$study_name"
	label 'getCellMeta'

	input:
	val(study_name)

	output:
	tuple val(study_name), path("**cell_level_meta.tsv.gz"), emit: cell_level_meta

	script:
	"""
	gemma-cli getSingleCellMetadata -e ${study_name} --use-raw-column-names --use-bioassay-ids \\
	--output-file ${study_name}_cell_level_meta.tsv.gz
	"""
}

process cleanCellMeta {
	tag "$study_name"

	input:
	tuple val(study_name), path(cell_level_meta)
	

	output:
	tuple val(study_name), path("**cleaned_cell_level_meta.tsv.gz"), emit: cleaned_cell_level_meta

	script:
	def remove_cols_pattern = params.remove_cols.join('|')
	"""
	gunzip -c ${cell_level_meta} > ${study_name}_cell_level_meta.tsv
	REMOVE_COLS="${remove_cols_pattern}"
    HEADER=\$(head -1 ${study_name}_cell_level_meta.tsv)
    KEEP_COLS=\$(awk -F'\\t' -v pat="\$REMOVE_COLS" '{for(i=1;i<=NF;i++){if(\$i !~ pat){printf "%d,",i}}}' <<< "\$HEADER" | sed 's/,\$//')
    cut -f\$KEEP_COLS ${study_name}_cell_level_meta.tsv > ${study_name}_cleaned_cell_level_meta.tsv
	gzip ${study_name}_cleaned_cell_level_meta.tsv
	"""

}

process getFullMatrix {
	tag "$study_name"
	label 'getFullMatrix'
	conda "/home/rschwartz//anaconda3/envs/scanpyenv"

	input:
	val(study_name)

	output:
	tuple val(study_name), path("**full_matrix.tsv.gz"), emit: full_matrix

	script:
	def gemma_cmd = params.use_staging ?  "gemma-cli-staging" : "gemma-cli" 
	"""
	${gemma_cmd} getSingleCellDataMatrix -e ${study_name} \\
						--format CELL_BROWSER \\
						--use-bioassay-ids \\
						--no-cursor-fetch \\
						--no-auto-flush \\
						--use-raw-column-names \\
						--output-file ${study_name}_full_matrix.tsv.gz 
	"""
}

process runCbScanpy {
	tag "$study_name"
	label 'runCbScanpy'

	conda "/home/rschwartz//anaconda3/envs/scanpyenv"

	//publishDir "${params.cb_build_dir}", mode: 'copy'

	input:
	tuple val(study_name), path(cell_level_meta), path(full_matrix)

	output:
	tuple val(study_name), path("${new_study_name}/"), emit: cb_dir_channel

	script:
	// new_study_name = study_name.replaceAll("\\.", "_")
	new_study_name = study_name.replaceAll("[^A-Za-z0-9_]", "_")
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

	//publishDir "${params.cb_outdir}", mode: 'copy'

	input:
	tuple val(study_name), path(cb_dir)

	output:
	path "${new_study_name}/"

	script:
	//new_study_name = study_name.replaceAll("\\.", "_")
	new_study_name = study_name.replaceAll("[^A-Za-z0-9_]", "_")
	""" 
	cbBuild -i ${cb_dir}/cellbrowser.conf \\
		-o ${params.cb_outdir} \\
	"""

}



process collect_cb_dirs {
	tag "$study_name"
	label 'collectCbDirs'

	// "${params.cb_build_dir}", mode: 'copy'

	input:
	val cb_dir_channel

	output:
	path "parent_dir/", emit: cb_parent_dir
	

	script:
	"""
	mkdir -p parent_dir
	for dir in ${cb_dir_channel}; do
		echo "Moving $dir to parent_dir"
		mv $dir parent_dir/
	done
	
	"""
}

workflow {
	// Define the study name
//Channel
   // .from(params.study_names.split(/\s+/))
   // .set { study_names }


 	study_names = (
            file(params.study_names).exists() && file(params.study_names).isFile()
                ? Channel.from(
                      file(params.study_names)
                          .readLines()
                          .collect { it.trim() }
                          .findAll { it }
                          .collect { it }
                  )
                : Channel.from(
                      params.study_names
                          .split(/\s+/)
                  )
        )

	// Call the process to get cell metadata
	getCellMeta(study_names)
	.set { cell_meta }

	cleanCellMeta(cell_meta)
	.set {clean_cell_meta}

	getFullMatrix(study_names)

	getFullMatrix.out.full_matrix
	.set { full_matrices }

	clean_cell_meta.combine(full_matrices, by: 0)
	.set { cbInput }

	// Run the cbScanpy process
	runCbScanpy(cbInput)

	runCbScanpy.out.cb_dir_channel
	.set { cb_dir_channel }

	runCbBuild(cb_dir_channel)
}

workflow.onComplete {
	// print where output files are located
	println "Cell Browser build complete. Output files are located in: ${params.cb_outdir}"
	// print workdir
	println "Work directory: ${workflow.workDir}"

	println "View the dataset at https://dev.gemma.msl.ubc.ca/cellbrowser/"
}

workflow.onError = {
println "Error: something went wrong, check the pipeline log at '.nextflow.log"
}