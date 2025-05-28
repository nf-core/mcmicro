process omextractor {
	tag ""
	label "process_single"

	container "docker.io/labsyspharm/bftools:latest"

	input:
	val(meta)
	path(image)

	output:
        tuple val(meta), path("ome.xml") emit: xml
	path "versions.yml", emit: versions

	when:
	task.ext.when == null || task.ext.when

	script:
	"""
	showinf -omexml-only -nopix -no-upgrade -option zeissczi.autostitch false -option zeissczi.attachments false $image > "ome.xml"

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
        	showinf: \$(showinf -version)
    	END_VERSIONS
	"""
}
