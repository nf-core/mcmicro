process omextractor {
	tag ""
	label "process_single"

	container "docker.io/labsyspharm/bftools:latest"

	input:
	path(image)

	output:
	path "ome.xml", emit: xml
	path "versions.yml", emit: versions

	when:
	task.ext.when == null || task.ext.when

	script:
	"""
	showinf -omexml-only -nopix -no-upgrade $image > "ome.xml"

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
        	showinf: \$(showinf -version)
    	END_VERSIONS
	"""
}
