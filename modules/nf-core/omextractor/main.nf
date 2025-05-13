process omextractor {
	tag ""
	label "process_single"

	container "labsyspharm/bftools"

	input:
	path(image)

	output:
	path "ome.xml", emit: xml

	when:
	task.ext.when == null || task.ext.when

	script:
	"""
	tiffcomment $image > "ome.xml"
	"""
}
