/*

Takes a list of paths and outputs a single path with all inputs symlinked into
it.  Used to stage https:// paths to a local directory for sample-based workflow
tests, which require a listable directory.

*/
process DIR_COMBINE {
    input:
    tuple val(meta), path(f)

    output:
    tuple val(meta), path("combined")

    publishDir enabled: false

    script:
    def relPaths = f.collect{ "'../$it'" }.join(' ')
    """
    mkdir combined
    cd combined
    ln -s $relPaths .
    """
}
