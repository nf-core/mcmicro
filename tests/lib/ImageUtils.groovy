@GrabResolver(name='scijava', root='https://maven.scijava.org/content/repositories/public/')
@Grab(group='ome', module='formats-bsd', version='7.2.0')

import loci.formats.ImageReader
import loci.formats.services.OMEXMLService
import loci.common.services.ServiceFactory
import ome.units.UNITS

static Map getImageMetadata(String path) {
    def factory = new ServiceFactory()
    def service = factory.getInstance(OMEXMLService)
    def om = service.createOMEXMLMetadata()
    def reader = new ImageReader()
    reader.setMetadataStore(om)
    reader.setFlattenedResolutions(false)
    reader.setId(path)
    assert om.imageCount == 1
    def metadata = [
        format: reader.format,
        type: om.getPixelsType(0).value,
        sizeX: om.getPixelsSizeX(0).numberValue,
        sizeY: om.getPixelsSizeY(0).numberValue,
        sizeZ: om.getPixelsSizeZ(0).numberValue,
        sizeC: om.getPixelsSizeC(0).numberValue,
        sizeT: om.getPixelsSizeT(0).numberValue,
        physicalSizeX: om.getPixelsPhysicalSizeX(0)?.value(UNITS.MICROMETER),
        physicalSizeY: om.getPixelsPhysicalSizeY(0)?.value(UNITS.MICROMETER),
        physicalSizeZ: om.getPixelsPhysicalSizeZ(0)?.value(UNITS.MICROMETER),
    ]
    return metadata
}
