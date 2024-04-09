@Grab(group='com.opencsv', module='opencsv', version='5.9')
import com.opencsv.CSVReader
import com.opencsv.CSVWriter

import com.askimed.nf.test.util.ObjectUtil
import java.nio.file.Paths

static String roundAndHashCsv(String path) {
    def name = Paths.get(path).getFileName().toString()
    def reader = new CSVReader(new File(path).newReader())
    def strWriter = new StringWriter()
    def writer = new CSVWriter(strWriter)
    def row
    while (row = reader.readNext()) {
        row = row.collect(this.&roundIfDouble)
        writer.writeNext(row as String[], false)
    }
    def csvContent = strWriter.toString()
    return name + ":rounded:md5," + ObjectUtil.getMd5(csvContent)
}

// Parse and round floating point values to 12 decimal places. Pass other
// strings through untouched.
static String roundIfDouble(String value) {
    value.contains(".")
    ? (Math.round(Double.parseDouble(value) * 1e12) / 1e12).toString()
    : value
}
