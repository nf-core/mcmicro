@Grab(group='com.opencsv', module='opencsv', version='5.9')
import com.opencsv.CSVReader
import com.opencsv.CSVWriter

import com.askimed.nf.test.util.ObjectUtil
import java.nio.file.Paths

static String roundAndHashCsv(String path, int precision = 12) {
    def name = Paths.get(path).getFileName().toString()
    def csvContent = roundCsv(path, precision)
    return name + ":rounded:md5," + ObjectUtil.getMd5(csvContent)
}

static String roundCsv(String path, int precision) {
    def reader = new CSVReader(new File(path).newReader())
    def strWriter = new StringWriter()
    def writer = new CSVWriter(strWriter)
    // Copy header row.
    writer.writeNext(reader.readNext(), false)
    def row
    while (row = reader.readNext()) {
        row = row.collect{ roundIfDouble(it, precision) }
        writer.writeNext(row as String[], false)
    }
    def csvContent = strWriter.toString()
    return csvContent
}

static Map summarizeCsv(String path) {
    def reader = new CSVReader(new File(path).newReader())
    def strWriter = new StringWriter()
    def headers = reader.readNext() as List
    def count = 0
    while (reader.readNext()) {
        count++
    }
    return [
        headers: headers,
        rowCount: count,
    ]
}

// Parse and round floating point values to the specified number of decimal
// digits of precision. Pass other strings through untouched.
static String roundIfDouble(String value, int precision) {
    assert precision > 0
    return value.contains(".")
    ? sprintf("%.${precision}g", Double.parseDouble(value))
    : value
}
