@Grab(group='com.opencsv', module='opencsv', version='5.9')
import com.opencsv.CSVReader

static List roundColumns(String path) {
    def reader = new CSVReader(new File(path).newReader())
    def headers = reader.readNext()
    def rows = [headers]
    def row
    while (row = reader.readNext()) {
        rows << row.collect(this.&roundIfDouble)
    }
    return rows
}

// Parse and round floating point values to 12 decimal places. Pass other
// strings through untouched.
static String roundIfDouble(String value) {
    value.contains(".")
    ? (Math.round(Double.parseDouble(value) * 1e12) / 1e12).toString()
    : value
}
