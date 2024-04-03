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

// Parse and round floating point values to 13 decimal places. Pass other
// strings through untouched.
static roundIfDouble(String value) {
    value.contains(".")
        ? Math.round(Double.parseDouble(value) * 1e13) / 1e13
        : value
}
