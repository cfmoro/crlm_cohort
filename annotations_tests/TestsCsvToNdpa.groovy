import com.opencsv.CSVParser
import com.opencsv.CSVParserBuilder
import com.opencsv.CSVReaderBuilder
import com.opencsv.CSVReader
import com.opencsv.CSVWriter
import groovy.io.FileType
import groovy.xml.MarkupBuilder

BASE_DIR = '/home/bibu/Workspace/crlm_cohort/annotations/'
TESTS_CSV_FILE = BASE_DIR + 'Annotation_tests_CRLM_cohort.csv'

def colors = ['D' : '#00ff00', 'R' : '#ff0000', 'R2' : '#ff00ff', 'P' : '#0000ff', '%' : '#ffff00']

// Open csv file for read
CSVParser parser = new CSVParserBuilder().withSeparator(',' as char).build()
CSVReader reader = new CSVReaderBuilder(new File(TESTS_CSV_FILE).newReader()).withCSVParser(parser).build()
	
// Iterate and retrieve into a list rows in csv file
def annotationRows = []
def rowNum = 0
while ((aRow = reader.readNext()) != null) {	
   rowNum += 1
   println "Parsing row number ${rowNum}: ${aRow}"
   
   // Skip header
   if(rowNum == 1) {	   
	   continue
   }
   
  def rowWithId = [rowNum]
  rowWithId.addAll(aRow)
  
  annotationRows << rowWithId      
}
println "Read rows: \n ${annotationRows}"

// Group rows by file name
def annotationsByFiles = annotationRows.groupBy { it -> it[1] }
println "Rows mapped by file name: \n" + annotationsByFiles

// Iterate over file names and corresponding rows and create the ndpa files
annotationsByFiles.each { key, value ->

	// Write ndpa compatible xml file
	def fileName = BASE_DIR + key + ".ndpi.ndpa"
	def writer = new FileWriter(fileName)
	def xml = new MarkupBuilder(writer)
	def yPos = -8984090
	xml.annotations() {
		
			value.each { annotRow ->
				
				def (annotId, annotLabel, annotValue) = [annotRow[0], annotRow[2], annotRow[3] as Integer] 
				def annotTitle = annotLabel == '%' ? annotValue + '%' : annotLabel
				def annotationDisplName = annotLabel == '%' ? 'AnnotatePointer' : 'AnnotateFreehandLine' 			
				def annotationType = annotLabel == '%' ? 'pointer' : 'freehand'
					 						
				ndpviewstate(id: annotId) {
					title(annotTitle)
					details('')
					coordformat('nanometers')
					lens('0.393133')
					x('0')
					y('0')
					z('0')
					showtitle('0')
					showhistogram('0')
					showlineprofile('0')
					annotation(type: annotationType, displayname: annotationDisplName, color : colors[annotLabel]) {  
						if(annotLabel != '%') {
							measuretype('1')
							closed('0')
							pointlist() {
								point() {
									x('0')
									y(yPos)
								}
								if(annotValue > 20) { // Annotations greater than 20 um draw correctly in NdpViewer
									point() {
										x('10000')
										y(yPos)
									}
									point() {
										x('20000')
										y(yPos)
									}
								} else {	// Small annotations < 20 um, don't draw correctly but can be parsed
									point() {
										x('10')
										y(yPos)
									}
									point() {
										x('20')
										y(yPos)
									}
								}								
								point() {
									x((annotValue) * 1000)
									y(yPos)
								}
							}
						} else {
							x1('15350931')
							y1('-2455415')
							x2('19850348')
							y2('-6299054')
						}	
					}
				}
			yPos += 500000
			}
	}		
		// Test some values in xml
		//def annotations = new XmlSlurper().parseText(writer.toString())
		//println annotations.ndpviewstate.first().@id
		//assert annotations.ndpviewstate.first().@id == '2'
		
	writer.close()
	println "Created annotations file: " + fileName
}
println 'Done'