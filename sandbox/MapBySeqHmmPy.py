def list_x_attributes(future_class_name, future_class_parents, future_class_attr):
  """
    Return a class object, with the list of its attribute turned
    into uppercase.
  """

  # pick up any attribute that doesn't start with '__' and uppercase it
  x_attr = []
  for name, val in future_class_attr.items():
      if name.startswith('x_') :
          x_attr.append(name)

  future_class_attr['_x_attr_'] = x_attr
  # let `type` do the class creation
  return type(future_class_name, future_class_parents, future_class_attr)


class Analyzer():

    def __init__(read_instance):
        assert(type(read_instance) == rawData)
    
        self.data = read_instance.provide_data()
        
#    def plot():

# class SpecificAnalyzer(Analyzer) :
#     pass

#######################################################

class filterMask:
    MAX_INT = 2^63-1
    
    def __init__(self, fieldNames):
		
        self.intFilterDict = {}
        self.intFilterDict['altCounts'] = [0,  self.MAX_INT]
        self.intFilterDict['refCounts'] = [0,  self.MAX_INT]
        self.intFilterDict['DeltaX'] = [0,  self.MAX_INT] # distance to the closest neighbours
        self.intFilterDict['totCounts'] = [0,  self.MAX_INT]
        
        self.allelesDict = {}
        self.allelesDict['refAllele'] = ['A', 'T']
        self.allelesDict['altAllele'] = ['G', 'C']
        
    
    def getSql(self):
        assert( len(self.intFilterDict)>0 )
        self.sqlString = [];
        for eName, eValue in self.intFilterDict.items():
            if eValue[0] > 0:
                self.sqlString.append('%s > %u' %(eName, eValue[0]) )
            if eValue[1] < self.MAX_INT:
                self.sqlString.append('%s < %u' %(eName, eValue[1]) )
        return ' AND '.join(self.sqlString)
#######################################################

class chromosomeData:    
    def getData(self, chrName):
        """retrieve data for each position on a given chromosome from a database"""
        curs = self.chromosomes[chrName].conn.cursor()
        for ff in chrDataFields.items():
            curs.execute('SELECT %s FROM stocks WHERE %s' % (ff, filterMask) )
            chrDataFields[ff] = curs.fetchall()
        
        
class rawData(chromosomeData):
    def __init__(self, db, chrNames, ):
        assert(type(db) == snpDbSqlite)
        for chrN in chrNames:
            self.chromosomes[chrN] = db.getSqlConnection(chrN)
		
        self.chrDataFields = {}
        for ff in fieldNames:
            self.chrDataFields[ff] = []

    
class snpDbSqlite:
    def __init__(self, dbPath):
    
    
############

        self.chromosome
        self.position
        self.refCounts
        self.altCounts
        
        self.properties = {};
