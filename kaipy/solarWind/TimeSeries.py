
class TimeSeries(dict):
    """
    TimeSeries is a standard Python Dictionary with some helper
    functions useful for data processing collections of time series.
    """
    def __init__(self, indict={}):
        dict.__init__(self)
        # not "self.__keydict" because I want it to be easily accessible by subclasses!
        self._keydict = {}

        for entry in indict:
            self[entry] = indict[entry]
                                    
    def append(self, key, name, units, data):
        """
        key: unique identifier of variable
        name: more descriptive name of variable
        data: time series data array
        units: units (grams, for example)
        """
        #if self.has_key(key):
        if key in self:
            raise KeyError('Error:  Dictionary already has key \"%s\".' % key)
        self[key] = { 'name': name, 'data': data, 'units': units }
        
    def getUnits(self, key):
        """
        >>> ts = TimeSeries({'v': {'name': 'Velocity', 'data': [1,2,3], 'units': 'm/s'}})
        >>> ts.getUnits('v')
        'm/s'
        """
        return self[key]['units']

    def setUnits(self, key, units):
        self[key]['units'] = units

    def getName(self, key):
        """
        >>> ts = TimeSeries({'v': {'name': 'Velocity', 'data': [1,2,3], 'units': 'm/s'}})
        >>> ts.getName('v')
        'Velocity'
        """
        return self[key]['name']

    def setName(self, key, name):
        self[key]['name'] = name

    def getData(self, key):
        """
        >>> ts = TimeSeries({'v': {'name': 'Velocity', 'data': [1,2,3], 'units': 'm/s'}})
        >>> ts.getData('v')
        [1, 2, 3]
        """
        return self[key]['data']

    def setData(self, key, data, index=None):
        try:
            self[key]['data'][index] = data
        except TypeError:
            self[key]['data'] = data

if __name__ == '__main__':
    import doctest
    doctest.testmod()
