"""Software specific Exceptions

Used for exceptions that can occur when using seqQscorer.

date:	2019-05-12
author:	Steffen Albrecht

"""

class WrongFeatureInputException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class WrongSettingException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class WrongOutputFileException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class IncorrectModelException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

