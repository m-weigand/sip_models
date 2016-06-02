#!/usr/bin/python
import numpy as np


class sip_response(np.ndarray):

    def __init__(self, *args, **kwargs):
        super(sip_response, self).__init__(*args, **kwargs)

        self._rcomplex = None
        self._ccomplex = None

    def rmag(self):
        print('PLA')

    def rpha(self):
        pass


obj = sip_response((3)) * 0
print obj
