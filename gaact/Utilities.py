__author__ = 'Donghui'

import os
import json
import Exceptions

VALID_ALLELES = ["A", "T", "C", "G"]


def contentsWithRegexpFromFolder(folder, regexp):
    contents = os.listdir(folder)
    paths = [x for x in contents if regexp.match(x)] if regexp else contents
    return paths


import logging


import gzip
class FileIterator(object):
    def __init__(self, path, header=None, compressed = False, ignore_until_header = False):
        self.path = path
        self.compressed = compressed
        self.header = header
        self.ignore_until_header = ignore_until_header
        if ignore_until_header and not header:
            raise Exceptions.InvalidArguments("File iterator received conflicting header information")

    def iterate(self,callback=None):
        if self.compressed:
            with gzip.open(self.path, 'rb') as file_object:
                self._iterateOverFile(file_object, callback)
        else:
            with open(self.path, 'rb') as file_object:
                self._iterateOverFile(file_object, callback)

    def _iterateOverFile(self, file_object, callback):
        if self.ignore_until_header:
            self._ignore_until_header(file_object)
        else:
            if self.header is not None:
                line = file_object.readline().strip("\n")
                if len(self.header) and line != self.header:
                    raise Exceptions.MalformedInputFile(self.path, "Unexpected header")

        self._processFile(file_object, callback)

    def _ignore_until_header(self, file_object):
        if self.ignore_until_header and self.header:
            skip = True
            while skip:
                l = file_object.readline()
                if not l:
                    raise Exceptions.InvalidArguments("Wrong header lookup in %s" % (self.path,))
                l = l.strip()
                if self.header in l:
                    skip = False

    def _processFile(self, file_object, callback):
        if callback is not None:
            for i,line in enumerate(file_object):
                callback(i, line)


def ensure_requisite_folders(path):
    folder = os.path.split(path)[0]
    if len(folder) and not os.path.exists(folder):
        os.makedirs(folder)
