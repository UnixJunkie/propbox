from __future__ import absolute_import, print_function


DEBUG = False

# Interface to the NCI services

import urllib2
from xml.etree import ElementTree
import sys

from . import calculate, collect_resolvers

class NCIName(object):
    __slots__ = ("classification", "name")
    def __init__(self, classification, name):
        self.classification = classification
        self.name = name
    def __repr__(self):
        return "NCIName(%r, %r)" % (self.classification, self.name)

@calculate()
def calc_nci_names(smiles):
    quoted_smiles = urllib2.quote(smiles)
    url = "http://cactus.nci.nih.gov/chemical/structure/%s/names/xml" % (quoted_smiles,)
    if DEBUG:
        sys.stderr.write("propbox.nci: urlopen(%r)\n" % (url,))
    f = urllib2.urlopen(url)
                        
    if DEBUG:
        content = f.read()
        sys.stderr.write("propbox.nci: RESPONSE FROM %r\n" % (url,))
        sys.stderr.write(content)
        sys.stderr.write("="*70 + "\n")
        etree = ElementTree.XML(content)
    else:
        etree = ElementTree.parse(f)
    
    names = []
    for item in etree.findall("data/item"):
        classification = item.attrib.get("classification", None)
        if classification is None:
            continue
        text = item.text

        names.append(NCIName(classification, text))


    if DEBUG:
        sys.stderr.write("propbox.nci: %r returns %r\n" % (url, names))
    
    return names

@calculate()
def calc_nci_iupac_name(nci_names):
    for name in nci_names:
        if name.classification == "pubchem_iupac_name":
            return name.name
    raise ValueError("No pubchem_iupac_name found")


resolver = collect_resolvers()

