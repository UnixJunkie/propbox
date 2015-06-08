from __future__ import absolute_import

import os
import sys

from rdkit import Chem

# A simple subset of what chemfp's rdkit_toolkit can do.

def _open(source, compression):
    if compression == "gz":
        import gzip
        if isinstance(source, basestring):
            infile = gzip.Gzip(filename)
            close = infile.close
        elif source is None:
            infile = gzip.GzipFile(sys.stdin)
            close = None
        else:
            infile = gzip.GzipFile(fileobj=source)
            close = None
            
    elif compression:
        raise ValueError("Unsupported compression %r" % (compression,))
    
    else:
        if isinstance(source, basestring):
            infile = open(source, "rb")
            close = infile.close
        elif source is None:
            infile = sys.stdin
            close = None
        else:
            infile = source
            close = None

    return infile, close
        


class Format(object):
    toolkit_name = "rdkit"
    is_input_format = True
    is_output_format = True
    is_available = True
    support_io = True
    def __init__(self, name, compression):
        self.name = name
        self.compression = compression
    def __repr__(self):
        if self.compression:
            name = "rdkit/%s.%s" % (self.name, self.compression)
        else:
            name = "rdkit/%s" % (self.name,)
        return "Format(%r)" % (name,)

_delimiter_settings = {
    "tab": "tab",
    "whitespace": "whitespace",
    "space": "space",
    "to-eol": "to-eol",
    "\t": "tab",
    " ": "space",
    }
    
def _get_delimiter(d, format, text_settings):
    v = text_settings.get("delimiter", None)
    if v is None:
        return
    try:
        delimiter = _delimiter_settings[v]
    except KeyError:
        raise ValueError("%r does not support the delimiter %r"
                         % (format, v))
    d["delimiter"] = v

def _get_boolean(d, format, text_settings, name):
    v = text_settings.get(name, None)
    if v is None:
        return
    if v in ("1", "True", "true"):
        value = True
    elif v in ("0", "False", "false"):
        value = False
    else:
        raise ValueError("%r setting %r is not a supported boolean: %r"
                         % (format, name, v))
    d[name] = value
    
class SmilesFormat(Format):
    def __init__(self, compression):
        super(SmilesFormat, self).__init__("smi", compression)
        
    def get_reader_args_from_text_settings(self, text_settings):
        d = {}
        _get_delimiter(d, self, text_settings)
        _get_boolean(d, self, text_settings, "has_header")
        _get_boolean(d, self, text_settings, "sanitize")
        return d

class SDFFormat(Format):
    def __init__(self, compression):
        super(SDFFormat, self).__init__("sdf", compression)
    
    def get_reader_args_from_text_settings(self, text_settings):
        d = {}
        _get_boolean(d, self, text_settings, "strictParsing")
        _get_boolean(d, self, text_settings, "removeHs")
        _get_boolean(d, self, text_settings, "sanitize")
        return d

_supported_formats = {
    "smi": (SmilesFormat, ""),
    "smi.gz": (SmilesFormat, "gz"),
    "sdf": (SDFFormat, ""),
    "sdf.gz": (SDFFormat, "gz"),
    }
    
def get_input_format_from_source(source=None, format=None):
    if format is not None:
        if not isinstance(format, basestring):
            return format
        try:
            klass, arg = _supported_formats[format]
            return klass(arg)
        except KeyError:
            raise ValueError("Unsupported format %r" % (format,))

    if source is None:
        return SmilesFormat("")
    elif isinstance(source, basestring):
        filename = source
    else:
        filename = getattr(source, "name", "<unknown>")

    filename, ext = os.path.splitext(filename)
    ext = ext.lower()
    if ext == ".gz":
        compression = "gz"
        filename, ext = os.path.splitext(filename)
    else:
        compression = ""

    if ext in ("smi", "ism", "isosmi"):
        return SmilesFormat(compression)
    elif ext in ("sdf", "sd", "mdl"):
        return SDFFormat(compression)
    else:
        # Don't know. Assume it's in SMILES format.
        return SmilesFormat(compression)


####

_delimiter_method = {
    " ": " ",
    "\t": "\t",
    "space": " ",
    "tab": "\t",
    "whitespace": "whitespace",
    "to-eol": "to-eol",
    None: "to-eol",
    }

def _read_smiles(infile, id_tag, reader_args):
    delimiter = reader_args.get("delimiter", None)
    delimiter = _delimiter_method[delimiter]
    
    titleLine = bool(reader_args.get("titleLine", False))
    sanitize = bool(reader_args.get("sanitize", True))

    lineno = 1
    if titleLine:
        next(infile)  # Will raise an exception if the file is empty
        lineno += 1

    if delimiter == "to-eol":
        for lineno, line in enumerate(infile, lineno):
            line = line.rstrip("\n\r")
            terms = line.split(None, 1)
            n = len(terms)
            if n == 0:
                raise ValueError("Line %d is empty" % (lineno,))
            elif n == 1:
                smiles = terms[0]
                id = None
            else:
                smiles = terms[0]
                id = terms[1]

            mol = Chem.MolFromSmiles(smiles, sanitize)
            if mol is not None:
                mol.SetProp("_input_smiles", smiles)
            yield id, mol

    elif delimiter == "whitespace":
        for lineno, line in enumerate(infile, lineno):
            terms = line.split()
            n = len(terms)
            if n == 0:
                raise ValueError("Line %d is empty" % (lineno,))
            elif n == 1:
                smiles = terms[0]
                id = None
            else:
                smiles = terms[0]
                id = terms[1]

            mol = Chem.MolFromSmiles(smiles, sanitize)
            if mol is not None:
                mol.SetProp("_input_smiles", smiles)
            yield id, mol

    else:
        for lineno, line in enumerate(infile, lineno):
            line = line.rstrip("\r\n")
            terms = line.split(delimiter)
            n = len(terms)
            if n == 0:
                raise ValueError("Line %d is empty" % (lineno,))
            elif n == 1:
                smiles = terms[0]
                id = None
            else:
                smiles = terms[0]
                id = terms[1]

            mol = Chem.MolFromSmiles(smiles, sanitize)
            if mol is not None:
                mol.SetProp("_input_smiles", smiles)
            yield id, mol
    

def _read_sdf(infile, id_tag, reader_args):
    sanitize = bool(reader_args.get("sanitize", True))
    removeHs = bool(reader_args.get("removeHs", True))
    strictParsing = bool(reader_args.get("strictParsing", True))

    reader = Chem.ForwardSDMolSupplier(infile)
    for mol in reader:
        if mol is None:
            id = None
        else:
            if id_tag is None:
                id = mol.GetProp("_Name")
                if not id:
                    id = None
            else:
                try:
                    id = mol.GetProp(id_tag)
                except KeyError:
                    id = None
        yield id, mol
                
def read_ids_and_molecules(source=None, format=None, id_tag=None, reader_args=None):
    format = get_input_format_from_source(format)
        
    if format.name in ("smi", "ism", "can", "usm"):
        f = _read_smiles
    elif format.name in ("sdf",):
        f = _read_sdf
    else:
        raise ValueError("Unsupported RDKit formt %r" % (format.name,))
    
    if reader_args is None:
        reader_args = {}

    infile, close = _open(source, format.compression)
        
    if close is None:
        return f(infile, id_tag, reader_args)
    else:
        return _with_close(f, infile, close, id_tag, reader_args) 

def _with_close(f, infile, close, id_tag, reader_args):
    try:
        for x in f(infile, id_tag, reader_args):
            yield x
    finally:
        close()
