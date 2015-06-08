from __future__ import print_function, division, absolute_import

import sys
import argparse
import csv
import re

from . import rdkit_descriptors
from . import rdkit_toolkit as T


dialects = ["tab", "space", "whitepace"] + csv.list_dialects() 


def parse_reader_arg(value):
    terms = value.split("=", 1)
    if len(terms) != 2:
        raise ValueError("reader arguments must be in the form NAME=VALUE")
    return terms

parser = argparse.ArgumentParser(
    description = "Compute molecular descriptors using RDKit")
parser.add_argument("--id-tag", default=None,
                    help="tag containing the identifier (default is the title)")

parser.add_argument("--format", choices=("smi", "sdf", "smi.gz", "sdf.gz"),
                    default="smi",
                    help="input file format, otherwise the format comes from "
                    "the structure file extension, or 'smi' if unknown")

parser.add_argument("--columns", "-P", default="id,smiles,MolWt",
                    help="comma-separated list of descriptors to compute (default: 'id,smiles,MolWt')")

parser.add_argument("--headers", "-H",
                    help="comma-separated list of headers (default: use the property names)")

parser.add_argument("--no-header", action="store_true",
                    help="Do not include the header in the output")

parser.add_argument("--output", "-o",
                    help="output filename (default: stdout)")

parser.add_argument("--dialect", choices=dialects, default="tab",
                    help="CSV output dialect, one of %s (default: 'tab')"
                    % (repr(dialects)[1:-1],))
                    
parser.add_argument("-R", dest="reader_args", type=parse_reader_arg, default=[], action="append",
                    help="add a reader argument")
                    
parser.add_argument("filename", nargs="?", default=None,
                    help="input structure file")


parser.add_argument("--list", action="store_true",
                    help="list the available property columns")

#### In the future I'll add the ability to specify the format for each field.
# This is a sketch to support things like: --fmt "{:2.3f},{:05d}".
# However, I think a per-column specifier like:
#   --fmt MW="{:2.3f}"
# might be better.


## _format_pat = re.compile(r"\{[^}]*\}")

## def split_format(fmt):
##     format_fields = []
##     i = 0
##     N = len(fmt)
##     while i < N:
##         m = _format_pat.match(fmt, i)
##         if m is None:
##             raise ValueError("Unsupported format starting with %r" % (fmt[i:],))
##         j = m.end()
##         format_fields.append(fmt[i:j])
##         if j == N:
##             break
##         if fmt[j:j+1] != ",":
##             raise ValueError("Expected comma at %r" % (fmt[j:],))
##         print("Found comma")
##         i = j+1
##     return format_fields

def main(argv=None):
    args = parser.parse_args(argv)

    resolver = rdkit_descriptors.resolver

    known_names = set(resolver.output_names)
    known_names.update(["id", "input_record", "input_format"])
    
    if args.list:
        for name in sorted(known_names, key=lambda s: s.lower()):
            sys.stdout.write(name + "\n")
        return

    columns = args.columns.split(",")
    if not columns:
        args.error("--columns my not be empty")

    headers = args.headers
    if headers is None:
        headers = columns
    else:
        headers = headers.split(",")
        if len(headers) != len(columns):
            args.error("--headers defines %d fields while --columns defines %d"
                       % (len(headers), len(columns)))
    
    for name in columns:
        if name not in known_names:
            parser.error("Unsupported column name %r" % (name,))
    

    format = T.get_input_format_from_source(args.filename, args.format)
    text_settings = dict(args.reader_args)
    reader_args = format.get_reader_args_from_text_settings(text_settings)
            
    properties = rdkit_descriptors.make_table_from_source(
        args.filename,
        format = format,
        id_tag = args.id_tag,
        reader_args = reader_args,
        )
    
    
    if args.output:
        properties.save(args.output, columns, headers=headers, dialect=args.dialect,
                        include_header = not args.no_header)
    else:
        properties.save(sys.stdout, columns, headers=headers, dialect=args.dialect,
                        include_header = not args.no_header)
                        
    
if __name__ == "__main__":
    main()
