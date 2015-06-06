from __future__ import print_function, division

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import SaltRemover

import pylru
from propbox2 import collect_resolvers, calculate, Propbox

@calculate()
def calc_input_mol(input_record, input_format):
    """Convert the input structure into an RDKit molecule

    Convert the `input_record` string, which is a structure
    record in format `input_format`, into an RDMol. The format
    must be one of 'smi', 'smiles', 'smistring', or 'sdf'.
    
    :param input_record: the input structure
    :type input_record: a string
    :param input_format: one of 'smiles' or 'sdf'
    :returns: `input_mol` as the RDKit molecule
    """
    if input_format in ("smi", "smiles", "smistring"):
        mol = Chem.MolFromSmiles(input_record)
    elif input_format in ("sdf",):
        sup = Chem.SDMolSupplier()
        sup.SetData(input_record)
        try:
            mol = next(sup)
        except StopIterationError:
            mol = None
    else:
        raise ValueError("Unsupported format %r" % (input_format,))
    if mol is None:
        raise ValueError("Could not parse record in %r format" % (input_format,))
    return mol

_get_salt_remover = pylru.FunctionCacheManager(
    SaltRemover.SaltRemover, 50)

_SALT_REMOVER_KEY = "1d3df49d-1a8e-4039-8cb1-7ecb38ba0ea1"
@calculate(include_table=True)
def calc_mol(input_mol, table):
    """Clean up the input_mol so it can be used for the other descriptors

    This removes the salts from the `input_mol` and places the resulting
    molecule in `mol`.

    Use the configuration option `SaltRemover` as the salt remover.
    If it does not exist or is None, use `SaltRemover_defnFilename`
    or `SaltRemover_defnData` to specify the salts to remove. Use the
    configuration option `SaltRemover_dontRemoveEverything` to specify
    if salt removal can leave the molecule empty (default is False).
    
    """
    data = table.get_cache_value(_SALT_REMOVER_KEY)
    if data is None:
        config = table.config
        remover = config.get("SaltRemover", None)
        if remover is None:
            filename = config.get("SaltRemover_defnFilename", None)
            data = config.get("SaltRemover_defnData", None)
            remover = _get_salt_remover(filename, data)
        dontRemoveEverything = config.get("SaltRemover_dontRemoveEverything", False)
        table.set_cache_value(_SALT_REMOVER_KEY, (remover, dontRemoveEverything))
    else:
        remover, dontRemoveEverything = data
    
    return remover.StripMol(input_mol, dontRemoveEverything)


@calculate()
def calc_MolWt(mol):
    return Descriptors.MolWt(mol)


resolver = collect_resolvers()  # get all of the calculators

# Break the resolver into indepenent stages.

_input_processing = set(["calc_mol", "calc_input_mol"])
input_resolver = Propbox(r for r in resolver if r.resolver_name in _input_processing)
mol_resolver = Propbox(r for r in resolver if r.resolver_name not in _input_processing)

if __name__ == "__main__":
    import propbox2
    smiles_list = ["c1ccccc1O", "CCO", "Q"]
    
    properties = propbox2.make_table_from_columns(
        resolver,
        {"input_record": smiles_list,
         "input_format": ["smiles"]*len(smiles_list),
         })
    
    print(properties.get_column_values("MolWt"))
    
    properties.save("/dev/stdout", ["id", "MolWt"], formatters=[str, "{:.2f}".format], dialect="excel2")
