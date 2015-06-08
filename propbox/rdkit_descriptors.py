from __future__ import print_function, division, absolute_import

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import SaltRemover

from . import pylru
from . import collect_resolvers, calculate, Propbox, make_table_from_columns, CalculateName
from . import simple_futures
import os

from collections import namedtuple

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

####

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

#### Some ways to get the SMILES

@calculate()
def calc_smiles(mol):
    return Chem.MolToSmiles(mol, True) # use isomeric canonical SMILES
    
@calculate()
def calc_cansmiles(mol):
    return Chem.MolToSmiles(mol, False) # use non-isomeric canonical SMILES
    
@calculate()
def calc_usmsmiles(mol):
    return Chem.MolToSmiles(mol, False, canonical=False) # use non-isomeric, non-canonical SMILES



### Descriptors from the Descriptors module

# Molecular weight

@calculate()
def calc_MolWt(mol):
    "The average molecular weight of the molecule"
    return Descriptors.MolWt(mol)

@calculate()
def calc_HeavyAtomMolWt(mol):
    "The average molecular weight of the non-hydrogen atoms in the molecule"
    return Descriptors.HeavyAtomMolWt(mol)

@calculate()
def calc_ExactMolWt(mol):
    "The exact molecular weight of the molecule, taking isotopes into account"
    return Descriptors.ExactMolWt(mol)

@calculate()
def calc_MolWt_version():
    return Descriptors.MolWt.version


# Electrons

@calculate()
def calc_NumValenceElectrons(mol):
    "The number of valence electrons in the molecule"
    return Descriptors.NumValenceElectrons(mol)

@calculate()
def calc_NumRadicalElectrons(mol):
    "The number of radical electrons in the molecule (says nothing about spin state)"
    return Descriptors.NumRadicalElectrons(mol)


## Charge descriptors

_charge_version = Descriptors.MaxPartialCharge.version

ChargeDescriptor = namedtuple("ChargeDescriptor", "minCharge maxCharge")

# An example of using propbox to manage intermediate values

@calculate(output_names=["_chargeDescriptors", "chargeDescriptorVersion"])
def calc__chargeDescriptors(mol):
    "A helper function to compute the most negative and most postive Gasteiger charges"
    Descriptors._ChargeDescriptors(mol)
    minCharge, maxCharge = mol._chargeDescriptors
    del mol._chargeDescriptors
    return (ChargeDescriptor(minCharge, maxCharge)), _charge_version

@calculate()
def calc_MaxPartialCharge(_chargeDescriptors):
    "The most positive Gasteiger partial charge"
    return _chargeDescriptors.maxCharge

@calculate()
def calc_MinPartialCharge(_chargeDescriptors):
    "The most negative Gasteiger partial charge"
    return _chargeDescriptors.minCharge

@calculate()
def calc_MaxAbsPartialCharge(_chargeDescriptors):
    "The largest of the absolute value of the minimum and maximum Gasteiger partial charges"
    return max(map(abs, _chargeDescriptors))

@calculate()
def calc_MinAbsPartialCharge(_chargeDescriptors):
    "The smallest of the absolute value of the minimum and maximum Gasteiger partial charges"
    return min(map(abs, _chargeDescriptors))


## 

def load_descriptor_module():
    rdkit_descriptors = Propbox()
    seen = [k[5:] for k in globals() if k.startswith("calc_")]
    
    for k, v in Descriptors.__dict__.items():
        if k in seen:
            continue
        if k[:1] == "_":
            continue
        if not hasattr(v, "version"):
            continue

        rdkit_descriptors.add_resolver(CalculateName(["mol"], k, v, v.__doc__))

    return rdkit_descriptors

descriptor_module = load_descriptor_module()



######

    
_resolver = resolver = collect_resolvers()  # get all of the calculators

# Break the resolver into independent stages.

_input_processing = set(["calc_mol", "calc_input_mol"])
input_resolver = Propbox(r for r in resolver if r.resolver_name in _input_processing)
mol_resolver = Propbox(r for r in resolver if r.resolver_name not in _input_processing)


def make_table_from_source(source, format=None, id_tag=None, reader_args={},
                           resolver=None, config=None):
    from . import rdkit_toolkit as T
    reader = T.read_ids_and_molecules(source, format, id_tag=id_tag,
                                      reader_args=reader_args)
    return make_table_from_ids_and_molecules(reader, resolver, config)

def make_table_from_ids_and_molecules(reader, resolver=None, config=None):
    if resolver is None:
        resolver = _resolver
        
    if config is None:
        config = {}
    
    has_smiles = None
    input_smiles = None

    id_list = []
    mol_list = []
    smiles_list = []

    new_future = simple_futures.new_future
    new_future_exception = simple_futures.new_future_exception
    
    for mol_index, (id, mol) in enumerate(reader, 1):
        if has_smiles is None:
            if mol is None:
                input_smiles = None
            else:
                has_smiles = mol.HasProp("_input_smiles")
                if has_smiles:
                    input_smiles = mol.GetProp("_input_smiles")
                    mol.ClearProp("_input_smiles")
        elif has_smiles:
            if mol is None:
                input_smiles = None
            else:
                input_smiles = mol.GetProp("_input_smiles")
                mol.ClearProp("_input_smiles")

        if id is None:
            id = "ID%d" % (mol_index,)
        id_list.append(new_future(id))
        
        if mol is None:
            mol_list.append(new_future_exception(ValueError("Cannot parse record %d with id %r"
                                                            % (mol_index, id))))
        else:
            mol_list.append(new_future(mol))
            
        smiles_list.append(new_future(input_smiles))

    if has_smiles:
        # Perhaps this should be a 'config' rather than a property?
        input_format = ["smi"] * len(smiles_list)
    else:
        input_format = [None] * len(smiles_list)
    columns = {
        "id": id_list,
        "input_mol": mol_list,
        "input_record": smiles_list,
        "input_format": input_format,
        }
    return make_table_from_columns(resolver, columns, config)
