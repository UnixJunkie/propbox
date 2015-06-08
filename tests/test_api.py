import unittest

import propbox

class Constant(propbox.Resolver):
    output_names = ["n"]
    def __init__(self, n=0):
        self.n = n
    
    def resolve_column(self, name, table):
        table.set_values("n", [self.n] * len(table))

class Length(propbox.Resolver):
    output_names = ["len"]
    def resolve_column(self, name, table):
        smiles_list = table.get_values("smiles")
        table.set_values("len", map(len, smiles_list))

class TestMakeTableFromRecords(unittest.TestCase):
    def test_with_no_data(self):
        table = propbox.make_table_from_records(Constant(5), [])
        values = table.get_values("n")
        self.assertEqual(values, [])
        self.assertEqual(len(table), 0)
        
    def test_with_one_record(self):
        table = propbox.make_table_from_records(Constant(5), [{}])
        values = table.get_values("n")
        self.assertEqual(values, [5])
        self.assertEqual(len(table), 1)

    def test_with_two_records(self):
        table = propbox.make_table_from_records(
            Length(), [{"smiles": "C"}, {"smiles": "CDE"}])

        self.assertEqual(len(table), 2)
        values = table.get_values("len")
        self.assertEqual(values, [1, 3])

    def test_with_missing_field(self):
        with self.assertRaisesRegexp(ValueError,
                                     r"Record 1 has 0 elements \(\[\]\), expecting 1 \(\['smiles'\]\)"):
            propbox.make_table_from_records(
                Constant(11), [{"smiles": "C"}, {}])
        
    def test_with_two_missing_fields(self):
        with self.assertRaisesRegexp(ValueError,
            r"Record 2 has 1 elements \(\['MW']\), expecting 3 \(\['MW', 'id', 'smiles']\)"):
            propbox.make_table_from_records(
                Constant(11), [{"smiles": "C", "id": "spam", "MW": 1.3},
                               {"smiles": "D", "id": "eggs", "MW": 2.3},
                               {"MW": 3.4}])
            
    def test_with_extra_field(self):
        with self.assertRaisesRegexp(ValueError,
                                     r"Record 1 has 1 elements \(\['smiles'\]\), expecting 0 \(\[\]\)"):
            propbox.make_table_from_records(
                Constant(11), [{}, {"smiles": "C"}])
        
    def test_with_two_extra_fields(self):
        with self.assertRaisesRegexp(ValueError,
            r"Record 2 has 3 elements \(\['MW', 'id', 'smiles']\), expecting 1 \(\['smiles']\)"):
            propbox.make_table_from_records(
                Constant(11), [{"smiles": "C"},
                               {"smiles": "D"},
                               {"smiles": "E", "MW": 3.4, "id": "eggs"}])

    def test_with_changed_field(self):
        with self.assertRaisesRegexp(ValueError,
            r"Record 1 contains the unexpected descriptor 'frowns'"):
            propbox.make_table_from_records(
                Constant(11), [{"smiles": "C"},
                               {"frowns": "D"}])

    def test_auto_id_creation(self):
        table = propbox.make_table_from_records(None, [{}]*10)
        ids = table.get_values("id")
        self.assertEqual(ids, ["ID1", "ID2", "ID3", "ID4", "ID5",
                               "ID6", "ID7", "ID8", "ID9", "ID10"])
        
class TestMakeTableFromColumns(unittest.TestCase):
    def test_with_no_data(self):
        table = propbox.make_table_from_columns(Constant(5), {})
        values = table.get_values("n")
        self.assertEqual(values, [])
        self.assertEqual(len(table), 0)
        
    def test_with_one_record(self):
        table = propbox.make_table_from_columns(Constant(5), {"smiles": ["C"]})
        values = table.get_values("n")
        self.assertEqual(values, [5])
        self.assertEqual(len(table), 1)

    def test_with_two_records(self):
        table = propbox.make_table_from_columns(
            Length(), {"smiles": ["C", "CDE"]})

        self.assertEqual(len(table), 2)
        values = table.get_values("len")
        self.assertEqual(values, [1, 3])

    def test_with_missing_field(self):
        with self.assertRaisesRegexp(
                ValueError,
                r"Column length mismatch: 'MW' has length 1 while 'smiles' has length 2"):
            propbox.make_table_from_columns(
                Constant(11), {"smiles": ["C", "CC"],
                               "MW": [1.2]})
        
    def test_with_extra_field(self):
        with self.assertRaisesRegexp(
                ValueError,
                r"Column length mismatch: 'MW' has length 3 while 'smiles' has length 2"):
            propbox.make_table_from_columns(
                Constant(11), {"smiles": ["C", "CC"],
                               "MW": [1.2, 3.4, 5.6]})
        
    def test_auto_id_creation(self):
        table = propbox.make_table_from_columns(None, {"blah": [None]*10})
        ids = table.get_values("id")
        self.assertEqual(ids, ["ID1", "ID2", "ID3", "ID4", "ID5",
                               "ID6", "ID7", "ID8", "ID9", "ID10"])
    
    
if __name__ == "__main__":
    unittest.main()
