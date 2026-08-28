import unittest

from build_tcell_library_v2 import (
    canonical_tiles,
    is_allowed_non_cns_control_source,
    normalize_iedb_tcell_rows,
    parse_epitope_coordinates,
)


class LiveEvidenceNormalizationTests(unittest.TestCase):
    def test_epitope_summary_coordinates_are_one_based_and_exact(self):
        self.assertEqual(parse_epitope_coordinates("PEPTIDE<br/>Protein name (627-641)<br/>EBV"), (627, 641))
        self.assertEqual(parse_epitope_coordinates("No coordinates"), (None, None))

    def test_iedb_records_collapse_duplicate_assays_not_distinct_peptides(self):
        base = {
            "linear_sequence": "TGGVYHFVKKHVHES",
            "linear_sequence_length": 15,
            "epitope_structure_defined": "Exact Epitope",
            "epitope_summary": "TGGVYHFVKKHVHES<br/>DNA polymerase (627-641)<br/>EBV",
            "parent_source_antigen_iri": "UNIPROT:P03198",
            "parent_source_antigen_name": "DNA polymerase (UniProt:P03198)",
            "source_organism_name": "Human herpesvirus 4 (Epstein Barr virus)",
            "host_organism_name": "Homo sapiens (human)",
            "qualitative_measure": "Positive",
            "mhc_class": "II",
            "mhc_allele_name": "HLA-DRB5*01:01",
            "assay_names": "proliferation",
            "pubmed_id": "12244309",
            "tcell_id": 10,
        }
        duplicate = dict(base, assay_names="cytokine release", tcell_id=11)
        other = dict(base, linear_sequence="GGVYHFVKKHVHES", linear_sequence_length=14)
        rows = normalize_iedb_tcell_rows([base, duplicate, other], "EBV")
        self.assertEqual(len(rows), 2)
        full = next(row for row in rows if row["sequence"] == "TGGVYHFVKKHVHES")
        self.assertEqual(full["source_start_1_based"], 627)
        self.assertEqual(full["source_end_1_based"], 641)
        self.assertEqual(full["supporting_record_count"], 2)
        self.assertEqual(full["source_accession"], "P03198")

    def test_missing_and_present_coordinates_can_be_sorted_together(self):
        base = {
            "linear_sequence": "TGGVYHFVKKHVHES",
            "linear_sequence_length": 15,
            "epitope_structure_defined": "Exact Epitope",
            "epitope_summary": "TGGVYHFVKKHVHES<br/>DNA polymerase (627-641)<br/>EBV",
            "parent_source_antigen_iri": "UNIPROT:P03198",
            "parent_source_antigen_name": "DNA polymerase (UniProt:P03198)",
            "source_organism_name": "Human herpesvirus 4 (Epstein Barr virus)",
            "host_organism_name": "Homo sapiens (human)",
            "qualitative_measure": "Positive",
            "mhc_class": "II",
            "mhc_allele_name": "HLA class II",
            "assay_names": "proliferation",
            "pubmed_id": "12244309",
        }
        missing = dict(base, linear_sequence="GGVYHFVKKHVHES", epitope_summary="No source coordinates")
        rows = normalize_iedb_tcell_rows([base, missing], "EBV")
        self.assertEqual(len(rows), 2)
        missing_row = next(row for row in rows if row["sequence"] == "GGVYHFVKKHVHES")
        self.assertEqual(missing_row["source_start_1_based"], "")
        self.assertEqual(missing_row["source_end_1_based"], "")

    def test_canonical_tiles_are_reproducible_and_inside_declared_region(self):
        sequence = "ACDEFGHIKLMNPQRSTVWY" * 10
        tiles = canonical_tiles(
            protein_symbol="ANO2",
            accession="Q9NQ90",
            sequence=sequence,
            region_start=79,
            region_end=168,
            count=4,
        )
        self.assertEqual(len(tiles), 4)
        self.assertEqual(tiles, canonical_tiles(
            protein_symbol="ANO2",
            accession="Q9NQ90",
            sequence=sequence,
            region_start=79,
            region_end=168,
            count=4,
        ))
        self.assertTrue(all(79 <= row["start"] <= row["end"] <= 168 for row in tiles))
        self.assertTrue(all(len(row["sequence"]) == 15 for row in tiles))

    def test_non_cns_control_filter_excludes_every_study_protein_by_accession(self):
        for accession in ("P02686.3", "P60201", "Q16653-2", "Q9NQ90", "P09543", "Q13875", "P20916", "Q02246", "O75508", "P37837", "P02511.2"):
            self.assertFalse(is_allowed_non_cns_control_source(accession, "opaque historical name"))
        self.assertFalse(is_allowed_non_cns_control_source("UNKNOWN", "MBP protein"))
        self.assertTrue(is_allowed_non_cns_control_source("Q05329.1", "Glutamate decarboxylase 2"))


if __name__ == "__main__":
    unittest.main()
