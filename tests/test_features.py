import contextlib
import io
import unittest

from moliety2.features import FEATURE_MODES, depict_input_smiles, process_smiles_main
from moliety2.patterns import daylight_smarts_patterns, smartsrx_patterns


class ProcessSmilesMainTests(unittest.TestCase):
    def test_input_depiction_is_separate_from_feature_outputs(self):
        input_image = depict_input_smiles("CCO")
        output_images, _message = process_smiles_main("CCO", "Functional Groups")

        self.assertIsInstance(input_image, str)
        self.assertNotIn(input_image, [image for image, _caption in output_images])

    def test_invalid_smiles_returns_empty_gallery(self):
        images, message = process_smiles_main("not-a-smiles", "Functional Groups")

        self.assertEqual(images, [])
        self.assertEqual(message, "Invalid SMILES.")

    def test_invalid_smiles_returns_no_input_depiction(self):
        self.assertIsNone(depict_input_smiles("not-a-smiles"))

    def test_invalid_mode_returns_empty_gallery(self):
        images, message = process_smiles_main("CCO", "Unknown Mode")

        self.assertEqual(images, [])
        self.assertEqual(message, "Invalid mode selected.")

    def test_all_modes_return_lists(self):
        smiles = "CC(=O)Oc1ccccc1C(=O)O"

        for mode in FEATURE_MODES:
            with self.subTest(mode=mode.name):
                images, message = process_smiles_main(smiles, mode.name)

                self.assertIsInstance(images, list)
                self.assertIsInstance(message, str)

    def test_no_result_modes_return_empty_lists_not_none(self):
        for mode_name in ("Chiral Centers", "Potential Stereogenic Centers"):
            with self.subTest(mode=mode_name):
                images, message = process_smiles_main("CCO", mode_name)

                self.assertEqual(images, [])
                self.assertIn("No ", message)

    def test_protonation_rejects_inverted_ph_range(self):
        images, message = process_smiles_main("CCO", "Protonation", 8.0, 7.0)

        self.assertEqual(images, [])
        self.assertEqual(message, "Minimum pH must be less than or equal to maximum pH.")


class PatternLoadingTests(unittest.TestCase):
    def test_daylight_smarts_loading_skips_invalid_patterns_quietly(self):
        daylight_smarts_patterns.cache_clear()

        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            patterns = daylight_smarts_patterns()

        self.assertGreater(len(patterns), 0)
        self.assertNotIn("SMARTS Parse Error", stderr.getvalue())

    def test_smartsrx_loading_skips_invalid_patterns_quietly(self):
        smartsrx_patterns.cache_clear()

        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            patterns = smartsrx_patterns()

        self.assertGreater(len(patterns), 0)
        self.assertIn("Acid > Acid_Aliphatic > Acid_SaturatedAliphatic", patterns)
        self.assertNotIn("SMARTS Parse Error", stderr.getvalue())

    def test_smartsrx_mode_finds_carboxylic_acid(self):
        images, message = process_smiles_main("CC(=O)O", "SMARTS-RX Moieties")

        self.assertGreater(len(images), 0)
        self.assertIn("Found ", message)


if __name__ == "__main__":
    unittest.main()
