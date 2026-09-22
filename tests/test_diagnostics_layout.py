"""keep diagnostic assets and entrypoints independent from example working directories"""

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class diagnostics_layout_tests(unittest.TestCase):
    def test_diagnostics_are_separate_from_examples(self):
        root = Path(__file__).resolve().parents[1]
        directory = root / "diagnostics"
        self.assertTrue((directory / "README.md").is_file())
        self.assertTrue((directory / "external_solver_inputs").is_dir())
        for name in ("diagnostics", "external_solver_inputs", "external_solver_compare.py", "run_comparison_tests.sh"):
            self.assertFalse((root / "examples" / "gravmag_sphere" / name).exists())

    def test_formatted_case_artifact_paths_exist(self):
        root = Path(__file__).resolve().parents[1]
        metadata = list((root / "diagnostics" / "external_solver_inputs").glob("*/metadata.json"))
        self.assertTrue(metadata)
        for path in metadata:
            values = json.loads(path.read_text())
            # metadata records the local artifact location; resolve its repository suffix portably
            for key in ("example_path", "baseline_xyz", "baseline_brtp"):
                recorded = values[key]
                relative = recorded.split("/GravMagPy/", 1)[-1]
                self.assertTrue((root / relative).is_file(), f"{path.name}: {key}")

    def test_entrypoint_help_outside_repository(self):
        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as directory:
            environment = dict(os.environ, MPLCONFIGDIR=directory)
            for script in sorted((root / "diagnostics").glob("*.py")):
                result = subprocess.run([sys.executable, str(script), "--help"], cwd=directory,
                                        env=environment, capture_output=True, text=True, timeout=30)
                self.assertEqual(result.returncode, 0, f"{script.name}: {result.stderr}")
            result = subprocess.run(["bash", str(root / "diagnostics" / "run_comparison_tests.sh"), "--help"],
                                    cwd=directory, capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)
