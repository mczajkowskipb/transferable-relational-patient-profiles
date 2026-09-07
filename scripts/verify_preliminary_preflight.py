#!/usr/bin/env python3
"""Run the existing zero-fixture tests and new unittest integration tests.

This dependency-light runner is limited to the current test layout. It fails
explicitly if a pytest fixture or parametrization is added; CI still uses pytest.
It generates no real-data evidence and does not freeze/unseal any protocol.
"""
from pathlib import Path
import importlib.util
import inspect
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
suite = unittest.TestSuite()
for path in sorted((ROOT / 'tests').glob('test_*.py')):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    suite.addTests(unittest.defaultTestLoader.loadTestsFromModule(module))
    for name, fn in inspect.getmembers(module, inspect.isfunction):
        if name.startswith('test_') and fn.__module__ == module.__name__:
            if inspect.signature(fn).parameters or getattr(fn, 'pytestmark', None):
                raise RuntimeError(f'{path.name}:{name} requires the pytest runner')
            suite.addTest(unittest.FunctionTestCase(fn, description=f'{path.name}:{name}'))
result = unittest.TextTestRunner(verbosity=2).run(suite)
sys.exit(0 if result.wasSuccessful() else 1)
