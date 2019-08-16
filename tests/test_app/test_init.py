"""testing the default import"""
from unittest import TestCase, main

from cogent3 import available_apps


class TestAvalableApps(TestCase):
    def test_available_apps(self):
        """available_apps returns a table"""
        from cogent3.util.table import Table

        apps = available_apps()
        self.assertIsInstance(apps, Table)
        self.assertTrue(apps.shape[0] > 10)


if __name__ == "__main__":
    main()
