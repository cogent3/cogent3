from unittest import TestCase, main

from numpy.testing import assert_allclose

from cogent3.draw import drawable


__author__ = "Gavin Huttley and Rahul Ghangas"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley", "Rahul Ghangas"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Alpha"


class Test_Shapes(TestCase):
    def test_shift(self):
        """Tests if the shape remains consistent with a shift of coordinate system along the line x=y"""
        rectangle = drawable.Rectangle([(5, 6)], 7, 0.25)
        rectangle.shift(5, 5)

        assert_allclose(rectangle.x, [10, 10, 11, 11, 10])
        assert_allclose(rectangle.y, [12, 12.25, 12.25, 12, 12])

    def test_linear_transform_shift(self):
        """Tests that the side lengths remain consistent with any arbitrary shift in the coordinate system"""
        rectangle = drawable.Rectangle([(10, 15)], 17, 3)
        rectangle.shift(2, 6)

        assert_allclose(rectangle.x, [12, 12, 17, 17, 12])
        assert_allclose(rectangle.y, [23, 26, 26, 23, 23])

    def test_height_property(self):
        """Tests if the height remains consistent after any arbitrary shift in the coordinate system"""
        rectangle = drawable.Rectangle([(5, 11)], 6, 7)
        rectangle.shift(5, 2)

        self.assertTrue(rectangle.height == 7, "Height is consistent with shift")


if __name__ == "__main__":
    main()
