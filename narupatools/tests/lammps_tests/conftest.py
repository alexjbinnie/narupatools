# This file is part of narupatools (https://github.com/alexjbinnie/narupatools).
# Copyright (c) Alex Jamieson-Binnie. All rights reserved.
#
# narupatools is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# narupatools is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with narupatools.  If not, see <http://www.gnu.org/licenses/>.


import os

import pytest
from testing import add_mark


def pytest_collection_modifyitems(items):
    add_mark(filename=__file__, mark=pytest.mark.lammps, items=items)


# Change the working directory so LAMMPS can find input files.
@pytest.fixture(scope="module", autouse=True)
def _change_test_dir(request):
    os.chdir(request.fspath.dirname)
    yield
    os.chdir(request.config.invocation_dir)
