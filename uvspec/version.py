# Generate UV-Vis spectra from Gaussian09 TDHF/TDDFT log files. 
# Copyright (C) 2013  Gaussian Toolkit
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Return the program version number.

The version number is stored as a 4-item tuple: (MAJOR, MINOR, PATCH, BUILD).
If BUILD is greater than zero, the current version is a beta release.

"""
from uvspec import VERSION as version


def get_version():
    """Return the version number as MAJOR.MINOR.PATCH.BUILD"""
    main = '.'.join(str(x) for x in version[:2])
    sub = 'b' + str(version[3]) if version[3] > 0 else ''
    return str(main + sub) 
