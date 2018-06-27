#!/bin/sh

#  rectangle.sh
#
#
#  Created by Eric Tovar on 2/26/2018.
#
# rough triangle areas 0.02, 0.005, 0.00125, 0.0003125, 0.000078125 using mesh sizes 0.2, 0.1, 0.05, 0.025, 0.0125
/Users/eric/Library/Mobile\ Documents/com~apple~CloudDocs/Texas\ A\&M/Math/spring17NumMthdsLab/triangle/triangle -npq28ea0.02 rectangle.poly
/Users/eric/Library/Mobile\ Documents/com~apple~CloudDocs/Texas\ A\&M/Math/spring17NumMthdsLab/triangle/triangle -rnpq28ea0.005  rectangle.1
/Users/eric/Library/Mobile\ Documents/com~apple~CloudDocs/Texas\ A\&M/Math/spring17NumMthdsLab/triangle/triangle -rnpq28ea0.00125  rectangle.2
/Users/eric/Library/Mobile\ Documents/com~apple~CloudDocs/Texas\ A\&M/Math/spring17NumMthdsLab/triangle/triangle -rnpq28ea0.0003125 rectangle.3
/Users/eric/Library/Mobile\ Documents/com~apple~CloudDocs/Texas\ A\&M/Math/spring17NumMthdsLab/triangle/triangle -rnpq28ea0.000078125 rectangle.4
#/Users/eric/Library/Mobile\ Documents/com~apple~CloudDocs/Texas\ A\&M/Math/spring17NumMthdsLab/triangle/triangle -rnpq28ea0.003125 rectangle.5
#/Users/eric/Library/Mobile\ Documents/com~apple~CloudDocs/Texas\ A\&M/Math/spring17NumMthdsLab/triangle/triangle -rnpq28ea0.0015625 rectangle.6
