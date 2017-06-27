#!/bin/sh

#///////////////////////////////////////////////////////////////////////////////
# <file type="public">
#   <license>
#     This file was originally written by Bradley S. Meyer.
#
#     This is free software; you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.
#
#     This software is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     Please see the src/README.txt in this distribution for more copyright
#     and license information.
#   </license>
#   <description>
#     <abstract>
#       This routine will generate a regression test script for the
#       example codes in this directory.
#     </abstract>
#     <keywords>
#       regression, test, example, codes
#     </keywords>
#   </description>
# 
#   <authors>
#     <current>
#       <author userid="mbradle" start_date="2011/05/04" />
#     </current>
#   </authors>
# </file>
#///////////////////////////////////////////////////////////////////////////////

xsltproc http://wnmatrix.sf.net/xsl_pub/regression_test_script_maker.xsl examples.xml > regression_execute.sh
chmod +x regression_execute.sh
regression_execute.sh
