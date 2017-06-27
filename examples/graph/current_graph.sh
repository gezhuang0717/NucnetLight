#///////////////////////////////////////////////////////////////////////////////
#  This file was originally written by Michael J. Bojazi and Bradley S. Meyer.
# 
#  This is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
# 
#  This software is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this software; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#  USA
# 
#///////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////
#//!
#//! \file current_graph.sh
#//! \brief A short shell script to generate integrated current graphs.
#//!
#///////////////////////////////////////////////////////////////////////////////

#/bin/bash

#///////////////////////////////////////////////////////////////////////////////
# Check input.
#///////////////////////////////////////////////////////////////////////////////

if [ ! $# -eq 6 -a $1 == "--example" ]
  then
    echo -e "\n$0 ../network/my_currents.xml \"\" \"\" \"[z >= 20 and z <= 30 and a - z >= 20 and a - z <= 30]\" linear my_graph.pdf\n"
    exit
fi

if [ ! $# -eq 6 ]
  then
    echo -e "\nUsage: $0 input_xml nuc_xpath reac_xpath induced_nuc_xpath scaling out_pdf\n"
    echo -e "  input_xml = network xml file\n"
    echo -e "  nuc_xpath = xpath to select nuclides for valid reactions\n"
    echo -e "  reac_xpath = xpath to select reactions\n"
    echo -e "  induced_nuc_xpath = xpath to induce a subgraph of flows\n"
    echo -e "  scaling = scaling of flow weights (\"linear\" or \"logarithmic\")\n"
    echo -e "  out_pdf = output pdf file\n"
    echo -e "For an example, type: $0 --example\n"
    exit
fi

if [[ $6 != *.pdf ]]
then
  echo
  echo "name of output pdf file must end in .pdf"
  echo
  exit
fi

base_string=`echo ${6:0:${#6}-4}`

#///////////////////////////////////////////////////////////////////////////////
# Process files.
#///////////////////////////////////////////////////////////////////////////////

./zone_current_graph "$1" "$2" "$3" "$4" "$5" ${base_string}.dot

neato -n -Txdot ${base_string}.dot | dot2tex -tmath --crop > ${base_string}.tex

pdflatex ${base_string}

rm -rf *.log
rm -rf *.aux

#//////////////////////////////////////////////////////////////////////////////
# Remove or comment out these lines if you want to retain the individual dot
# and/or tex files.
#///////////////////////////////////////////////////////////////////////////////

rm -rf ${base_string}.tex
rm -rf ${base_string}.dot
