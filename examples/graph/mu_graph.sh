#///////////////////////////////////////////////////////////////////////////////
#  Copyright (c) 2012-2013 Clemson University.
# 
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
#//! \file flow_graph.sh
#//! \brief A short shell script to generate mu graphs.
#//!
#///////////////////////////////////////////////////////////////////////////////

#/bin/bash

#///////////////////////////////////////////////////////////////////////////////
# Edit these parameters.
#///////////////////////////////////////////////////////////////////////////////

FILE_BASE=graph.dot
INPUT_DIR=dot_files
OUTPUT_DIR=pdf_files

#///////////////////////////////////////////////////////////////////////////////
# End of edit.
#///////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////
# Check input.
#///////////////////////////////////////////////////////////////////////////////

if [ ! $# -eq 6 -a $1 == "--example" ]
  then
    echo -e "\n$0 ../network/my_output.xml \"[z >= 82 and z <= 84 and a - z >= 124 and a - z <= 128]\" \"[position() >= last() - 10]\" 1. 1.e-25 my_graph.pdf\n"
    exit
fi

if [ ! $# -eq 6 ]
  then
    echo -e "\nUsage: $0 input_xml nuc_xpath zone_xpath cutoff out_pdf\n"
    echo -e "  input_xml = network xml file\n"
    echo -e "  nuc_xpath = xpath to select nuclides\n"
    echo -e "  zone_xpath = xpath to select zones\n"
    echo -e "  scale = scaling for chemical potentials\n"
    echo -e "  cutoff = minimum abundance\n"
    echo -e "  out_pdf = output pdf file\n"
    echo -e "For an example, type: $0 --example\n"
    exit
fi

#!/bin/bash

ORIGINAL_DIR=`pwd`

if [ -d "$INPUT_DIR" ]; then
  rm -fr $INPUT_DIR
fi
if [ -d "$OUTPUT_DIR" ]; then
  rm -fr $OUTPUT_DIR
fi

mkdir $INPUT_DIR 
mkdir $OUTPUT_DIR

./zone_mu_graph "$1" "$2" "$3" "$4" "$5" $INPUT_DIR/$FILE_BASE

#///////////////////////////////////////////////////////////////////////////////
# Find all the input dot files and put them into an array.
#///////////////////////////////////////////////////////////////////////////////

inputFilesArray=($(find $INPUT_DIR -type f))

#///////////////////////////////////////////////////////////////////////////////
# Now, find the length of the array.
#///////////////////////////////////////////////////////////////////////////////

numfiles=${#inputFilesArray[@]}

#///////////////////////////////////////////////////////////////////////////////
# Loop over all the files in the array.
#///////////////////////////////////////////////////////////////////////////////

for(( i=0; i<numfiles; i++ )); 
do  
  neato -n -Txdot "${inputFilesArray[$i]}" | dot2tex -tmath --crop > $OUTPUT_DIR/$(basename "${inputFilesArray[$i]}").tex 
  pdflatex -output-directory $OUTPUT_DIR $OUTPUT_DIR/$(basename "${inputFilesArray[$i]}").tex 
  rm $OUTPUT_DIR/*.log $OUTPUT_DIR/*.tex $OUTPUT_DIR/*.aux; 
done  

outputFilesArray=("$(ls $OUTPUT_DIR | sort -n -t _ -k 5)")

cd $OUTPUT_DIR

#///////////////////////////////////////////////////////////////////////////////
# Use ghostscript to combine all the separate pdf files into a single pdf file.
#///////////////////////////////////////////////////////////////////////////////

gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=${ORIGINAL_DIR}/$6 ${outputFilesArray[@]}

cd "$ORIGINAL_DIR"

#///////////////////////////////////////////////////////////////////////////////
# Remove these lines if you want to retain the individual dot and tex files.
#///////////////////////////////////////////////////////////////////////////////

rm -rf $INPUT_DIR 
rm -fr $OUTPUT_DIR
