#!/bin/bash

xyz_files=$1
chem_networks=$2
output_file=$3
chem_input=$4
chem_binary=$5

counter=0
back_to_here=`pwd`

echo 'XYZ file path : '$xyz_files
echo 'ChemNetWorks binary file path : '$chem_networks
echo 'PWD : '$back_to_here

create_clean_directory(){
    dir_name=$1
    if [ -d "$dir_name" ]; then
        echo "Removing $dir_name"
        rm -rf "$dir_name"
    elif [ -f "$dir_name" ]; then
        echo "File with this name already exists, not a directory."
        exit
    fi
    if mkdir "$dir_name"; then
        echo "Clean directory created: $dir_name"
        return 0
    else
        echo "Creating directory failed: $dir_name"
        return 1
    fi
}

for xyz_file in $(ls $xyz_files* | sort)
do
    if (($counter%10==0))
    then
        folder=$output_file$counter
        create_clean_directory $folder
        cp $chem_networks $folder
        cp $xyz_file $folder
        cp $chem_input $folder
        cd $folder
        ./$chem_binary `(basename $chem_input)` `(basename $xyz_file)`
        cd $back_to_here
    fi
    counter=`echo "$counter+1" | bc`
    echo 'Counter :' $counter
done
