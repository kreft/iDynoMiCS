#!/bin/bash

# Dowloads the Wiki, makes some changes to the files that are downloaded, and 
# then creates a pdf version of the tutorial
# Kieran Alden - 20th March 2013

# Firstly, use git to clone the wiki
# Make sure you're in a directory where you want the wiki to go
# Obviously git should be installed and configured correctly
git clone https://github.com/kreft/iDynoMiCS.wiki.git

# Now change into the iDynomics wiki directory
cd iDynoMiCS.wiki/

# We're going to remove the navigation on the bottom of each page - for better reading of PDF
# So search for all lines that begin 'Back'
for file in *.md
do
	# get the filename without an extension
	filename="${file%.*}"

	# Remove the 'Back to Part of the Github tutorial
	sed ' /\[Back/ c\ ' $filename.md>$filenameTemp.md
	# Remove the temp file and reinstate as the former file
	rm $filename.md
	mv $filenameTemp.md $filename.md

done

## Remove the menu from the intro page and replace (as these are github links)
sed '/### [0-9]./ c\ ' iDynomics-Tutorial.md>iDynomics-Tutorial2.md

sed '/## Tutorial Sections  ##/ a\\n### 1. Required Installations \n### 2. Running a Simulation in iDynoMiCS \n### 3. Protocol Files \n### 4. Result Files \n### 5. Matlab Routines to Analyse 2D Simulations \n### 6. Converting XML Output into Plain Columns of Numbers \n### 7. R Routines to Analyse 2D Simulations' iDynomics-Tutorial2.md>iDynomics-Tutorial.md

## Now for the sake of space and layout, we're going to join some pages.
# Firstly, Required software, Java, and Getting iDynomics
cat Required-\&-Optional-Software.md ../spacer.txt Getting-Java.md ../spacer.txt Tutorial\:-Getting-iDynoMiCS.md > Software_Java_iDyno.md

# Now ImageMagick and Matlab
cat ImageMagick.md Matlab.md > Matlab_ImageMagick.md

# Create the PDF for each file out of the wiki - this is done by converting markdown to pdf using gimli
# To install Gimli (on linux)
# sudo apt-get install aptitude
# sudo apt get install libxslt-dev libxml2-dev
# sudo aptitude install rubygems wkhtmltopdf ruby-1.9.1-dev
# sudo gem install gimli
gimli

# Now we're going to join all the PDFs in the order required for the tutorial
# First sudo apt-get install pdftk
pdftk iDynomics-Tutorial.pdf Software_Java_iDyno.pdf Matlab_ImageMagick.pdf PovRay.pdf Eclipse.pdf Running-a-Simulation.pdf RunIdyno.py.pdf iDynoMiCS-using-Eclipse.pdf Protocol-Files.pdf Defining-the-Simulation.pdf Starting-from-Previous-State-Files.pdf Defining-Solutes-and-Particle-Types.pdf Defining-the-World.pdf The-Bulk.pdf The-Computational-Domain.pdf Reactions.pdf Solute-Fields-and-Pressure-Solvers.pdf agentGrid.pdf Species.pdf Bacterium.pdf BactEPS.pdf BactAdaptable.pdf ParticulateEPS.pdf Result-Files.pdf agent_State-and-agent_Sum-Files.pdf env_State-and-env_Sum-Files.pdf POV-Ray-Result-Files.pdf Matlab-Routines.pdf Converting-XML-output.pdf R-Routines.pdf cat output tutorial.pdf

# Move the completed tutorial into its own directory
mkdir iDyno_Tutorial
mv tutorial.pdf iDyno_Tutorial/iDyno_Tutorial.pdf

# Now cleanup - delete the PDFs that gimli produced and the markdown files
rm *.pdf
rm *.md

cd ..
