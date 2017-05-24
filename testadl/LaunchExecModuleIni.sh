#!/bin/bash

FILE=$1

SWMOD_PATH=/path/to/swmod

. ${SWMOD_PATH}/bin/swmod.sh init

swmod load root@6.05.02
swmod load clhep@2.1.3.1
swmod load gerda@master

echo "execModuleIni -o ${FILE%-*}-Tier2.root $2/configfiles/inifile.ini ${FILE}"
execModuleIni -o $2/${FILE%-*}-Tier2.root $2/configfiles/inifile.ini $2/${FILE}
