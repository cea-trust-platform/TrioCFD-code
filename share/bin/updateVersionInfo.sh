#!/usr/bin/bash

# This script shall be used to update the trust version that must be used with
# the TrioCFD version of the current commit.
# This script shall be called by the integrator after each merge :
# git merge --no commit --no-ff ....
# ./updateVersion.sh
# git merge --continue

INFOFILE="./src/Version_info"

currentVersion=$(gawk '/TrioCFD v.* merge [0-9]+$/ {print $NF}' ${INFOFILE})
echo -e "Current version = ${currentVersion}"

if [[ ${currentVersion} == "" ]]
then
  echo -e "Error: the current version has not been found !"
  exit 1
fi

newVersion=$((${currentVersion} + 1))
echo -e "New version = ${newVersion}"

echo -e "The \${TRUST_ROOT} variable is ${TRUST_ROOT}"
if [[ ${TRUST_ROOT} == "" ]]
then
  echo -e "Error: the \${TRUST_ROOT} variable is not defined !"
  exit 1
fi

cd ${TRUST_ROOT}
shaID1=$(git log -n 1 --format=format:%H)
echo -e "The tust current sha ID1 is ${shaID1}"
cd -

echo -e "File ${INFOFILE} update"
sed -ri "s/(TrioCFD v.* merge) ${currentVersion}/\1 ${newVersion}/gmi" ${INFOFILE}
sed -ri "s/(Based on trust Id SHA1 :) .*/\1 ${shaID1}/gmi" ${INFOFILE}
git add -u ${INFOFILE}

exit 0
