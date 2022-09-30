foamListTimes -rm
if [ -d "postProcessing" ]; then
  echo "postProcessing does exist. Let's remove it"
  rm -r postProcessing
fi
hz_bioCycle > /dev/null
echo "Ended :smilyface:"