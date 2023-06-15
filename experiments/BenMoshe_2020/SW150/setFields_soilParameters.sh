cd constant/
rm -r soilParameters
tar -xvf soilParameters.tar.gz 
mv soilParameters/* ../0.000/
setFields
cd 0.000
mv K_0 alpha n_vangenucthen Sw_s Sw_r ../constant/soilParameters
cd ..