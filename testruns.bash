

echo Do first run
lamboot
make clean && make && make runmpi
lamhalt

echo analyse
cp params.f90 output
mkdir analysis
python start_python.py output/ analysis/

echo save results
mv output runone_sph
mv analysis runone_sph_analysis
mkdir output

echo create next run
sed 's/osoft=0.0/osoft=0.1309/g' params.f90 > params.f90

echo do second run
lamboot
make clean && make && make runmpi
lamhalt

echo analyse
cp params.f90 output
mkdir analysis
python start_python.py output/ analysis/

echo save results
mv output runtwo_col
mv analysis runtwo_col_analysis
mkdir output
