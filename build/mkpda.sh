#run from build directory!
#start from scratch
mv pda deleteme
rm -fR pda-cran
#rsync and exclude directories that CRAN doesn't like
mkdir pda-cran
rsync -av ../R ./pda-cran/
rsync -av ../man ./pda-cran/
rsync -av ../src ./pda-cran/ --exclude=*.sw* --exclude=*.so --exclude=*.o --exclude=*.rds 
cp ../DESCRIPTION ./pda-cran/
cp ../NAMESPACE  ./pda-cran
#create the source package
R CMD build pda-cran
#extract
tar -xzvf pda_1.0.tar.gz
#check the package
R CMD check --as-cran pda
#build the binary
R CMD INSTALL --build pda_1.0.tar.gz
