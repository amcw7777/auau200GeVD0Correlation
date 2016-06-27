 #!/usr/bin/env csh
# CHANGE the following to suit your environment
setenv BASEDIR `pwd`
 
### ROOT
#setenv ROOTSYS ${BASEDIR}/root/
 
### FastJet
setenv FASTJETDIR ${BASEDIR}/fastjet3
 
# ### PYTHIA8
# setenv PYTHIA8DIR ${BASEDIR}/pythia8
# setenv PYTHIA8DATA ${PYTHIA8DIR}/xmldoc
 
### TStarJetPicoDst structure
setenv STARPICOPATH ${BASEDIR}/eventStructuredAu
 
 
############## Done with indivivual settings.
 
###### Update paths
if(! $?LD_LIBRARY_PATH) setenv LD_LIBRARY_PATH
if(! $?DYLD_LIBRARY_PATH) setenv DYLD_LIBRARY_PATH
 
setenv PATH ${FASTJETDIR}/bin:./bin:${PATH}
setenv LD_LIBRARY_PATH ${FASTJETDIR}/lib:${STARPICOPATH}:${LD_LIBRARY_PATH}                                                                                       
#setenv DYLD_LIBRARY_PATH ${FASTJETDIR}/lib:${STARPICOPATH}:${PYTHIA8DIR}/lib:${DYLD_LIBRARY_PATH}
 
if ($?TERM == 0 || $?prompt == 0) exit 0
 
echo ''
echo 'Setup FastJet and picoDST libraries'
echo '====================================='
echo ''
echo "<I>---------------Info--------------------<I>"
echo "Setting up the following environments: "
echo "ROOT: " $ROOTSYS
#echo "PYTHIA8: " $PYTHIA8DIR
echo "FastJet: " $FASTJETDIR
echo "STARPICOPATH: " $STARPICOPATH
 
echo "<I>---------------Info--------------------<I>"
echo ""
