# Positive Resampler for Negative Weights: CMSSW plugin

This is the github repo for a positive resampler within CMSSW framework. To use it, execute shell commands as showing below

```shell
export SCRAM_ARCH="amd64_slcx_gccxxx" # Please note the unique_ptr/auto_ptr issue for different SCRAM_ARCH
cmsrel CMSSW_X_Y_Z # Developed with CMSSW_11_2_0
cd CMSSW_X_Y_Z/src
cmsenv
# add clone this repo~
```