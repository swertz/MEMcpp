##### Compute Matrix Element weights ####

How to install:
 * Fork this repo and/or clone it (username=swertz or AlexandreMertens): 
```
$ git clone https://github.com/(username)/MEMcpp
```
 
 * Setup environment (works on both ingrid-ui1 and lxplus). Starting from a clean environment, do:
```
$ cd MEMcpp
$ source init.sh
```

 * Build TTbar:
```
$ make ttbar -j8 
```

For now, the master branch allows to compute gg>tt~ weights in the e+/mu- channel, on parton-level events (no ISR correction). Usage:
```
$ examples/ME_ttbar_CUBA.exe /home/fynu/swertz/tests_MEM/MEMcpp/data/ttbar.root output.root /home/fynu/swertz/tests_MEM/binnedTF/TF_generator/Control_plots_hh_TF.root 0 0
```

Some comments:
* The ttbar.root file contains Delphes-parsed LHE evens.
* The transfer functions are binned transfer functions in electrons, muons and jets, built on a Delphes HH sample by Miguel.
* The two last arguments of the program call are start and end event numbers (0 0 computes the weight on the first event only)
* Sourcing init.sh will link to Sébastien's Delphes install. You can change your environment to link to your own install.
* Delphes is only used in main() to read the input datafile, and nowhere else.
